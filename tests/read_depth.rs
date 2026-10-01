//! `mung infercnv` on a read-depth matrix whose rows are genomic bins
//! (`chr:start-end`, as `faba depth` writes them): no annotation is needed,
//! bins keep their names, and a gain shows as a positive log-ratio.

use data_beans::sparse_io::SparseIoBackend;
use data_beans::sparse_io::{create_sparse_from_triplets, open_sparse_matrix_by_path};
use std::path::PathBuf;
use std::process::{Command, Output};

const BIN: i64 = 1_000_000;

/// A temporary directory, removed when dropped (also when a test fails).
struct TmpDir(PathBuf);

impl TmpDir {
    fn new(name: &str) -> Self {
        let d = std::env::temp_dir().join(format!("mung-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        Self(d)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_string()
    }
}

impl Drop for TmpDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// `n` bins of [`BIN`] bp on `chr`, starting at `offset`.
fn bins(chr: &str, n: i64, offset: i64) -> Vec<Box<str>> {
    (0..n)
        .map(|i| format!("{chr}:{}-{}", offset + i * BIN, offset + (i + 1) * BIN).into())
        .collect()
}

/// A backend of `n_cells` cells over `rows`; rows on `gain` get twice the
/// depth.
fn write_backend(path: &str, rows: &[Box<str>], prefix: &str, n_cells: usize, gain: Option<&str>) {
    let mut triplets = Vec::new();
    for j in 0..n_cells {
        for (i, name) in rows.iter().enumerate() {
            let base = 20.0 + ((i * 7 + j * 3) % 5) as f32;
            let doubled = gain.is_some_and(|g| name.starts_with(&format!("{g}:")));
            triplets.push((i as u64, j as u64, if doubled { 2.0 * base } else { base }));
        }
    }
    let shape = (rows.len(), n_cells, triplets.len());
    let mut data =
        create_sparse_from_triplets(&triplets, shape, Some(path), Some(&SparseIoBackend::Zarr))
            .unwrap();
    data.register_row_names_vec(rows);
    let cells: Vec<Box<str>> = (0..n_cells)
        .map(|j| format!("{prefix}{j}@D1").into())
        .collect();
    data.register_column_names_vec(&cells);
}

/// `mung` offline, with an empty cache and no user config.
fn mung(dir: &TmpDir, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_mung"))
        .args(args)
        .env("MUNG_OFFLINE", "1")
        .env("MUNG_CACHE_DIR", dir.path("cache"))
        .env("MUNG_CONFIG_DIR", dir.path("config"))
        .output()
        .unwrap()
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

/// Ref and query over `rows` (chr2 gained in the query); infercnv's output
/// row names and per-row mean log-ratio.
fn profile(dir: &TmpDir, rows: &[Box<str>], extra: &[&str]) -> (Vec<Box<str>>, Vec<f32>) {
    let (r, q, out) = (
        dir.path("ref.zarr"),
        dir.path("query.zarr"),
        dir.path("out"),
    );
    write_backend(&r, rows, "r", 30, None);
    write_backend(&q, rows, "q", 30, Some("chr2"));
    let mut args = vec![
        "infercnv",
        "--ref",
        &r,
        "--window",
        "5",
        "--no-center",
        "--out",
        &out,
    ];
    args.extend_from_slice(extra);
    args.push(&q);
    let run = mung(dir, &args);
    assert!(run.status.success(), "{}", stderr(&run));
    let cnv = open_sparse_matrix_by_path(&format!("{out}.zarr.zip")).unwrap();
    let names = cnv.row_names().unwrap();
    let x = cnv
        .read_columns_dmatrix((0..cnv.num_columns().unwrap()).collect())
        .unwrap();
    let means = (0..names.len()).map(|i| x.row(i).mean()).collect();
    (names, means)
}

fn mean_on(names: &[Box<str>], means: &[f32], chr: &str) -> f32 {
    let idx: Vec<usize> = (0..names.len())
        .filter(|&i| names[i].starts_with(&format!("{chr}:")))
        .collect();
    assert!(!idx.is_empty(), "no {chr} rows in {names:?}");
    idx.iter().map(|&i| means[i]).sum::<f32>() / idx.len() as f32
}

#[test]
fn a_depth_matrix_needs_no_annotation_and_shows_the_gain() {
    let dir = TmpDir::new("depth");
    let rows = [bins("chr1", 20, 0), bins("chr2", 20, 0)].concat();
    let (names, means) = profile(&dir, &rows, &[]);
    assert!(
        names.iter().all(|n| n.starts_with("chr")),
        "bins keep their names"
    );
    let (m1, m2) = (
        mean_on(&names, &means, "chr1"),
        mean_on(&names, &means, "chr2"),
    );
    assert!(m2 - m1 > 0.3, "chr2 gain not seen: chr1 {m1}, chr2 {m2}");
}

#[test]
fn sex_chromosome_bins_follow_exclude_chr() {
    let dir = TmpDir::new("chrx");
    let rows = [
        bins("chr1", 20, 0),
        bins("chr2", 20, 0),
        bins("chrX", 20, 0),
    ]
    .concat();
    let (names, _) = profile(&dir, &rows, &[]);
    assert!(
        !names.iter().any(|n| n.starts_with("chrX:")),
        "chrX excluded by default"
    );
    let (names, _) = profile(&dir, &rows, &["--exclude-chr", "none"]);
    assert!(
        names.iter().any(|n| n.starts_with("chrX:")),
        "chrX kept on request"
    );
}

#[test]
fn bins_on_unplaced_contigs_need_no_annotation() {
    let dir = TmpDir::new("contigs");
    let mut rows = [bins("chr1", 20, 0), bins("chr2", 20, 0)].concat();
    rows.push("chrUn_KI270302v1:0-2274".into());
    rows.push("chr1_KI270706v1_random:0-175055".into());
    let (names, _) = profile(&dir, &rows, &[]);
    assert!(names.iter().any(|n| n.starts_with("chr1:")));
}

#[test]
fn files_binned_on_different_grids_are_refused() {
    let dir = TmpDir::new("grids");
    let (r, q) = (dir.path("ref.zarr"), dir.path("query.zarr"));
    write_backend(&r, &bins("chr1", 10, 0), "r", 5, None);
    write_backend(&q, &bins("chr1", 10, BIN / 2), "q", 5, None);
    let run = mung(
        &dir,
        &["infercnv", "--ref", &r, "--out", &dir.path("out"), &q],
    );
    assert!(!run.status.success());
    assert!(stderr(&run).contains("different grids"), "{}", stderr(&run));
}

#[test]
fn gene_rows_offline_and_uncached_say_what_to_do() {
    let dir = TmpDir::new("genes");
    let q = dir.path("query.zarr");
    // One interval among gene rows still makes a gene axis.
    let rows: Vec<Box<str>> = vec!["GENE1".into(), "GENE2".into(), "chr1:0-1000".into()];
    write_backend(&q, &rows, "q", 3, None);
    let run = mung(&dir, &["infercnv", "--out", &dir.path("out"), &q]);
    assert!(!run.status.success());
    let err = stderr(&run);
    assert!(
        err.contains("--gff") && err.contains("mung data fetch"),
        "{err}"
    );
}
