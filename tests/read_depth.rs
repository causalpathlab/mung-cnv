//! `mung infercnv` on a read-depth matrix whose rows are genomic bins
//! (`chr:start-end`, as `faba read-depth` writes them): no annotation is
//! needed, and a gain shows as a positive log-ratio on its chromosome.

use data_beans::sparse_io::SparseIoBackend;
use data_beans::sparse_io::{create_sparse_from_triplets, open_sparse_matrix_by_path};
use std::path::{Path, PathBuf};
use std::process::Command;

const BIN: i64 = 1_000_000;
const N_BINS_PER_CHR: i64 = 20;

fn bins() -> Vec<Box<str>> {
    ["chr1", "chr2"]
        .iter()
        .flat_map(|c| {
            (0..N_BINS_PER_CHR).map(move |i| format!("{c}:{}-{}", i * BIN, (i + 1) * BIN))
        })
        .map(Into::into)
        .collect()
}

/// A backend of `n_cells` cells over [`bins`]; `gain` doubles chr2's depth.
fn write_backend(path: &Path, prefix: &str, n_cells: usize, gain: bool) {
    let rows = bins();
    let mut triplets = Vec::new();
    for j in 0..n_cells {
        for (i, name) in rows.iter().enumerate() {
            let base = 20.0 + ((i * 7 + j * 3) % 5) as f32;
            let x = if gain && name.starts_with("chr2:") {
                2.0 * base
            } else {
                base
            };
            triplets.push((i as u64, j as u64, x));
        }
    }
    let shape = (rows.len(), n_cells, triplets.len());
    let mut data = create_sparse_from_triplets(
        &triplets,
        shape,
        Some(path.to_str().unwrap()),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    data.register_row_names_vec(&rows);
    let cells: Vec<Box<str>> = (0..n_cells)
        .map(|j| format!("{prefix}{j}@D1").into())
        .collect();
    data.register_column_names_vec(&cells);
}

fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("mung-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn mung_offline(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_mung"))
        .args(args)
        .env("MUNG_OFFLINE", "1")
        .env("MUNG_CACHE_DIR", dir.join("cache"))
        .env("MUNG_CONFIG_DIR", dir.join("config"))
        .output()
        .unwrap()
}

#[test]
fn a_depth_matrix_needs_no_annotation_and_shows_the_gain() {
    let dir = tmp("depth");
    let (r, q, out) = (
        dir.join("ref.zarr"),
        dir.join("query.zarr"),
        dir.join("out"),
    );
    write_backend(&r, "r", 30, false);
    write_backend(&q, "q", 30, true);
    let s = |p: &Path| p.to_str().unwrap().to_string();
    let run = mung_offline(
        &dir,
        &[
            "infercnv",
            "--ref",
            &s(&r),
            "--window",
            "5",
            "--no-center",
            "--out",
            &s(&out),
            &s(&q),
        ],
    );
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );

    let cnv = open_sparse_matrix_by_path(&format!("{}.zarr.zip", s(&out))).unwrap();
    let rows = cnv.row_names().unwrap();
    let x = cnv
        .read_columns_dmatrix((0..cnv.num_columns().unwrap()).collect())
        .unwrap();
    let mean_on = |chr: &str| {
        let idx: Vec<usize> = (0..rows.len())
            .filter(|&i| rows[i].starts_with(chr))
            .collect();
        assert!(!idx.is_empty(), "no {chr} rows in {rows:?}");
        idx.iter().map(|&i| x.row(i).mean()).sum::<f32>() / idx.len() as f32
    };
    let (m1, m2) = (mean_on("1:"), mean_on("2:"));
    assert!(m2 - m1 > 0.3, "chr2 gain not seen: chr1 {m1}, chr2 {m2}");
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn gene_rows_offline_and_uncached_say_what_to_do() {
    let dir = tmp("genes");
    let q = dir.join("query.zarr");
    let rows: Vec<Box<str>> = vec!["GENE1".into(), "GENE2".into()];
    let mut data = create_sparse_from_triplets(
        &[(0, 0, 1.0), (1, 1, 2.0)],
        (2, 2, 2),
        Some(q.to_str().unwrap()),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    data.register_row_names_vec(&rows);
    data.register_column_names_vec(&["a".into(), "b".into()]);
    let out = dir.join("out");
    let run = mung_offline(
        &dir,
        &[
            "infercnv",
            "--out",
            out.to_str().unwrap(),
            q.to_str().unwrap(),
        ],
    );
    assert!(!run.status.success());
    let err = String::from_utf8_lossy(&run.stderr);
    assert!(
        err.contains("--gff") && err.contains("mung data fetch"),
        "{err}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
