//! `mung` (Malignancy Unmixing on Normalized Genomes with CNV estimation) —
//! copy-number profiles from single-cell expression backends
//! (package: `mung-cnv`; lib: `cnv`).

use anyhow::Context;
use clap::{Args, Parser, Subcommand};
use cnv::cell_profile::{run_cell_profiles, CellProfileConfig};
use cnv::clone_bayes::{cells_table_beside, DEFAULT_MIN_PURITY};
use cnv::clone_call::{call_clones_with_burden, write_clone_table, CloneCallConfig, CloneEngine};
use cnv::gene_annotation::{self, Config};
use cnv::gene_loci::{check_one_grid, interval_loci, GeneLocusIndex};
use data_beans::aux::data_loading::{read_data_on_shared_rows, ReadSharedRowsArgs};
use data_beans::aux::feature_names::FeatureNameKind;
use data_beans::convert::try_open_or_convert;
use data_beans::sparse_io_vector::SparseIoVec;
use log::{info, warn};
use std::sync::Arc;

#[derive(Parser, Debug)]
#[command(
    name = "mung",
    version,
    about = "mung — Malignancy Unmixing on Normalized Genomes with CNV estimation",
    long_about = "Reads `data-beans` backends (.zarr.zip / .zarr / .h5) and writes\n\
                  copy-number profiles and clone strata for `--cnv-clones` consumers.\n\
                  \n\
                  Subcommands:\n  \
                  \x20 infercnv — per-cell inferCNV log-ratio profiles on a genomic-interval axis\n  \
                  \x20 clones   — donor-private CNV strata (`{out}.clones.parquet`) for\n               \
                  senna / pinto `--cnv-clones`"
)]
struct Cli {
    #[arg(short = 'v', long, global = true, help = "Increase output verbosity")]
    verbose: bool,

    #[arg(
        long = "n-threads",
        visible_aliases = ["threads", "num-threads"],
        global = true,
        value_name = "N",
        help = "Limit the number of CPU threads (rayon global pool)"
    )]
    n_threads: Option<usize>,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    #[command(
        about = "Per-cell inferCNV profiles against a normal reference",
        long_about = "For every query cell: depth-normalise, take the log-ratio against the\n\
                      per-gene mean of the reference (normal) cells, clip, average over a\n\
                      sliding window of genes within each chromosome, and median-centre.\n\
                      \n\
                      Output is a new data-beans backend `{out}.zarr.zip` whose rows are\n\
                      genomic intervals (`chr:start-end`; one gene each, or `--bin-size`\n\
                      tiles — prefer `1000000` / 1 Mb on large cohorts) and whose columns\n\
                      are the query cells. `{out}.features.parquet`\n\
                      maps each row back to its genes; `{out}.cells.parquet` has per-cell\n\
                      depth and mean |log-ratio| (CNV burden).\n\
                      \n\
                      With no `--ref`, the query cohort mean is the baseline: any CNV shared\n\
                      by every cell becomes invisible.\n\
                      \n\
                      Example:\n  \
                      mung infercnv \\\n    \
                      --ref Control1.zarr.zip Control2.zarr.zip \\\n    \
                      --out sample1.cnv sample1.zarr.zip"
    )]
    Infercnv(InferCnvArgs),
    #[command(
        about = "Donor-private CNV clone strata from inferCNV profiles",
        long_about = "Gate cells on a genomic sketch + CNV burden into stratum 0\n\
                      (mixable) vs donor-private clones for `--cnv-clones`.\n\
                      \n\
                      Default `--engine bayes`: a cell is a clone member only if\n\
                      posterior malignancy is high and it sits in a peaky donor-mix\n\
                      component. False clones hard-partition collapse and under-\n\
                      integrate batch δ — costly; the Bayes gate is the primary\n\
                      bar and works without `--ref`. Legacy `--engine mixture`\n\
                      keeps the cluster-score BIC path.\n\
                      \n\
                      Generative model (bayes): each cell is mixable (m=0) or\n\
                      malignant (m=1); burden is low/high LogNormal given m;\n\
                      sketch ~ N(μ0,σ²) if mixable, else N(μ_z,σ²) for clone z\n\
                      with peaky donor mix φ_z. Report p_malig = E[m], then MAP z.\n\
                      \n\
                      Sketch: `--bin-size` defaults to `0` (chromosome means — fast\n\
                      for Bayes). Override with Mb tiles (e.g. `1000000`) if you want;\n\
                      that still applies to inferCNV when run first.\n\
                      \n\
                      Reads an existing CNV backend (`--from`), or runs `infercnv`\n\
                      first on `--ref` / QUERY and then calls clones. Burden comes\n\
                      from `{prefix}.cells.parquet` when present, else mean |CNV|.\n\
                      \n\
                      Writes `{out}.clones.parquet` (cell, donor, cluster, stratum,\n\
                      …, burden, p_malig). Pass that file as `--cnv-clones` to\n\
                      senna / pinto so collapse cannot mix across clone boundaries."
    )]
    Clones(CloneArgs),
    /// A subcommand's flags as JSON, for front ends (`senna run`) that
    /// build a form from them and start `mung` as a separate program.
    /// Reference data mung downloads: where it is cached, and fetching it
    /// ahead of time for machines that will run offline.
    #[command(subcommand)]
    Data(DataCmd),
    #[command(hide = true)]
    Describe {
        /// The subcommand to describe, e.g. `clones`.
        command: String,
    },
}

#[derive(Subcommand, Debug)]
enum DataCmd {
    /// Print the annotation config in use, its species, and which are cached
    Where,
    /// Download gene annotations named in the config into the cache
    Fetch {
        #[arg(long, num_args = 1.., help = "Species to fetch; default: all in the config")]
        species: Vec<Box<str>>,
        #[arg(long, help = "Replace what the cache already holds")]
        force: bool,
    },
}

#[derive(Args, Debug, Clone)]
struct InferCnvArgs {
    #[arg(
        required = true,
        value_name = "QUERY",
        help = "Query backend files (.zarr.zip / .zarr / .h5); cells to profile"
    )]
    query: Vec<Box<str>>,

    #[arg(
        short,
        long,
        help = "Output prefix; writes {out}.zarr.zip, {out}.features.parquet, {out}.cells.parquet"
    )]
    out: Box<str>,

    #[arg(
        long,
        default_value_t = 0,
        help = "Average smoothed genes into fixed genomic tiles of this many bp (0 = one row per gene)",
        long_help = "Average smoothed genes into fixed genomic tiles of this many bp.\n\
                     Default 0 keeps classic inferCNV gene-level rows.\n\
                     For large cohorts, prefer `--bin-size 1000000` (1 Mb):\n\
                     after a 101-gene window the signal is already ~Mb-scale,\n\
                     and autosomes compress to ~3k rows instead of ~15–20k genes."
    )]
    bin_size: i64,

    #[command(flatten)]
    opts: InferOpts,
}

/// The inferCNV flags `infercnv` and `clones` share.
#[derive(Args, Debug, Clone)]
struct InferOpts {
    #[arg(
        long,
        num_args = 1..,
        value_name = "REF",
        help = "Reference (normal) backend files; their per-gene mean is the diploid baseline"
    )]
    r#ref: Vec<Box<str>>,

    #[arg(
        short = 'g',
        long = "gff",
        help = "Gene annotation file (GFF/GTF)",
        long_help = "Path to gene annotation file in GFF/GTF format, with `gene`\n\
                     features carrying gene_id and gene_name. Places each gene on\n\
                     the genome. Not needed for rows named as genomic intervals\n\
                     (`chr:start-end`, e.g. `faba depth` bins). Optional: without\n\
                     it, the `--species` entry of the annotation config is used,\n\
                     downloaded once into the cache (`mung data where`). Ignored\n\
                     by `clones --from`."
    )]
    gff: Option<Box<str>>,

    #[arg(
        long,
        help = "Species whose annotation to use without --gff",
        long_help = "Species in the annotation config (`mung data where`) whose\n\
                     gene annotation to use when --gff is not given. Default: the\n\
                     config's `default`."
    )]
    species: Option<Box<str>>,

    #[arg(
        long,
        default_value_t = 101,
        help = "Smoothing window in genes (odd; ≤1 disables)"
    )]
    window: usize,

    #[arg(
        long,
        default_value_t = 3.0,
        help = "Clip per-gene log-ratio at ±clip before smoothing (≤0 disables)"
    )]
    clip: f32,

    #[arg(
        long,
        default_value_t = 1e4,
        help = "Depth-normalisation target: ln(1 + scale·x/depth)"
    )]
    scale: f32,

    #[arg(
        long,
        default_value_t = 0.1,
        help = "Keep genes whose mean raw count over reference cells is ≥ this"
    )]
    min_mean_expr: f32,

    #[arg(long, default_value_t = 1000, help = "Cells per streamed block")]
    block_size: usize,

    #[arg(long, help = "Do not median-centre each cell after smoothing")]
    no_center: bool,

    #[arg(
        long,
        value_delimiter = ',',
        default_value = "chrX,chrY,chrM",
        value_name = "CHR[,CHR..]",
        help = "Chromosomes to drop from the gene axis (inferCNV default X,Y,M; `none` keeps all)"
    )]
    exclude_chr: Vec<Box<str>>,

    #[arg(long, help = "Preload backend columns into memory")]
    preload: bool,
}

#[derive(Args, Debug, Clone)]
struct CloneArgs {
    #[arg(
        value_name = "QUERY",
        required_unless_present = "from",
        help = "Query backends (expression); ignored when `--from` is set"
    )]
    query: Vec<Box<str>>,

    #[arg(
        long,
        help = "Existing inferCNV backend (`.zarr.zip`); skip the expression pass"
    )]
    from: Option<Box<str>>,

    #[arg(
        short,
        long,
        help = "Output prefix; writes {out}.clones.parquet (and inferCNV artifacts when not `--from`)"
    )]
    out: Box<str>,

    #[arg(
        long,
        default_value_t = 0,
        help = "Genomic tile size in bp for inferCNV rows and the clone sketch (0 = inferCNV default)",
        long_help = "Same `--bin-size` as `mung infercnv` (default 0 = classic inferCNV).\n\
                     When running inferCNV first: average smoothed genes into fixed genomic\n\
                     tiles of this many bp (0 = one row per gene). Prefer `1000000` (1 Mb)\n\
                     on large cohorts for the CNV backend itself.\n\
                     For the clone sketch (also with `--from`): default `0` = one dim per\n\
                     chromosome (recommended for Bayes — Gibbs scales with C).\n\
                     `>0` = one dim per tile of the interval midpoint."
    )]
    bin_size: i64,

    #[command(flatten)]
    infer: InferOpts,

    #[arg(
        long,
        value_enum,
        default_value_t = CloneEngine::Bayes,
        help = "Clone gate: bayes (default) or legacy mixture",
        long_help = "bayes (default): m ~ Bern; b|m LogNormal; sketch|m,z Gaussian;\n\
                     clone z has peaky donor mix φ. mixture: legacy cluster-score\n\
                     BIC gate. False clones under-integrate batch δ — prefer bayes."
    )]
    engine: CloneEngine,

    #[arg(
        long,
        default_value_t = 8,
        help = "Finite K_max (Bayes components / mixture k-means K)",
        long_help = "Finite K_max on the genomic sketch (`--bin-size` dims).\n\
                     Bayes: empty components OK. Mixture: overclustering dumps\n\
                     weak clusters to stratum 0. False clones under-integrate δ."
    )]
    k_max: usize,
    #[arg(
        long,
        help = "Donor-mix purity floor (Bayes default 0.8; mixture omit = mixture only)"
    )]
    min_purity: Option<f32>,
    #[arg(
        long,
        default_value_t = 50,
        help = "Minimum cells to keep a clone component / cluster"
    )]
    min_cells: usize,
    #[arg(
        long = "segmental-z",
        visible_alias = "spatial-z",
        help = "Optional segmental-CN z floor (mixture engine only)",
        long_help = "Optional floor on the segmental (genomic) z-score of a cluster's\n\
                     chromosome-sketch L1 vs a size-matched null. Mixture only.\n\
                     `--spatial-z` is a deprecated alias."
    )]
    spatial_z: Option<f32>,
    #[arg(
        long,
        default_value_t = 32,
        help = "Null permutations for segmental-CN z (mixture only)"
    )]
    n_perm: usize,
    #[arg(long, default_value_t = 100, help = "Bayes Gibbs sweeps after burn-in")]
    n_sweeps: usize,
    #[arg(long, default_value_t = 40, help = "Bayes Gibbs burn-in sweeps")]
    n_burnin: usize,
    #[arg(
        long,
        default_value_t = 0.5,
        help = "Stratum 0 if posterior malignancy is below this (Bayes)"
    )]
    p_malig_threshold: f32,
    #[arg(long, default_value_t = 1, help = "Random seed")]
    seed: u64,
}

/// How rows are placed on the genome.
enum RowAxis {
    /// Gene rows, looked up in the annotation.
    Genes(GeneLocusIndex),
    /// Genomic-interval rows, placed by their own names.
    Intervals,
}

fn run_infercnv(
    query: &[Box<str>],
    out: &str,
    bin_size: i64,
    args: &InferOpts,
) -> anyhow::Result<()> {
    let has_ref = !args.r#ref.is_empty();
    if !has_ref {
        warn!("no --ref given: using the query cohort mean as baseline; shared CNVs will be invisible");
    }

    // Reference backends first, then query, in one shared-row load.
    let mut files: Vec<Box<str>> = args.r#ref.clone();
    files.extend(query.iter().cloned());

    // The row axis decides how rows are placed, so it is read first: genes
    // need the annotation (found now, before any data is loaded), while a
    // matrix whose every row is a genomic interval (faba depth) is placed by
    // its own names and aligned across files by exact name.
    let mut raw_rows: Vec<Box<str>> = Vec::new();
    for f in &files {
        raw_rows.extend(try_open_or_convert(f)?.row_names()?);
    }
    let axis = match interval_loci(&raw_rows) {
        Some(bins) => {
            check_one_grid(&bins)?;
            if args.gff.is_some() || args.species.is_some() {
                warn!("rows are genomic intervals: --gff / --species are not used");
            }
            let mut widths: Vec<i64> = bins.iter().map(|b| b.stop - b.start + 1).collect();
            widths.sort_unstable();
            let width = widths.get(widths.len() / 2).copied().unwrap_or(0);
            info!(
                "rows are genomic intervals of about {width} bp: --window {} smooths over about {} bp",
                args.window,
                width.saturating_mul(args.window as i64)
            );
            RowAxis::Intervals
        }
        None => {
            let gff = gene_annotation::resolve(args.gff.as_deref(), args.species.as_deref())?;
            RowAxis::Genes(
                GeneLocusIndex::from_gff(&gff).with_context(|| format!("reading {gff}"))?,
            )
        }
    };

    let loaded = read_data_on_shared_rows(ReadSharedRowsArgs {
        data_files: files,
        preload: args.preload,
        // The reference/query split below is positional per file, and every
        // query cell gets a profile.
        keep_empty_barcodes: true,
        // Bins keep their names as written (data-beans' locus canonical form
        // would rename them and merge overlapping ones).
        feature_kind: matches!(axis, RowAxis::Intervals).then_some(FeatureNameKind::Exact),
        ..Default::default()
    })?;
    let data = loaded.data;

    let per_backend = data.num_columns_by_data()?;
    let n_ref_cols: usize = per_backend.iter().take(args.r#ref.len()).sum();
    let n_total = data.num_columns();
    let (ref_cols, query_cols): (Vec<usize>, Vec<usize>) = if has_ref {
        ((0..n_ref_cols).collect(), (n_ref_cols..n_total).collect())
    } else {
        ((0..n_total).collect(), (0..n_total).collect())
    };
    info!(
        "loaded {} genes; {} reference cells, {} query cells",
        data.num_rows(),
        ref_cols.len(),
        query_cols.len()
    );

    let row_names = data.row_names()?;
    let loci: Vec<Option<_>> = match &axis {
        RowAxis::Genes(index) => index.resolve_all(&row_names),
        RowAxis::Intervals => interval_loci(&row_names)
            .context("aligned rows are no longer genomic intervals")?
            .into_iter()
            .map(Some)
            .collect(),
    };

    let cfg = CellProfileConfig {
        window: args.window,
        clip: args.clip,
        scale: args.scale,
        min_mean_expr: args.min_mean_expr,
        bin_size,
        block_size: args.block_size,
        center: !args.no_center,
        exclude_chr: args
            .exclude_chr
            .iter()
            .filter(|c| !c.is_empty() && !c.eq_ignore_ascii_case("none"))
            .cloned()
            .collect(),
    };
    let outs = run_cell_profiles(&data, &ref_cols, &query_cols, &loci, &cfg, out)?;
    info!(
        "done: {} ({} intervals × {} cells)",
        outs.backend, outs.n_rows, outs.n_cells
    );
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    data_beans::aux::logging::init_logger(cli.verbose);

    if let Some(n) = cli.n_threads {
        anyhow::ensure!(n >= 1, "--n-threads must be >= 1");
        if let Err(e) = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
        {
            warn!("--n-threads {n} ignored: the rayon global pool already exists ({e})");
        }
    }

    match &cli.command {
        Commands::Infercnv(a) => run_infercnv(&a.query, &a.out, a.bin_size, &a.opts),
        Commands::Clones(args) => run_clones(args),
        Commands::Data(cmd) => run_data(cmd),
        Commands::Describe { command } => {
            println!("{}", describe(command)?);
            Ok(())
        }
    }
}

/// Subcommand `name`'s arguments as JSON: what a front end needs to show
/// each flag, fill it in, and lay out the command line.
fn describe(name: &str) -> anyhow::Result<String> {
    use clap::{ArgAction, CommandFactory};
    let mut cli = Cli::command();
    cli.build();
    let cmd = cli
        .find_subcommand(name)
        .filter(|c| !c.is_hide_set())
        .ok_or_else(|| anyhow::anyhow!("mung has no `{name}` command"))?;
    let text = |s: Option<&clap::builder::StyledStr>| s.map(ToString::to_string);
    let args: Vec<serde_json::Value> = cmd
        .get_arguments()
        .map(|a| {
            let action = match a.get_action() {
                ArgAction::SetTrue => "set_true",
                ArgAction::SetFalse => "set_false",
                ArgAction::Set => "set",
                ArgAction::Append => "append",
                ArgAction::Count => "count",
                ArgAction::Help | ArgAction::HelpShort | ArgAction::HelpLong => "help",
                ArgAction::Version => "version",
                _ => "other",
            };
            let num_args = a.get_num_args();
            // What a value must parse as, so a front end can check it.
            let t = a.get_value_parser().type_id();
            let is = |ids: &[std::any::TypeId]| ids.iter().any(|x| t == *x);
            use std::any::TypeId as T;
            let value_type = if is(&[T::of::<usize>(), T::of::<u64>(), T::of::<u32>()]) {
                "unsigned"
            } else if is(&[T::of::<i64>(), T::of::<i32>()]) {
                "integer"
            } else if is(&[T::of::<f64>(), T::of::<f32>()]) {
                "number"
            } else {
                "text"
            };
            serde_json::json!({
                "id": a.get_id().as_str(),
                "long": a.get_long(),
                "short": a.get_short().map(String::from),
                "positional": a.is_positional(),
                "help": text(a.get_help()),
                "long_help": text(a.get_long_help()),
                "action": action,
                "value_type": value_type,
                "values": a
                    .get_possible_values()
                    .iter()
                    .filter(|v| !v.is_hide_set())
                    .map(|v| v.get_name().to_string())
                    .collect::<Vec<_>>(),
                "default": a
                    .get_default_values()
                    .iter()
                    .map(|v| v.to_string_lossy().into_owned())
                    .collect::<Vec<_>>(),
                "delimiter": a.get_value_delimiter().map(String::from),
                "num_args_min": num_args.map(|r| r.min_values()),
                "num_args_max": num_args.and_then(|r| {
                    (r.max_values() != usize::MAX).then_some(r.max_values())
                }),
                "hidden": a.is_hide_set(),
                "required": a.is_required_set(),
                "global": a.is_global_set(),
            })
        })
        .collect();
    Ok(serde_json::json!({
        "describe": 1,
        "program": "mung",
        "version": env!("CARGO_PKG_VERSION"),
        "command": name,
        "about": text(cmd.get_about()),
        "args": args,
    })
    .to_string())
}

fn run_data(cmd: &DataCmd) -> anyhow::Result<()> {
    let config = Config::load()?;
    match cmd {
        DataCmd::Where => {
            println!("config: {}", config.source());
            for a in &config.annotations {
                let state = match a.cached() {
                    Some(p) if p.is_file() => p.display().to_string(),
                    _ => "not cached".into(),
                };
                let mark = if config.default.as_deref() == Some(&a.species) {
                    " (default)"
                } else {
                    ""
                };
                println!("  {}{mark}: {}; {state}", a.species, a.label);
            }
        }
        DataCmd::Fetch { species, force } => {
            let picked: Vec<_> = if species.is_empty() {
                config.annotations.iter().collect()
            } else {
                species
                    .iter()
                    .map(|s| config.pick(Some(s)))
                    .collect::<anyhow::Result<_>>()?
            };
            for a in picked {
                info!("{}: {}", a.species, a.ensure_cached(*force)?.display());
            }
        }
    }
    Ok(())
}

fn run_clones(args: &CloneArgs) -> anyhow::Result<()> {
    let backend_path = if let Some(from) = args.from.as_ref() {
        from.to_string()
    } else {
        run_infercnv(&args.query, &args.out, args.bin_size, &args.infer)?;
        format!("{}.zarr.zip", args.out)
    };

    let opened = try_open_or_convert(&backend_path)?;
    let mut data = SparseIoVec::new();
    data.push(Arc::from(opened), None)?;

    let min_purity = args.min_purity.or(match args.engine {
        CloneEngine::Bayes => Some(DEFAULT_MIN_PURITY),
        CloneEngine::Mixture => None,
    });
    let cfg = CloneCallConfig {
        bin_size: args.bin_size,
        k_max: args.k_max,
        engine: args.engine,
        min_purity,
        min_cells: args.min_cells,
        spatial_z: args.spatial_z,
        n_perm: args.n_perm,
        seed: args.seed,
        n_sweeps: args.n_sweeps,
        n_burnin: args.n_burnin,
        p_malig_threshold: args.p_malig_threshold,
        ..Default::default()
    };
    // Only the backend's own sibling table: a `{out}.cells.parquet` left by
    // an earlier run on a different dataset must not feed `--from`.
    let beside = cells_table_beside(&backend_path);
    let rows = call_clones_with_burden(&data, &cfg, None, &[beside.as_str()])?;
    let n_clone = rows.iter().filter(|r| r.stratum > 0).count();
    let n_strata = rows.iter().map(|r| r.stratum).max().unwrap_or(0);
    let path = format!("{}.clones.parquet", args.out);
    write_clone_table(&rows, &path)?;
    info!(
        "wrote {path} ({} cells, {n_clone} in {} donor-private clone(s); engine={:?})",
        rows.len(),
        n_strata,
        args.engine
    );
    Ok(())
}
