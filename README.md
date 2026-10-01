# mung-cnv

**M**alignancy **U**nmixing on **N**ormalized **G**enomes with **CNV** estimation

`mung` builds inferCNV-style copy-number profiles from single-cell expression
backends and calls donor-private clone strata for consumers such as
[`senna`](https://github.com/causalpathlab/legume-rs) / `pinto` (`--cnv-clones`).

| Layer | Name |
|-------|------|
| crates.io package | `mung-cnv` |
| CLI binary | `mung` |
| Rust library | `cnv` (`use cnv::…`) |

## Installation

### Prerequisites

- **Rust** (stable; MSRV 1.91+) — [rustup](https://rustup.rs)
- A C toolchain and OpenBLAS (pulled in via `legume-numeric` / `ndarray`):

  ```sh
  # Debian / Ubuntu
  sudo apt-get install build-essential pkg-config libopenblas-dev

  # macOS (Homebrew)
  brew install openblas
  ```

### Install from crates.io

```sh
cargo install mung-cnv
```

Optional features:

```sh
cargo install mung-cnv --features hdf5   # .h5 / .h5ad backends (needs libhdf5)
cargo install mung-cnv --features cuda   # NVIDIA CUDA via legume-numeric / data-beans
cargo install mung-cnv --features metal  # Apple Metal
```

### Install from GitHub

```sh
cargo install --git https://github.com/causalpathlab/mung-cnv.git
```

### Build from a local clone

```sh
git clone https://github.com/causalpathlab/mung-cnv.git
cd mung-cnv
cargo build --release
cargo install --path .
```

## Usage

```sh
mung infercnv --gff gencode.v46.gtf.gz \
  --ref Control1.zarr.zip Control2.zarr.zip \
  --out aml001.cnv AML001.zarr.zip

mung clones --from aml001.cnv.zarr.zip --out aml001
# writes aml001.clones.parquet for senna / pinto --cnv-clones
```

`mung --help` and `mung <subcommand> --help` document all flags.
`mung describe <subcommand>` prints the same flags as JSON for front ends
such as `senna run`, which start `mung` as a separate program.

## Library

Dependents (e.g. senna) take the published package and keep the `cnv` import path:

```toml
cnv = { version = "0.2", package = "mung-cnv" }
```

```rust
use cnv::clone_call::read_clone_table;
```

## Ecosystem

- [`data-beans`](https://crates.io/crates/data-beans) — sparse backends
- [`legume-numeric`](https://crates.io/crates/legume-numeric) — matrix / MCMC
- [`legume-genomic-types`](https://crates.io/crates/legume-genomic-types) — GFF loci
- [`legume-rs`](https://github.com/causalpathlab/legume-rs) — senna, pinto, and friends

## License

MIT
