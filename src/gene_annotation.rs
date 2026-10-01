//! The gene annotation (GTF) placing genes on the genome. A file named with
//! `--gff` wins; otherwise `--species` (or the config's `default`) picks an
//! entry of the annotation config, whose file is downloaded once into the
//! cache unless `MUNG_OFFLINE` is set. `mung data fetch` fills the cache
//! ahead of time for machines that will run offline.
//!
//! The config (`annotations.json`) names each species' source, release,
//! assembly and URL, so no address lives in the code. It is the user's copy
//! in `~/.config/mung/` (`MUNG_CONFIG_DIR` overrides) when there is one,
//! else the `data/annotations.json` this release was built with.

use anyhow::{Context, Result};
use log::info;
use serde_json::Value;
use std::fs;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::time::Duration;

const CONFIG: &str = "annotations.json";
/// The config this release ships, used when the user has none.
const SHIPPED: &str = include_str!("../data/annotations.json");
const OFFLINE_ENV: &str = "MUNG_OFFLINE";
const CACHE_DIR_ENV: &str = "MUNG_CACHE_DIR";
const CONFIG_DIR_ENV: &str = "MUNG_CONFIG_DIR";

/// One downloadable annotation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    pub species: Box<str>,
    pub url: Box<str>,
    /// Source, release and assembly, as the config gives them, for the log.
    pub label: Box<str>,
}

impl Annotation {
    /// Where the cache keeps it: under its species, by the URL's file name.
    #[must_use]
    pub fn cached(&self) -> Option<PathBuf> {
        let name = self.url.rsplit('/').next().unwrap_or(&self.url);
        cache_dir().map(|d| d.join(self.species.as_ref()).join(name))
    }

    /// The cached copy, downloaded first when missing (or always, with
    /// `force`).
    pub fn ensure_cached(&self, force: bool) -> Result<PathBuf> {
        let to = self.cached().context("no cache directory")?;
        if force || !to.is_file() {
            download_gzip(&self.url, &to)?;
        }
        Ok(to)
    }
}

/// The annotation config: every species it names, and the default one.
#[derive(Debug, Clone)]
pub struct Config {
    pub default: Option<Box<str>>,
    pub annotations: Vec<Annotation>,
    /// Where it was read from, or `None` for the shipped copy.
    pub path: Option<PathBuf>,
}

impl Config {
    /// The user's config, else the shipped one.
    pub fn load() -> Result<Self> {
        let user = std::env::var_os(CONFIG_DIR_ENV)
            .map(PathBuf::from)
            .or_else(|| dirs::home_dir().map(|h| h.join(".config").join("mung")))
            .map(|d| d.join(CONFIG))
            .filter(|p| p.is_file());
        match user {
            Some(p) => {
                let text =
                    fs::read_to_string(&p).with_context(|| format!("reading {}", p.display()))?;
                Self::parse(&text, Some(p))
            }
            None => Self::parse(SHIPPED, None),
        }
    }

    /// Where the config came from, for messages.
    #[must_use]
    pub fn source(&self) -> String {
        self.path
            .as_ref()
            .map_or_else(|| format!("built-in {CONFIG}"), |p| p.display().to_string())
    }

    fn parse(text: &str, path: Option<PathBuf>) -> Result<Self> {
        let mut config = Self {
            default: None,
            annotations: Vec::new(),
            path,
        };
        let v: Value =
            serde_json::from_str(text).with_context(|| format!("parsing {}", config.source()))?;
        let entries = v["annotations"]
            .as_object()
            .with_context(|| format!("{}: no `annotations` object", config.source()))?;
        for (species, e) in entries {
            let url = e["url"]
                .as_str()
                .with_context(|| format!("{}: `{species}` has no `url`", config.source()))?;
            let label = ["source", "release", "assembly"]
                .iter()
                .filter_map(|k| e[*k].as_str())
                .collect::<Vec<_>>()
                .join(" ");
            config.annotations.push(Annotation {
                species: species.as_str().into(),
                url: url.into(),
                label: label.into(),
            });
        }
        config.default = v["default"].as_str().map(Into::into);
        Ok(config)
    }

    /// The entry for `species`, else the config's default.
    pub fn pick(&self, species: Option<&str>) -> Result<&Annotation> {
        let name = species
            .or(self.default.as_deref())
            .with_context(|| format!("{} has no `default`; pass --species", self.source()))?;
        self.annotations
            .iter()
            .find(|a| a.species.as_ref() == name)
            .with_context(|| {
                let known: Vec<&str> = self
                    .annotations
                    .iter()
                    .map(|a| a.species.as_ref())
                    .collect();
                format!(
                    "{} has no `{name}` annotation (it has: {}); pass --gff",
                    self.source(),
                    known.join(", ")
                )
            })
    }
}

fn cache_dir() -> Option<PathBuf> {
    std::env::var_os(CACHE_DIR_ENV)
        .map(PathBuf::from)
        .or_else(|| dirs::cache_dir().map(|d| d.join("mung")))
}

/// The annotation to read: `gff` when given, else `species`' (or the
/// default's) cached copy, downloaded first if missing and `MUNG_OFFLINE`
/// is not set.
pub fn resolve(gff: Option<&str>, species: Option<&str>) -> Result<Box<str>> {
    if let Some(gff) = gff {
        return Ok(gff.into());
    }
    let config = Config::load()?;
    let a = config.pick(species)?;
    let cached = a.cached().filter(|p| p.is_file());
    let path = match cached {
        Some(p) => p,
        None => {
            anyhow::ensure!(
                std::env::var_os(OFFLINE_ENV).is_none(),
                "{OFFLINE_ENV} is set and the `{}` annotation is not cached; pass --gff or run `mung data fetch`",
                a.species
            );
            a.ensure_cached(false).with_context(|| {
                format!(
                    "downloading the `{}` annotation; pass --gff instead",
                    a.species
                )
            })?
        }
    };
    info!(
        "gene annotation: {} ({}) at {}",
        a.species,
        a.label,
        path.display()
    );
    Ok(path.to_string_lossy().into())
}

/// Download gzip file `url` to `to`, streamed to a side file and renamed
/// once whole, so a half-finished download is never read. A response that
/// does not start as gzip is refused before anything is written.
fn download_gzip(url: &str, to: &Path) -> Result<()> {
    info!("downloading {url}");
    let response = ureq::AgentBuilder::new()
        .timeout_connect(Duration::from_secs(30))
        .timeout_read(Duration::from_secs(60))
        .build()
        .get(url)
        .call()
        .with_context(|| format!("requesting {url}"))?;
    let mut body = std::io::BufReader::new(response.into_reader());
    anyhow::ensure!(
        body.fill_buf()?.starts_with(&[0x1f, 0x8b]),
        "{url} is not gzip"
    );
    fs::create_dir_all(to.parent().context("no cache directory")?)?;
    let tmp = to.with_extension("part");
    let mut file = fs::File::create(&tmp)?;
    std::io::copy(&mut body, &mut file).context("reading the response")?;
    drop(file);
    fs::rename(&tmp, to)?;
    info!("cached {}", to.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEXT: &str = r#"{
        "default": "sp1",
        "annotations": {
            "sp1": {"source": "SRC", "release": "1", "url": "https://example.org/a/sp1.gtf.gz"},
            "sp2": {"url": "https://example.org/b/sp2.gtf.gz"}
        }
    }"#;

    fn config() -> Config {
        Config::parse(TEXT, None).unwrap()
    }

    #[test]
    fn picks_the_default_or_the_named_species() {
        let c = config();
        assert_eq!(c.pick(None).unwrap().species.as_ref(), "sp1");
        assert_eq!(c.pick(Some("sp2")).unwrap().species.as_ref(), "sp2");
        assert_eq!(c.pick(None).unwrap().label.as_ref(), "SRC 1");
    }

    #[test]
    fn an_unknown_species_names_the_known_ones() {
        let err = config().pick(Some("sp3")).unwrap_err().to_string();
        assert!(err.contains("sp1") && err.contains("sp2"), "{err}");
    }

    #[test]
    fn an_entry_without_a_url_is_an_error() {
        let bad = r#"{"annotations": {"sp1": {"release": "1"}}}"#;
        assert!(Config::parse(bad, None).is_err());
    }

    #[test]
    fn the_cache_keys_on_species_and_file_name() {
        let p = config().pick(None).unwrap().cached().unwrap();
        assert!(
            p.ends_with(Path::new("sp1").join("sp1.gtf.gz")),
            "{}",
            p.display()
        );
    }

    #[test]
    fn a_named_file_wins() {
        let got = resolve(Some("genes.gtf"), None).unwrap();
        assert_eq!(got.as_ref(), "genes.gtf");
    }

    #[test]
    fn the_shipped_config_has_a_default() {
        let c = Config::parse(SHIPPED, None).unwrap();
        assert!(c.pick(None).is_ok());
    }
}
