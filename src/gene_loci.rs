//! Resolve data-beans row names to genomic loci.
//!
//! A row axis made only of genomic intervals (`chr:start-end`, 0-based
//! half-open, as `faba depth` names its bins) places each row by its own
//! coordinates ([`interval_loci`]). Any other row is a
//! gene, looked up in a GFF/GTF. Gene row names come in a few shapes:
//! `{ensg}_{symbol}`, a bare symbol, a bare Ensembl id, or any of these with
//! a `/modality/...` suffix (faba).
//!
//! Matching goes through the workspace's one canonical gene matcher,
//! [`data_beans::utilities::name_matching::GeneIndex`], built over the GFF
//! vocabulary as `{ensg}_{symbol}`: case-insensitive; exact → symbol →
//! Ensembl id → decomposed `ENSG_SYMBOL[/aux]` → HGNC alias table (an
//! old symbol finds its renamed gene) → flexible fallback. The
//! Ensembl version suffix is dropped on both sides first.
//!
//! Same source of truth as `faba` (`genomic_data::gff::GffRecordMap`,
//! `gene` features stretched over all their records).

use data_beans::utilities::name_matching::GeneIndex;
use genomic_data::coordinates::{parse_interval, PeakCoord};
use genomic_data::gff::{GeneId, GeneSymbol, GffRecordMap};
use rayon::prelude::*;

/// One gene's locus, from the GFF `gene` feature. Coordinates are GFF
/// 1-based inclusive `[start, stop]`.
#[derive(Debug, Clone)]
pub struct GeneLocus {
    pub chromosome: Box<str>,
    pub start: i64,
    pub stop: i64,
    /// Transcription start: `start` on `+`, `stop` on `-`.
    pub tss: i64,
    pub gene_id: Box<str>,
    pub symbol: Box<str>,
}

impl GeneLocus {
    /// BED-style interval name `chr:start0-end` (0-based half-open) — the
    /// `chr:start-end` grammar `genomic_data::coordinates::parse_region`
    /// reads back.
    pub fn interval_name(&self) -> Box<str> {
        format!("{}:{}-{}", self.chromosome, self.start - 1, self.stop).into()
    }

    /// `{gene_id}_{symbol}`, faba's gene key; an interval row's own name.
    pub fn gene_key(&self) -> Box<str> {
        if self.symbol.is_empty() {
            return self.gene_id.clone();
        }
        format!("{}_{}", self.gene_id, self.symbol).into()
    }

    /// The locus of interval `r` (0-based half-open), named `name` and
    /// placed at its midpoint.
    fn from_coord(name: &str, r: PeakCoord) -> Self {
        Self {
            tss: r.start + 1 + (r.end - r.start - 1) / 2,
            start: r.start + 1,
            stop: r.end,
            gene_id: name.into(),
            symbol: "".into(),
            chromosome: r.chr,
        }
    }
}

/// Drop a faba `/modality/...` suffix and surrounding blanks.
fn strip_modality(row_name: &str) -> &str {
    row_name.split('/').next().unwrap_or(row_name).trim()
}

/// Each row's locus when every row is named as a genomic interval
/// (`chr:start-end` or `chr_start_end`, 0-based half-open, with an optional
/// faba `/modality/...` suffix); `None` when any row is not, and the axis is
/// genes. A row keeps its own name and its chromosome as written.
pub fn interval_loci(row_names: &[Box<str>]) -> Option<Vec<GeneLocus>> {
    row_names
        .iter()
        .map(|n| {
            let name = strip_modality(n);
            Some(GeneLocus::from_coord(name, parse_interval(name)?))
        })
        .collect()
}

/// An error naming two distinct bins that overlap: the files were binned on
/// different grids, and their rows cannot be aligned one to one.
pub fn check_one_grid(loci: &[GeneLocus]) -> anyhow::Result<()> {
    let mut by_pos: Vec<&GeneLocus> = loci.iter().collect();
    by_pos.sort_by(|a, b| (&a.chromosome, a.start, a.stop).cmp(&(&b.chromosome, b.start, b.stop)));
    by_pos.dedup_by(|a, b| a.chromosome == b.chromosome && a.start == b.start && a.stop == b.stop);
    for w in by_pos.windows(2) {
        if w[0].chromosome == w[1].chromosome && w[1].start <= w[0].stop {
            anyhow::bail!(
                "bins {} and {} overlap: the inputs were binned on different grids; \
                 bin every file with the same resolution",
                w[0].gene_id,
                w[1].gene_id
            );
        }
    }
    Ok(())
}

/// Lookup from GFF through the canonical [`GeneIndex`] matcher: a row name
/// in any of the supported shapes resolves to one [`GeneLocus`].
pub struct GeneLocusIndex {
    /// One per GFF gene, in coordinate order (first wins on duplicate keys).
    loci: Vec<GeneLocus>,
    /// Built over `{ensg}_{symbol}` keys, parallel to `loci`.
    index: GeneIndex,
}

/// Drop a faba `/modality/...` suffix and the Ensembl version on the
/// leading segment (`ENSG00000000003.15_TSPAN6` → `ENSG00000000003_TSPAN6`),
/// so the version never blocks a match.
fn normalize_query(row_name: &str) -> String {
    let core = strip_modality(row_name);
    let (head, rest) = match core.split_once('_') {
        Some((h, r)) => (h, Some(r)),
        None => (core, None),
    };
    let head = if head.len() >= 4 && head[..4].eq_ignore_ascii_case("ensg") {
        head.split('.').next().unwrap_or(head)
    } else {
        head
    };
    match rest {
        Some(r) => format!("{head}_{r}"),
        None => head.to_string(),
    }
}

impl GeneLocusIndex {
    /// Parse the GFF/GTF (`gene` features only) and build the index.
    pub fn from_gff(gff_file: &str) -> anyhow::Result<Self> {
        let map = GffRecordMap::from(gff_file)?;
        Ok(Self::from_record_map(&map))
    }

    pub fn from_record_map(map: &GffRecordMap) -> Self {
        let mut loci = Vec::new();
        let mut keys: Vec<Box<str>> = Vec::new();
        for rec in map.records() {
            let gene_id: Box<str> = match &rec.gene_id {
                GeneId::Ensembl(id) => id.clone(),
                GeneId::Missing => continue,
            };
            let symbol: Box<str> = match &rec.gene_name {
                GeneSymbol::Symbol(s) => s.clone(),
                GeneSymbol::Missing => gene_id.clone(),
            };
            let tss = match rec.strand {
                genomic_data::sam::Strand::Forward => rec.start,
                genomic_data::sam::Strand::Backward => rec.stop,
            };
            keys.push(format!("{gene_id}_{symbol}").into());
            loci.push(GeneLocus {
                chromosome: rec.seqname.clone(),
                start: rec.start,
                stop: rec.stop,
                tss,
                gene_id,
                symbol,
            });
        }
        let index = GeneIndex::build(&keys);
        log::info!("GFF: {} genes indexed as {{ensg}}_{{symbol}}", loci.len());
        Self { loci, index }
    }

    /// Resolve one row name through the canonical matcher.
    pub fn resolve(&self, row_name: &str) -> Option<&GeneLocus> {
        let q = normalize_query(row_name);
        if q.is_empty() {
            return None;
        }
        self.index.match_gene(&q).map(|i| &self.loci[i])
    }

    /// Resolve every row name (in parallel); `None` where a row has no GFF gene.
    pub fn resolve_all(&self, row_names: &[Box<str>]) -> Vec<Option<GeneLocus>> {
        let out: Vec<Option<GeneLocus>> = row_names
            .par_iter()
            .map(|n| self.resolve(n).cloned())
            .collect();
        let matched = out.iter().filter(|x| x.is_some()).count();
        log::info!(
            "GFF: matched {}/{} row names to gene loci",
            matched,
            row_names.len()
        );
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use genomic_data::gff::{FeatureType, GeneType, GffRecord, TranscriptId};
    use genomic_data::sam::Strand;

    fn rec(id: &str, sym: &str, chr: &str, start: i64, stop: i64, strand: Strand) -> GffRecord {
        GffRecord {
            seqname: chr.into(),
            feature_type: FeatureType::Gene,
            start,
            stop,
            strand,
            gene_id: GeneId::Ensembl(id.into()),
            gene_name: GeneSymbol::Symbol(sym.into()),
            gene_type: GeneType::CodingGene,
            transcript_id: TranscriptId::Missing,
        }
    }

    fn index() -> GeneLocusIndex {
        let records = vec![
            rec(
                "ENSG1",
                "GENE1",
                "chr17",
                7_661_779,
                7_687_538,
                Strand::Backward,
            ),
            rec(
                "ENSG2",
                "GENE2",
                "chr8",
                127_735_434,
                127_742_951,
                Strand::Forward,
            ),
        ];
        let map = GffRecordMap::from_map(
            genomic_data::gff::build_gene_map(&records, Some(&FeatureType::Gene)).unwrap(),
        );
        GeneLocusIndex::from_record_map(&map)
    }

    #[test]
    fn resolves_all_name_shapes() {
        let idx = index();
        assert_eq!(idx.resolve("GENE1").unwrap().gene_id.as_ref(), "ENSG1");
        assert_eq!(idx.resolve("ENSG1").unwrap().symbol.as_ref(), "GENE1");
        assert_eq!(idx.resolve("ENSG1.12").unwrap().symbol.as_ref(), "GENE1");
        assert_eq!(idx.resolve("ENSG1_GENE1").unwrap().symbol.as_ref(), "GENE1");
        assert_eq!(
            idx.resolve("ENSGX_GENE2/count/spliced")
                .unwrap()
                .gene_id
                .as_ref(),
            "ENSG2"
        );
        // Canonical matcher: case-insensitive, with or without the ENSG id.
        assert_eq!(idx.resolve("gene1").unwrap().gene_id.as_ref(), "ENSG1");
        assert_eq!(
            idx.resolve("ENSG1.12_GENE1").unwrap().symbol.as_ref(),
            "GENE1"
        );
        assert!(idx.resolve("NOPE").is_none());
    }

    fn names(xs: &[&str]) -> Vec<Box<str>> {
        xs.iter().map(|&x| x.into()).collect()
    }

    #[test]
    fn an_all_interval_axis_places_rows_by_name() {
        let loci =
            interval_loci(&names(&["chr2:1000-2000/depth", "chrUn_KI1:0-10", "2_5_9"])).unwrap();
        let l = &loci[0];
        assert_eq!(
            (l.chromosome.as_ref(), l.start, l.stop, l.tss),
            ("chr2", 1001, 2000, 1500)
        );
        assert_eq!(l.interval_name().as_ref(), "chr2:1000-2000");
        assert_eq!(l.gene_key().as_ref(), "chr2:1000-2000");
        assert_eq!(loci[1].chromosome.as_ref(), "chrUn_KI1");
        assert_eq!((loci[2].chromosome.as_ref(), loci[2].start), ("2", 6));
    }

    #[test]
    fn one_gene_row_makes_a_gene_axis() {
        for other in [
            "GENE1",
            "ENSG1_GENE1",
            "GENE-1",
            "chr2:2000-1000",
            "chr2:15",
        ] {
            assert!(
                interval_loci(&names(&["chr1:0-10", other])).is_none(),
                "{other}"
            );
        }
    }

    #[test]
    fn the_midpoint_does_not_overflow() {
        let l = &interval_loci(&names(&["chr1:0-9223372036854775807"])).unwrap()[0];
        assert!(l.tss > 0);
    }

    #[test]
    fn bins_on_different_grids_are_refused() {
        let same = interval_loci(&names(&["chr1:0-10", "chr1:10-20", "chr1:0-10"])).unwrap();
        assert!(check_one_grid(&same).is_ok());
        let shifted = interval_loci(&names(&["chr1:0-10", "chr1:5-15"])).unwrap();
        let err = check_one_grid(&shifted).unwrap_err().to_string();
        assert!(
            err.contains("chr1:0-10") && err.contains("chr1:5-15"),
            "{err}"
        );
    }

    #[test]
    fn tss_follows_strand_and_interval_is_bed() {
        let idx = index();
        let g1 = idx.resolve("GENE1").unwrap();
        assert_eq!(g1.tss, 7_687_538);
        assert_eq!(g1.interval_name().as_ref(), "chr17:7661778-7687538");
        let g2 = idx.resolve("GENE2").unwrap();
        assert_eq!(g2.tss, 127_735_434);
        assert_eq!(g2.gene_key().as_ref(), "ENSG2_GENE2");
    }
}
