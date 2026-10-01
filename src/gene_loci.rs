//! Resolve data-beans row names to genomic loci.
//!
//! A row named as a genomic interval (`chr:start-end`, as `faba read-depth`
//! names its bins, or `chr_start_end` once data-beans has aligned them;
//! 0-based half-open) is its own locus. Any other row is a
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
use genomic_data::coordinates::parse_peak_coordinates;
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

    /// `{gene_id}_{symbol}`, faba's gene key; an interval's own name.
    pub fn gene_key(&self) -> Box<str> {
        if self.symbol.is_empty() {
            return self.gene_id.clone();
        }
        format!("{}_{}", self.gene_id, self.symbol).into()
    }

    /// The locus a row named as a genomic interval stands for, placed at its
    /// midpoint: `chr:start-end` as `faba read-depth` writes it, or
    /// `chr_start_end` as data-beans aligns it (0-based half-open), with an
    /// optional faba `/modality/...` suffix. `None` when the name is not an
    /// interval.
    pub fn from_interval_name(row_name: &str) -> Option<Self> {
        let core: Box<str> = row_name.split('/').next()?.trim().into();
        let r = parse_peak_coordinates(std::slice::from_ref(&core)).pop()??;
        if r.start < 0 {
            return None;
        }
        Some(Self {
            tss: (r.start + 1 + r.end) / 2,
            start: r.start + 1,
            stop: r.end,
            gene_id: core,
            symbol: "".into(),
            chromosome: r.chr,
        })
    }
}

/// Every row's locus: interval rows stand for themselves, and the rest are
/// looked up in the GFF that `gff` loads, called only when some row is not an
/// interval (so a read-depth matrix needs no annotation).
pub fn resolve_rows(
    row_names: &[Box<str>],
    gff: impl FnOnce() -> anyhow::Result<GeneLocusIndex>,
) -> anyhow::Result<Vec<Option<GeneLocus>>> {
    let mut loci: Vec<Option<GeneLocus>> = row_names
        .par_iter()
        .map(|n| GeneLocus::from_interval_name(n))
        .collect();
    let n_genes = loci.iter().filter(|l| l.is_none()).count();
    log::info!(
        "{} of {} rows are genomic intervals",
        row_names.len() - n_genes,
        row_names.len()
    );
    if n_genes > 0 {
        let index = gff()?;
        loci.par_iter_mut()
            .zip(row_names.par_iter())
            .filter(|(l, _)| l.is_none())
            .for_each(|(l, n)| *l = index.resolve(n).cloned());
        let matched = loci.iter().filter(|l| l.is_some()).count();
        log::info!(
            "GFF: matched {matched}/{} row names to loci",
            row_names.len()
        );
    }
    Ok(loci)
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
    let core = row_name.split('/').next().unwrap_or(row_name).trim();
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

    #[test]
    fn interval_rows_are_their_own_loci() {
        let l = GeneLocus::from_interval_name("chr2:1000-2000/depth").unwrap();
        assert_eq!(
            (l.chromosome.as_ref(), l.start, l.stop),
            ("chr2", 1001, 2000)
        );
        assert_eq!(l.tss, 1500);
        assert_eq!(l.interval_name().as_ref(), "chr2:1000-2000");
        assert_eq!(l.gene_key().as_ref(), "chr2:1000-2000");
        for bad in [
            "GENE1",
            "ENSG1_GENE1",
            "chr2:2000-1000",
            "chr2:x-10",
            ":0-10",
        ] {
            assert!(GeneLocus::from_interval_name(bad).is_none(), "{bad}");
        }
    }

    #[test]
    fn the_gff_is_read_only_for_gene_rows() {
        let bins: Vec<Box<str>> = vec!["chr1:0-100".into(), "chr1:100-200".into()];
        let loci = resolve_rows(&bins, || anyhow::bail!("no GFF needed")).unwrap();
        assert!(loci.iter().all(Option::is_some));

        let mixed: Vec<Box<str>> = vec!["chr1:0-100".into(), "GENE1".into(), "NOPE".into()];
        let loci = resolve_rows(&mixed, || Ok(index())).unwrap();
        assert_eq!(loci[0].as_ref().unwrap().chromosome.as_ref(), "chr1");
        assert_eq!(loci[1].as_ref().unwrap().gene_id.as_ref(), "ENSG1");
        assert!(loci[2].is_none());
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
