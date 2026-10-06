use crate::pair_eligibility::{FragmentObservationStats, PairCandidate, PairRejectionCounters};
use serde::Serialize;

/// Deterministic identity of one natural reservoir bucket.
///
/// The bank stores exactly one payload kind, `natural`, so there is no
/// match/mismatch stream and no event category: the key is purely positional.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct ReservoirKey {
    pub region_id: i64,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct ChromStats {
    pub chrom: String,
    pub reads: u64,
    pub skipped: u64,
    pub evicted: u64,
    pub invalid_fragments: u64,
    /// Per-reason pair-eligibility accounting for this chromosome, persisted so
    /// a pilot does not have to scrape the extraction log.
    pub pair_rejections: PairRejectionCounters,
    /// Span- and quality-stratified read-versus-reference histograms for the
    /// fragments this chromosome accepted. Persisted to JSON only, because the
    /// CSV carries one scalar per column.
    pub observations: FragmentObservationStats,
}

impl ChromStats {
    pub fn new(chrom: impl Into<String>) -> Self {
        Self {
            chrom: chrom.into(),
            ..Default::default()
        }
    }

    /// Accumulate another chromosome's counts into a running total.
    pub fn merge(&mut self, other: &Self) {
        self.reads += other.reads;
        self.skipped += other.skipped;
        self.evicted += other.evicted;
        self.invalid_fragments += other.invalid_fragments;
        self.pair_rejections.merge(&other.pair_rejections);
        self.observations.merge(&other.observations);
    }
}

/// Fully materialized fragment ready for serialization into the bank.
///
/// Deliberately contains only metadata required after materialization. In
/// particular, it must not retain the source `PairCandidate`, because that
/// would keep both heap-backed BAM `Record`s alive for the entire materialize
/// stage and make memory usage grow with the number of completed fragments.
#[derive(Clone, Debug)]
pub struct MaterializedFragment {
    pub key: ReservoirKey,
    pub anchor_rpos0: i64,
    pub ref_aln: Vec<u8>,
    pub r1_aln: Vec<u8>,
    pub r2_aln: Vec<u8>,
    pub r1_qual: Vec<u8>,
    pub r2_qual: Vec<u8>,
    pub r1_cycle: Vec<u8>,
    pub r2_cycle: Vec<u8>,
    pub fragment_span_bp: u16,
    pub orientation_code: u8,
    pub atac_overlap_fractions: [f32; 3],
}

#[derive(Clone, Debug)]
pub struct ValidatedFragmentCandidate {
    pub key: ReservoirKey,
    pub pair: PairCandidate,
    pub anchor_rpos0: i64,
}
