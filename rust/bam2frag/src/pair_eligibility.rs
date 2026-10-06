//! Shared Phase-1 paired-fragment eligibility and candidate construction.
//!
//! Genomic coordinates are zero-based. Intervals are half-open except the two
//! named 5-prime endpoints, which are genomic base coordinates. The molecular
//! span is inclusive: `right_5p_genomic - left_5p_genomic + 1`.

use crate::{
    config::{ResolvedConfig, MOLECULAR_WINDOW_BP},
    reference::ReferenceCache,
};
use anyhow::Result;
use blake2::{Blake2b512, Digest};
use rust_htslib::bam::{record::Cigar, Record};
use serde::Serialize;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairOrientation {
    R1ForwardR2Reverse,
    R1ReverseR2Forward,
}

#[derive(Clone, Debug)]
pub struct PairCandidate {
    pub read1: Record,
    pub read2: Record,
    pub chrom: String,
    pub left_5p_genomic: i64,
    pub right_5p_genomic: i64,
    pub fragment_span_bp: u16,
    pub window_start0: i64,
    pub window_end0: i64,
    pub orientation: PairOrientation,
    pub fragment_hash: u64,
    /// Read-versus-reference tallies measured while validating this pair.
    ///
    /// Sequence disagreement is no longer grounds for rejection, so these are
    /// descriptive rather than justifying: they let a pilot quantify how much
    /// disagreement the reservoir actually admits, and at which spans and base
    /// qualities, without a second pass over the BAM.
    pub observation: FragmentObservation,
}

/// Number of validated-span bins used by the observational histograms.
pub const SPAN_BIN_COUNT: usize = 5;

/// Number of base-quality bins used by the observational histograms.
///
/// The edges are chosen for four-bin Illumina quality binning: a no-call bin, a
/// low bin, a middle bin, and the high bin.
pub const QUALITY_BIN_COUNT: usize = 4;

/// Number of per-fragment reference-mismatch-fraction histogram buckets.
///
/// Buckets 0..=9 are one percentage point wide; bucket 10 collects everything
/// at or above ten percent.
pub const MISMATCH_FRACTION_BIN_COUNT: usize = 11;

/// Map a validated inclusive molecular span to its observational bin.
pub fn span_bin(span_bp: i64) -> usize {
    match span_bp {
        ..=99 => 0,
        100..=149 => 1,
        150..=199 => 2,
        200..=249 => 3,
        _ => 4,
    }
}

/// Map a Phred base quality to its observational bin.
pub fn quality_bin(quality: u8) -> usize {
    match quality {
        0..=2 => 0,
        3..=22 => 1,
        23..=35 => 2,
        _ => 3,
    }
}

/// Per-fragment read-versus-reference tallies over the validated span.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct FragmentObservation {
    /// Validated inclusive molecular span, in base pairs.
    pub span_bp: i64,
    /// Span positions carrying at least one projected query base.
    pub covered_positions: u32,
    /// Positions where exactly one mate disagrees with the reference, either
    /// because only that mate covers the position or because the mates differ.
    pub unilateral_disagreements: u32,
    /// Positions where both mates disagree with the reference but not with the
    /// same allele.
    pub discordant_disagreements: u32,
    /// Positions where both mates agree with each other on a non-reference
    /// A/C/G/T allele.
    pub concordant_mismatches: u32,
    /// Disagreeing observations by [`quality_bin`] of the disagreeing base.
    pub disagreements_by_quality_bin: [u32; QUALITY_BIN_COUNT],
}

impl FragmentObservation {
    /// Total span positions where any mate disagrees with the reference.
    pub fn mismatching_positions(&self) -> u32 {
        self.unilateral_disagreements + self.discordant_disagreements + self.concordant_mismatches
    }

    /// Fraction of the validated span that disagrees with the reference.
    ///
    /// This is the quantity gated by `max_reference_mismatch_fraction`. A zero
    /// span yields zero rather than a division by zero; `construct_pair_candidate`
    /// has already rejected such pairs as [`PairRejection::Invalid5pSpan`].
    pub fn mismatch_fraction(&self) -> f64 {
        if self.span_bp <= 0 {
            return 0.0;
        }
        f64::from(self.mismatching_positions()) / self.span_bp as f64
    }

    /// Histogram bucket for this fragment's reference-mismatch fraction.
    pub fn mismatch_fraction_bin(&self) -> usize {
        let scaled = (self.mismatch_fraction() * 100.0).floor();
        if scaled >= (MISMATCH_FRACTION_BIN_COUNT - 1) as f64 {
            MISMATCH_FRACTION_BIN_COUNT - 1
        } else {
            scaled.max(0.0) as usize
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairRejection {
    NonInwardOrientation,
    Ambiguous5p,
    Clipping,
    Protruding3p,
    UnpairedOrMateUnmapped,
    ImproperPair,
    InvalidFlags,
    DifferentChromosomes,
    Mapq,
    FragmentSpan,
    IncompleteFragmentCoverage,
    /// The validated span disagrees with hg38 at more than
    /// `max_reference_mismatch_fraction` of its positions.
    ///
    /// This is a structural guard against mismapping, chimeras and unmasked
    /// repeats, not a sequence-quality judgement: a real cfDNA molecule sits
    /// far below the threshold, so only fragments that cannot plausibly come
    /// from this locus are rejected.
    ReferenceMismatchBurden,
    IndelCigar,
    Invalid5pSpan,
    IncompleteReferenceWindow,
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct PairRejectionCounters {
    pub pairs_accepted: u64,
    pub pairs_rejected_non_inward_orientation: u64,
    pub pairs_rejected_ambiguous_5p: u64,
    pub pairs_rejected_5p_clipping: u64,
    pub pairs_rejected_protruding_3p: u64,
    pub pairs_rejected_unpaired_or_mate_unmapped: u64,
    pub pairs_rejected_improper_pair: u64,
    pub pairs_rejected_invalid_flags: u64,
    pub pairs_rejected_different_chromosomes: u64,
    pub pairs_rejected_mapq: u64,
    pub pairs_rejected_fragment_span: u64,
    pub pairs_rejected_incomplete_fragment_coverage: u64,
    pub pairs_rejected_reference_mismatch_burden: u64,
    pub pairs_rejected_indel_cigar: u64,
    pub pairs_rejected_invalid_5p_span: u64,
    pub pairs_rejected_incomplete_reference_window: u64,
    /// Span positions of accepted pairs carrying at least one projected base.
    ///
    /// The four `accepted_*` counters below are observational, not rejections.
    /// They describe the disagreement the reservoir now admits, and their
    /// denominator is `accepted_covered_positions`.
    pub accepted_covered_positions: u64,
    pub accepted_unilateral_disagreements: u64,
    pub accepted_discordant_disagreements: u64,
    pub accepted_concordant_mismatches: u64,
    /// Accepted pairs where one of the eight paired 5-prime end-motif bases is
    /// not an unambiguous A/C/G/T.
    ///
    /// This is **not** a rejection, and not a quality judgement: an `N` is
    /// simply not a base, so a reader has nothing to reconstruct there.
    /// The fragment stays in the natural reservoir and remains valid for
    /// ordinary span reconstruction. Fragments whose motif bases *are* all
    /// A/C/G/T contribute their lowest motif quality to
    /// `FragmentObservationStats::end_motif_min_quality_histogram` instead, so
    /// the end-motif target threshold can be chosen downstream from the data
    /// rather than fixed at extraction.
    pub accepted_without_end_motif_bases: u64,
}

/// Stable serialized field order for [`PairRejectionCounters`].
///
/// `extraction_stats.csv` needs a fixed header and the `csv` crate cannot
/// flatten a nested struct, so the column order lives here and is asserted
/// against the derived JSON field names by `tests/extraction_stats.rs`.
pub const PAIR_COUNTER_FIELDS: [&str; 21] = [
    "pairs_accepted",
    "pairs_rejected_non_inward_orientation",
    "pairs_rejected_ambiguous_5p",
    "pairs_rejected_5p_clipping",
    "pairs_rejected_protruding_3p",
    "pairs_rejected_unpaired_or_mate_unmapped",
    "pairs_rejected_improper_pair",
    "pairs_rejected_invalid_flags",
    "pairs_rejected_different_chromosomes",
    "pairs_rejected_mapq",
    "pairs_rejected_fragment_span",
    "pairs_rejected_incomplete_fragment_coverage",
    "pairs_rejected_reference_mismatch_burden",
    "pairs_rejected_indel_cigar",
    "pairs_rejected_invalid_5p_span",
    "pairs_rejected_incomplete_reference_window",
    "accepted_covered_positions",
    "accepted_unilateral_disagreements",
    "accepted_discordant_disagreements",
    "accepted_concordant_mismatches",
    "accepted_without_end_motif_bases",
];

impl PairRejectionCounters {
    /// Accumulate another worker's counters, field by field.
    ///
    /// Both extraction commands fan out over chromosomes, so per-worker
    /// counters must be summed before they are reported or written.
    pub fn merge(&mut self, other: &Self) {
        self.pairs_accepted += other.pairs_accepted;
        self.pairs_rejected_non_inward_orientation += other.pairs_rejected_non_inward_orientation;
        self.pairs_rejected_ambiguous_5p += other.pairs_rejected_ambiguous_5p;
        self.pairs_rejected_5p_clipping += other.pairs_rejected_5p_clipping;
        self.pairs_rejected_protruding_3p += other.pairs_rejected_protruding_3p;
        self.pairs_rejected_unpaired_or_mate_unmapped +=
            other.pairs_rejected_unpaired_or_mate_unmapped;
        self.pairs_rejected_improper_pair += other.pairs_rejected_improper_pair;
        self.pairs_rejected_invalid_flags += other.pairs_rejected_invalid_flags;
        self.pairs_rejected_different_chromosomes += other.pairs_rejected_different_chromosomes;
        self.pairs_rejected_mapq += other.pairs_rejected_mapq;
        self.pairs_rejected_fragment_span += other.pairs_rejected_fragment_span;
        self.pairs_rejected_incomplete_fragment_coverage +=
            other.pairs_rejected_incomplete_fragment_coverage;
        self.pairs_rejected_reference_mismatch_burden +=
            other.pairs_rejected_reference_mismatch_burden;
        self.pairs_rejected_indel_cigar += other.pairs_rejected_indel_cigar;
        self.pairs_rejected_invalid_5p_span += other.pairs_rejected_invalid_5p_span;
        self.pairs_rejected_incomplete_reference_window +=
            other.pairs_rejected_incomplete_reference_window;
        self.accepted_covered_positions += other.accepted_covered_positions;
        self.accepted_unilateral_disagreements += other.accepted_unilateral_disagreements;
        self.accepted_discordant_disagreements += other.accepted_discordant_disagreements;
        self.accepted_concordant_mismatches += other.accepted_concordant_mismatches;
        self.accepted_without_end_motif_bases += other.accepted_without_end_motif_bases;
    }

    /// Return every counter in [`PAIR_COUNTER_FIELDS`] order.
    pub fn as_named_counts(&self) -> [(&'static str, u64); 21] {
        [
            ("pairs_accepted", self.pairs_accepted),
            (
                "pairs_rejected_non_inward_orientation",
                self.pairs_rejected_non_inward_orientation,
            ),
            (
                "pairs_rejected_ambiguous_5p",
                self.pairs_rejected_ambiguous_5p,
            ),
            (
                "pairs_rejected_5p_clipping",
                self.pairs_rejected_5p_clipping,
            ),
            (
                "pairs_rejected_protruding_3p",
                self.pairs_rejected_protruding_3p,
            ),
            (
                "pairs_rejected_unpaired_or_mate_unmapped",
                self.pairs_rejected_unpaired_or_mate_unmapped,
            ),
            (
                "pairs_rejected_improper_pair",
                self.pairs_rejected_improper_pair,
            ),
            (
                "pairs_rejected_invalid_flags",
                self.pairs_rejected_invalid_flags,
            ),
            (
                "pairs_rejected_different_chromosomes",
                self.pairs_rejected_different_chromosomes,
            ),
            ("pairs_rejected_mapq", self.pairs_rejected_mapq),
            (
                "pairs_rejected_fragment_span",
                self.pairs_rejected_fragment_span,
            ),
            (
                "pairs_rejected_incomplete_fragment_coverage",
                self.pairs_rejected_incomplete_fragment_coverage,
            ),
            (
                "pairs_rejected_reference_mismatch_burden",
                self.pairs_rejected_reference_mismatch_burden,
            ),
            (
                "pairs_rejected_indel_cigar",
                self.pairs_rejected_indel_cigar,
            ),
            (
                "pairs_rejected_invalid_5p_span",
                self.pairs_rejected_invalid_5p_span,
            ),
            (
                "pairs_rejected_incomplete_reference_window",
                self.pairs_rejected_incomplete_reference_window,
            ),
            (
                "accepted_covered_positions",
                self.accepted_covered_positions,
            ),
            (
                "accepted_unilateral_disagreements",
                self.accepted_unilateral_disagreements,
            ),
            (
                "accepted_discordant_disagreements",
                self.accepted_discordant_disagreements,
            ),
            (
                "accepted_concordant_mismatches",
                self.accepted_concordant_mismatches,
            ),
            (
                "accepted_without_end_motif_bases",
                self.accepted_without_end_motif_bases,
            ),
        ]
    }

    pub fn record(&mut self, result: &std::result::Result<PairCandidate, PairRejection>) {
        match result {
            Ok(pair) => {
                self.pairs_accepted += 1;
                let observation = &pair.observation;
                self.accepted_covered_positions += u64::from(observation.covered_positions);
                self.accepted_unilateral_disagreements +=
                    u64::from(observation.unilateral_disagreements);
                self.accepted_discordant_disagreements +=
                    u64::from(observation.discordant_disagreements);
                self.accepted_concordant_mismatches += u64::from(observation.concordant_mismatches);
            }
            Err(reason) => match reason {
                PairRejection::NonInwardOrientation => {
                    self.pairs_rejected_non_inward_orientation += 1
                }
                PairRejection::Ambiguous5p => self.pairs_rejected_ambiguous_5p += 1,
                PairRejection::Clipping => self.pairs_rejected_5p_clipping += 1,
                PairRejection::Protruding3p => self.pairs_rejected_protruding_3p += 1,
                PairRejection::UnpairedOrMateUnmapped => {
                    self.pairs_rejected_unpaired_or_mate_unmapped += 1
                }
                PairRejection::ImproperPair => self.pairs_rejected_improper_pair += 1,
                PairRejection::InvalidFlags => self.pairs_rejected_invalid_flags += 1,
                PairRejection::DifferentChromosomes => {
                    self.pairs_rejected_different_chromosomes += 1
                }
                PairRejection::Mapq => self.pairs_rejected_mapq += 1,
                PairRejection::FragmentSpan => self.pairs_rejected_fragment_span += 1,
                PairRejection::IncompleteFragmentCoverage => {
                    self.pairs_rejected_incomplete_fragment_coverage += 1
                }
                PairRejection::ReferenceMismatchBurden => {
                    self.pairs_rejected_reference_mismatch_burden += 1
                }
                PairRejection::IndelCigar => self.pairs_rejected_indel_cigar += 1,
                PairRejection::Invalid5pSpan => self.pairs_rejected_invalid_5p_span += 1,
                PairRejection::IncompleteReferenceWindow => {
                    self.pairs_rejected_incomplete_reference_window += 1
                }
            },
        }
    }
}

/// Span- and quality-stratified accounting for accepted fragments.
///
/// The scalar totals in [`PairRejectionCounters`] answer "how much disagreement
/// is admitted"; these arrays answer "is it concentrated anywhere". The span
/// arrays exist because a filter whose acceptance rate varies with fragment
/// length biases the stored length distribution, which is the one signal the
/// downstream length baselines depend on. They are persisted as JSON only,
/// because `extraction_stats.csv` carries one scalar per column.
#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct FragmentObservationStats {
    pub accepted_by_span_bin: [u64; SPAN_BIN_COUNT],
    pub covered_positions_by_span_bin: [u64; SPAN_BIN_COUNT],
    pub disagreements_by_span_bin: [u64; SPAN_BIN_COUNT],
    pub concordant_mismatches_by_span_bin: [u64; SPAN_BIN_COUNT],
    pub disagreements_by_quality_bin: [u64; QUALITY_BIN_COUNT],
    pub mismatch_fraction_histogram: [u64; MISMATCH_FRACTION_BIN_COUNT],
    /// Lowest base quality across the eight paired 5-prime end-motif bases, by
    /// [`quality_bin`], for accepted pairs whose motif bases are all A/C/G/T.
    ///
    /// Extraction sets no end-motif quality threshold. This histogram is what
    /// a reader uses to choose one, so the decision is made from the
    /// observed distribution rather than fixed before the data is seen.
    pub end_motif_min_quality_histogram: [u64; QUALITY_BIN_COUNT],
}

impl FragmentObservationStats {
    /// Accumulate another worker's histograms, bin by bin.
    pub fn merge(&mut self, other: &Self) {
        for index in 0..SPAN_BIN_COUNT {
            self.accepted_by_span_bin[index] += other.accepted_by_span_bin[index];
            self.covered_positions_by_span_bin[index] += other.covered_positions_by_span_bin[index];
            self.disagreements_by_span_bin[index] += other.disagreements_by_span_bin[index];
            self.concordant_mismatches_by_span_bin[index] +=
                other.concordant_mismatches_by_span_bin[index];
        }
        for index in 0..QUALITY_BIN_COUNT {
            self.disagreements_by_quality_bin[index] += other.disagreements_by_quality_bin[index];
            self.end_motif_min_quality_histogram[index] +=
                other.end_motif_min_quality_histogram[index];
        }
        for index in 0..MISMATCH_FRACTION_BIN_COUNT {
            self.mismatch_fraction_histogram[index] += other.mismatch_fraction_histogram[index];
        }
    }

    /// Fold one accepted fragment's observation into the histograms.
    pub fn record_accepted(&mut self, observation: &FragmentObservation) {
        let bin = span_bin(observation.span_bp);
        self.accepted_by_span_bin[bin] += 1;
        self.covered_positions_by_span_bin[bin] += u64::from(observation.covered_positions);
        self.disagreements_by_span_bin[bin] += u64::from(observation.mismatching_positions());
        self.concordant_mismatches_by_span_bin[bin] += u64::from(observation.concordant_mismatches);
        for index in 0..QUALITY_BIN_COUNT {
            self.disagreements_by_quality_bin[index] +=
                u64::from(observation.disagreements_by_quality_bin[index]);
        }
        self.mismatch_fraction_histogram[observation.mismatch_fraction_bin()] += 1;
    }

    /// Fold one accepted fragment's lowest end-motif base quality into the
    /// histogram. Callers skip this when a motif base is not A/C/G/T.
    pub fn record_end_motif_quality(&mut self, quality: u8) {
        self.end_motif_min_quality_histogram[quality_bin(quality)] += 1;
    }
}

#[derive(Clone, Copy)]
struct Observation {
    base: u8,
    quality: u8,
}

struct Projection {
    observations: Vec<(i64, Observation)>,
    first_aligned0: i64,
    last_aligned0: i64,
}

fn project(record: &Record) -> std::result::Result<Projection, PairRejection> {
    let sequence = record.seq().as_bytes();
    if sequence.len() != record.qual().len() || sequence.is_empty() {
        return Err(PairRejection::Ambiguous5p);
    }
    let mut reference_pos0 = record.pos();
    let mut query_idx = 0usize;
    let mut observations = Vec::with_capacity(sequence.len());
    for cigar in record.cigar().iter() {
        match cigar {
            Cigar::Ins(length) | Cigar::Del(length) if *length > 0 => {
                return Err(PairRejection::IndelCigar)
            }
            Cigar::SoftClip(length) | Cigar::HardClip(length) if *length > 0 => {
                return Err(PairRejection::Clipping)
            }
            Cigar::Match(length) | Cigar::Equal(length) | Cigar::Diff(length) => {
                for _ in 0..*length {
                    let Some((&base, &quality)) =
                        sequence.get(query_idx).zip(record.qual().get(query_idx))
                    else {
                        return Err(PairRejection::Ambiguous5p);
                    };
                    observations.push((
                        reference_pos0,
                        Observation {
                            base: base.to_ascii_uppercase(),
                            quality,
                        },
                    ));
                    reference_pos0 += 1;
                    query_idx += 1;
                }
            }
            Cigar::RefSkip(length) => reference_pos0 += *length as i64,
            Cigar::Pad(_) | Cigar::HardClip(_) | Cigar::SoftClip(_) => {}
            Cigar::Ins(_) | Cigar::Del(_) => {}
        }
    }
    if query_idx != sequence.len() || observations.is_empty() {
        return Err(PairRejection::Ambiguous5p);
    }
    Ok(Projection {
        first_aligned0: observations.first().unwrap().0,
        last_aligned0: observations.last().unwrap().0,
        observations,
    })
}

fn stable_fragment_hash(chrom: &str, left_5p: i64, right_5p: i64, qname: &[u8]) -> u64 {
    let mut hasher = Blake2b512::new();
    hasher.update(chrom.as_bytes());
    hasher.update(left_5p.to_le_bytes());
    hasher.update(right_5p.to_le_bytes());
    hasher.update(qname);
    let digest = hasher.finalize();
    u64::from_le_bytes(digest[..8].try_into().unwrap())
}

/// Validate a BAM pair and construct an event-free molecular candidate.
pub fn construct_pair_candidate(
    first: &Record,
    second: &Record,
    chrom: &str,
    chrom_len_bp: i64,
    cfg: &ResolvedConfig,
    reference: &mut ReferenceCache,
) -> Result<std::result::Result<PairCandidate, PairRejection>> {
    let (read1, read2) = if first.is_first_in_template() && second.is_last_in_template() {
        (first, second)
    } else if second.is_first_in_template() && first.is_last_in_template() {
        (second, first)
    } else {
        return Ok(Err(PairRejection::InvalidFlags));
    };

    if read1.qname() != read2.qname()
        || read1.is_secondary()
        || read2.is_secondary()
        || read1.is_supplementary()
        || read2.is_supplementary()
        || read1.is_duplicate()
        || read2.is_duplicate()
        || read1.is_quality_check_failed()
        || read2.is_quality_check_failed()
    {
        return Ok(Err(PairRejection::InvalidFlags));
    }
    if !read1.is_paired()
        || !read2.is_paired()
        || read1.is_unmapped()
        || read2.is_unmapped()
        || read1.is_mate_unmapped()
        || read2.is_mate_unmapped()
    {
        return Ok(Err(PairRejection::UnpairedOrMateUnmapped));
    }
    if !read1.is_proper_pair() || !read2.is_proper_pair() {
        return Ok(Err(PairRejection::ImproperPair));
    }
    if read1.tid() < 0
        || read1.tid() != read2.tid()
        || read1.mtid() != read2.tid()
        || read2.mtid() != read1.tid()
    {
        return Ok(Err(PairRejection::DifferentChromosomes));
    }
    if read1.is_mate_reverse() != read2.is_reverse()
        || read2.is_mate_reverse() != read1.is_reverse()
    {
        return Ok(Err(PairRejection::InvalidFlags));
    }
    if read1.mapq() < cfg.min_mapq || read2.mapq() < cfg.min_mapq {
        return Ok(Err(PairRejection::Mapq));
    }
    if read1.is_reverse() == read2.is_reverse() {
        return Ok(Err(PairRejection::NonInwardOrientation));
    }

    let projection1 = match project(read1) {
        Ok(value) => value,
        Err(reason) => return Ok(Err(reason)),
    };
    let projection2 = match project(read2) {
        Ok(value) => value,
        Err(reason) => return Ok(Err(reason)),
    };
    let (forward_projection, reverse_projection) = if read1.is_reverse() {
        (&projection2, &projection1)
    } else {
        (&projection1, &projection2)
    };
    let left_5p = forward_projection.first_aligned0;
    let right_5p = reverse_projection.last_aligned0;
    if left_5p < 0 || right_5p < left_5p {
        return Ok(Err(PairRejection::Invalid5pSpan));
    }
    if forward_projection.last_aligned0 > right_5p || reverse_projection.first_aligned0 < left_5p {
        return Ok(Err(PairRejection::Protruding3p));
    }
    let span = right_5p - left_5p + 1;
    if span <= 0 {
        return Ok(Err(PairRejection::Invalid5pSpan));
    }
    if span < cfg.min_len || span > cfg.max_len {
        return Ok(Err(PairRejection::FragmentSpan));
    }

    let window_width_bp = MOLECULAR_WINDOW_BP as i64;
    let left_context_bp = (window_width_bp - span) / 2;
    let window_start0 = left_5p - left_context_bp;
    let window_end0 = window_start0 + window_width_bp;
    if window_start0 < 0 || window_end0 > chrom_len_bp {
        return Ok(Err(PairRejection::IncompleteReferenceWindow));
    }

    let reference_window = reference.fetch_span_exact(chrom, window_start0, window_end0)?;
    let molecular_offset = (left_5p - window_start0) as usize;
    let reference_span = &reference_window[molecular_offset..molecular_offset + span as usize];
    let mut coverage = vec![[None; 2]; span as usize];
    for (mate_idx, projection) in [forward_projection, reverse_projection]
        .into_iter()
        .enumerate()
    {
        for (position0, observation) in &projection.observations {
            if *position0 >= left_5p && *position0 <= right_5p {
                coverage[(*position0 - left_5p) as usize][mate_idx] = Some(*observation);
            }
        }
    }

    // Sequence disagreement is recorded, never rejected. Both read planes and
    // both quality planes are stored, so the model sees each mate's observation
    // together with the confidence the instrument assigned it, and can weigh
    // them itself. Discarding a whole molecule over one base would delete its
    // length and both end motifs, and would do so more often the longer the
    // fragment, biasing the stored length distribution.
    let mut observation = FragmentObservation {
        span_bp: span,
        ..FragmentObservation::default()
    };
    for (offset, observations) in coverage.iter().enumerate() {
        let reference_base = reference_span[offset].to_ascii_uppercase();
        match observations {
            // A position no mate ever reached is uninterpretable, not merely
            // uncertain, so it still rejects the pair.
            [None, None] => return Ok(Err(PairRejection::IncompleteFragmentCoverage)),
            [Some(single), None] | [None, Some(single)] => {
                observation.covered_positions += 1;
                if single.base != reference_base {
                    observation.unilateral_disagreements += 1;
                    observation.disagreements_by_quality_bin[quality_bin(single.quality)] += 1;
                }
            }
            [Some(left), Some(right)] => {
                observation.covered_positions += 1;
                let left_mismatch = left.base != reference_base;
                let right_mismatch = right.base != reference_base;
                match (left_mismatch, right_mismatch) {
                    (false, false) => {}
                    // Exactly one mate disagrees with the reference: the other
                    // corroborates it, so this is almost always a base-calling
                    // error in the disagreeing mate.
                    (true, false) => {
                        observation.unilateral_disagreements += 1;
                        observation.disagreements_by_quality_bin[quality_bin(left.quality)] += 1;
                    }
                    (false, true) => {
                        observation.unilateral_disagreements += 1;
                        observation.disagreements_by_quality_bin[quality_bin(right.quality)] += 1;
                    }
                    (true, true) => {
                        if left.base == right.base && b"ACGT".contains(&left.base) {
                            observation.concordant_mismatches += 1;
                        } else {
                            observation.discordant_disagreements += 1;
                        }
                        observation.disagreements_by_quality_bin[quality_bin(left.quality)] += 1;
                    }
                }
            }
        }
    }

    // Structural guard against mismapping, chimeras and unmasked repeats. A
    // molecule genuinely from this locus disagrees with hg38 at well under one
    // percent of its positions; anything above the configured fraction cannot
    // be interpreted against this reference window at all.
    if observation.mismatch_fraction() > cfg.max_reference_mismatch_fraction {
        return Ok(Err(PairRejection::ReferenceMismatchBurden));
    }

    let orientation = if read1.is_reverse() {
        PairOrientation::R1ReverseR2Forward
    } else {
        PairOrientation::R1ForwardR2Reverse
    };
    Ok(Ok(PairCandidate {
        read1: read1.clone(),
        read2: read2.clone(),
        chrom: chrom.to_string(),
        left_5p_genomic: left_5p,
        right_5p_genomic: right_5p,
        fragment_span_bp: span as u16,
        window_start0,
        window_end0,
        orientation,
        fragment_hash: stable_fragment_hash(chrom, left_5p, right_5p, read1.qname()),
        observation,
    }))
}
