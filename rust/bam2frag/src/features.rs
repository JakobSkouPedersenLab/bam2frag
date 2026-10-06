use crate::{
    atac::TissueAtacMasks,
    config::ResolvedConfig,
    model::{MaterializedFragment, ValidatedFragmentCandidate},
    pair_eligibility::{
        FragmentObservationStats, PairCandidate, PairRejection, PairRejectionCounters,
    },
    reference::ReferenceCache,
};
use anyhow::{bail, Result};
use rust_htslib::bam::{record::Cigar, Record};
use std::collections::HashMap;

fn ref_query_map(r: &Record, start: i64, end: i64) -> HashMap<i64, usize> {
    let seq_len = r.seq_len();
    let mut out = HashMap::new();
    let mut rp = r.pos();
    let mut qp = 0usize;
    for c in r.cigar().iter() {
        match c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for _ in 0..*n {
                    if rp >= start && rp < end {
                        out.insert(rp, qp);
                    }
                    rp += 1;
                    qp += 1
                }
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => qp += *n as usize,
            Cigar::Del(n) | Cigar::RefSkip(n) => rp += *n as i64,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    out.retain(|_, q| *q < seq_len);
    out
}

fn base(r: &Record, q: usize) -> u8 {
    r.seq()
        .as_bytes()
        .get(q)
        .copied()
        .unwrap_or(b'N')
        .to_ascii_uppercase()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CatalogAlleleSupport {
    Alt,
    Ref,
    Ambiguous,
}

/// Determine whether a paired fragment supports one catalogue SNV allele.
///
/// `pos0` is zero-based and REF/ALT use ordinary left-anchored VCF notation.
/// The tumour experiment is biallelic-SNV only, so any record whose REF or ALT
/// is not a single base is reported as `Ambiguous` and can never create a
/// group. Support requires both mates to cover the locus, to agree, and to
/// carry base quality of at least `min_bq`.
pub fn catalog_allele_support(
    a: &Record,
    b: &Record,
    pos0: i64,
    obs_ref: &str,
    obs_alt: &str,
    min_bq: u8,
) -> CatalogAlleleSupport {
    let rf = obs_ref.as_bytes();
    let al = obs_alt.as_bytes();
    if rf.len() != 1 || al.len() != 1 || rf == al {
        return CatalogAlleleSupport::Ambiguous;
    }

    let m1 = ref_query_map(a, pos0, pos0 + 1);
    let m2 = ref_query_map(b, pos0, pos0 + 1);
    let (Some(&q1), Some(&q2)) = (m1.get(&pos0), m2.get(&pos0)) else {
        return CatalogAlleleSupport::Ambiguous;
    };
    if a.qual().get(q1).copied().unwrap_or(0) < min_bq
        || b.qual().get(q2).copied().unwrap_or(0) < min_bq
    {
        return CatalogAlleleSupport::Ambiguous;
    }
    let b1 = base(a, q1);
    let b2 = base(b, q2);
    if b1 != b2 {
        return CatalogAlleleSupport::Ambiguous;
    }
    if b1 == al[0].to_ascii_uppercase() {
        return CatalogAlleleSupport::Alt;
    }
    if b1 == rf[0].to_ascii_uppercase() {
        return CatalogAlleleSupport::Ref;
    }
    CatalogAlleleSupport::Ambiguous
}

#[inline]
fn encode_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'N' => 0,
        b'A' => 1,
        b'C' => 2,
        b'G' => 3,
        b'T' => 4,
        _ => 0,
    }
}

fn build_read_track_into(
    r: &Record,
    ws: i64,
    w: usize,
    bases: &mut Vec<u8>,
    qualities: &mut Vec<u8>,
    cycles: &mut Vec<u8>,
) {
    bases.clear();
    bases.resize(w, 0);
    qualities.clear();
    qualities.resize(w, 0);
    cycles.clear();
    cycles.resize(w, 0);
    let seq = r.seq().as_bytes();
    let mut rp = r.pos();
    let mut qp = 0usize;
    for c in r.cigar().iter() {
        match c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                for _ in 0..*n {
                    let i = rp - ws;
                    if i >= 0 && (i as usize) < w {
                        bases[i as usize] =
                            encode_base(seq.get(qp).copied().unwrap_or(b'N').to_ascii_uppercase());
                        qualities[i as usize] = r.qual().get(qp).copied().unwrap_or(0);
                        let cycle0 = if r.is_reverse() {
                            seq.len() - 1 - qp
                        } else {
                            qp
                        };
                        cycles[i as usize] = (cycle0 + 1) as u8;
                    }
                    rp += 1;
                    qp += 1
                }
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => qp += *n as usize,
            Cigar::Del(n) => rp += *n as i64,
            Cigar::RefSkip(n) => rp += *n as i64,
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
}

#[derive(Debug, Default)]
pub struct MaterializeScratch {
    pub ref_aln: Vec<u8>,
    pub r1_aln: Vec<u8>,
    pub r2_aln: Vec<u8>,
    pub r1_qual: Vec<u8>,
    pub r2_qual: Vec<u8>,
    pub r1_cycle: Vec<u8>,
    pub r2_cycle: Vec<u8>,
    pub reference_span: Vec<u8>,
}

impl MaterializeScratch {
    pub fn with_width(width: usize) -> Self {
        Self {
            ref_aln: Vec::with_capacity(width),
            r1_aln: Vec::with_capacity(width),
            r2_aln: Vec::with_capacity(width),
            r1_qual: Vec::with_capacity(width),
            r2_qual: Vec::with_capacity(width),
            r1_cycle: Vec::with_capacity(width),
            r2_cycle: Vec::with_capacity(width),
            reference_span: Vec::with_capacity(width),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MaterializedScalars {
    pub fragment_span_bp: u16,
    pub orientation_code: u8,
    pub atac_overlap_fractions: [f32; 3],
}

fn encode_reference(reference: &[u8], output: &mut Vec<u8>) -> Result<()> {
    output.clear();
    output.reserve(reference.len());
    for base in reference {
        let base = base.to_ascii_uppercase();
        match base {
            b'A' | b'C' | b'G' | b'T' | b'N' => output.push(encode_base(base)),
            // Valid IUPAC ambiguity codes carry no unambiguous reference base
            // for the model, so represent them using the existing N code.
            b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => {
                output.push(encode_base(b'N'))
            }
            _ => bail!("reference window contains unsupported base byte {base}"),
        }
    }
    Ok(())
}

/// Materialize the Phase-1 validated molecular candidate without reconstructing
/// endpoints from legacy outer-CIGAR coordinates.
pub fn materialize_validated_into_scratch(
    candidate: ValidatedFragmentCandidate,
    cfg: &ResolvedConfig,
    reference: &mut ReferenceCache,
    scratch: &mut MaterializeScratch,
    atac_masks: Option<&TissueAtacMasks>,
) -> Result<MaterializedScalars> {
    let pair = candidate.pair;
    let width = cfg.typed.feature_sets.fragment_sequence.seq_window_w;
    scratch.reference_span =
        reference.fetch_span_exact(&pair.chrom, pair.window_start0, pair.window_end0)?;
    if scratch.reference_span.len() != width {
        bail!(
            "complete reference window has {} bp, expected {width}",
            scratch.reference_span.len()
        );
    }
    if pair.read1.seq().len() > u8::MAX as usize || pair.read2.seq().len() > u8::MAX as usize {
        bail!("u8 cycle planes support reads up to 255 cycles");
    }
    encode_reference(&scratch.reference_span, &mut scratch.ref_aln)?;
    build_read_track_into(
        &pair.read1,
        pair.window_start0,
        width,
        &mut scratch.r1_aln,
        &mut scratch.r1_qual,
        &mut scratch.r1_cycle,
    );
    build_read_track_into(
        &pair.read2,
        pair.window_start0,
        width,
        &mut scratch.r2_aln,
        &mut scratch.r2_qual,
        &mut scratch.r2_cycle,
    );
    let orientation_code = match pair.orientation {
        crate::pair_eligibility::PairOrientation::R1ForwardR2Reverse => 0,
        crate::pair_eligibility::PairOrientation::R1ReverseR2Forward => 1,
    };
    Ok(MaterializedScalars {
        fragment_span_bp: pair.fragment_span_bp,
        orientation_code,
        atac_overlap_fractions: atac_masks
            .map(|masks| masks.overlap_fractions(pair.left_5p_genomic, pair.right_5p_genomic + 1))
            .unwrap_or([0.0; 3]),
    })
}

/// Materialize one Phase-1 validated candidate into an owned packed fragment row.
///
/// This is the owned counterpart of [`materialize_validated_into_scratch`] and
/// exists so that non-streaming callers share the exact validated
/// materializer used by the natural extractor rather than reconstructing
/// endpoints themselves.
pub fn materialize_validated(
    candidate: ValidatedFragmentCandidate,
    cfg: &ResolvedConfig,
    reference: &mut ReferenceCache,
    atac_masks: Option<&TissueAtacMasks>,
) -> Result<MaterializedFragment> {
    let key = candidate.key.clone();
    let anchor_rpos0 = candidate.anchor_rpos0;
    let width = cfg.typed.feature_sets.fragment_sequence.seq_window_w;
    let mut scratch = MaterializeScratch::with_width(width);
    let scalars =
        materialize_validated_into_scratch(candidate, cfg, reference, &mut scratch, atac_masks)?;
    Ok(MaterializedFragment {
        key,
        anchor_rpos0,
        ref_aln: scratch.ref_aln,
        r1_aln: scratch.r1_aln,
        r2_aln: scratch.r2_aln,
        r1_qual: scratch.r1_qual,
        r2_qual: scratch.r2_qual,
        r1_cycle: scratch.r1_cycle,
        r2_cycle: scratch.r2_cycle,
        fragment_span_bp: scalars.fragment_span_bp,
        orientation_code: scalars.orientation_code,
        atac_overlap_fractions: scalars.atac_overlap_fractions,
    })
}

/// Decode one base code back to its IUPAC letter.
pub fn decode_base(code: u8) -> u8 {
    match code {
        1 => b'A',
        2 => b'C',
        3 => b'G',
        4 => b'T',
        _ => b'N',
    }
}

/// Return the first double-stranded, high-quality REF-versus-read difference
/// inside the validated molecular span of an already-materialized packed fragment row.
///
/// The natural reservoir admits unilateral and discordant disagreement, so this
/// function performs the consensus test itself rather than inheriting it from
/// pair eligibility: a difference qualifies only when both mates cover the
/// position, agree on the same non-reference base, and both carry base quality
/// of at least `min_bq`. That is the same locus-local standard
/// [`catalog_allele_support`] applies at a catalogue site, so the
/// `other_mismatch` control group is held to the bar as `tumor_ALT`.
///
/// Reading the stored planes keeps the tumour experiment's mismatch definition
/// identical to the bytes written into the fragment bank.
///
/// Returns `(window_index, reference_code, observed_code)`.
pub fn first_consensus_snv(fragment: &MaterializedFragment, min_bq: u8) -> Option<(usize, u8, u8)> {
    let width = fragment.ref_aln.len();
    let span = fragment.fragment_span_bp as usize;
    if span == 0 || span > width {
        return None;
    }
    let left_cut_idx = (width - span) / 2;
    for offset in left_cut_idx..left_cut_idx + span {
        let reference_code = fragment.ref_aln[offset];
        // Code zero means absent or an ambiguous base; neither is evidence, so
        // a position either mate fails to cover cannot carry a consensus SNV.
        let r1_code = fragment.r1_aln[offset];
        let r2_code = fragment.r2_aln[offset];
        if r1_code == 0 || r2_code == 0 {
            continue;
        }
        if r1_code != r2_code || r1_code == reference_code {
            continue;
        }
        if fragment.r1_qual[offset] < min_bq || fragment.r2_qual[offset] < min_bq {
            continue;
        }
        return Some((offset, reference_code, r1_code));
    }
    None
}

/// Number of sequenced cycles forming one 5-prime end motif.
const MOTIF_BASES_PER_END: usize = 4;

/// Return the lowest base quality across a validated pair's eight paired
/// 5-prime end-motif bases, or `None` when any of them is not an unambiguous
/// A/C/G/T.
///
/// Stored cycle `k` is query index `k-1` on a forward read and `len-k` on a
/// reverse read, matching `build_read_track_into`, so this reads the same eight
/// bases a downstream model's corruption masks.
///
/// Extraction applies no base-quality threshold, here or anywhere else: this
/// reports the quality so a reader can decide end-motif target
/// eligibility for itself. `None` is not a quality judgement — an `N` is simply
/// not a base, so there is nothing to reconstruct.
pub fn end_motif_min_quality(read1: &Record, read2: &Record) -> Option<u8> {
    let mut lowest = u8::MAX;
    for read in [read1, read2] {
        let sequence = read.seq().as_bytes();
        let qualities = read.qual();
        if sequence.len() < MOTIF_BASES_PER_END || qualities.len() != sequence.len() {
            return None;
        }
        for cycle0 in 0..MOTIF_BASES_PER_END {
            let query_idx = if read.is_reverse() {
                sequence.len() - 1 - cycle0
            } else {
                cycle0
            };
            if !(1..=4).contains(&encode_base(sequence[query_idx])) {
                return None;
            }
            lowest = lowest.min(qualities[query_idx]);
        }
    }
    Some(lowest)
}

/// Record one shared eligibility decision together with the accepted fragment's
/// observational histograms, for either extraction command's QC counters.
///
/// Nothing here rejects a fragment. Sequence disagreement and end-motif base
/// quality are both descriptive: an accepted pair's read-versus-reference
/// tallies and its lowest end-motif base quality are folded into `observations`
/// so a pilot can see the distributions the reservoir admits, and the Fragment
/// a reader chooses its own end-motif target threshold from the stored planes.
pub fn record_pair_decision(
    counters: &mut PairRejectionCounters,
    observations: &mut FragmentObservationStats,
    decision: &std::result::Result<PairCandidate, PairRejection>,
    read1: &Record,
    read2: &Record,
) {
    counters.record(decision);
    if let Ok(pair) = decision {
        observations.record_accepted(&pair.observation);
        match end_motif_min_quality(read1, read2) {
            Some(quality) => observations.record_end_motif_quality(quality),
            None => counters.accepted_without_end_motif_bases += 1,
        }
    }
}
