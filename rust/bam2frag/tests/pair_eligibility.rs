mod common;

use bam2frag::{
    config::MOLECULAR_WINDOW_BP,
    pair_eligibility::{construct_pair_candidate, PairRejection, PairRejectionCounters},
    reference::ReferenceCache,
};
use common::natural_resolved_config;
use rust_htslib::bam::{
    record::{Cigar, CigarString},
    Record,
};
use std::fs;
use tempfile::TempDir;

const CHROM_LEN: i64 = 2_000;

fn reference_cache() -> (TempDir, ReferenceCache) {
    let temp = tempfile::tempdir().unwrap();
    let fasta = temp.path().join("ref.fa");
    let sequence = "A".repeat(CHROM_LEN as usize);
    fs::write(&fasta, format!(">chr1\n{sequence}\n")).unwrap();
    fs::write(
        temp.path().join("ref.fa.fai"),
        format!(
            "chr1\t{}\t6\t{}\t{}\n",
            sequence.len(),
            sequence.len(),
            sequence.len() + 1
        ),
    )
    .unwrap();
    (temp, ReferenceCache::open(&fasta).unwrap())
}

fn record(
    qname: &[u8],
    pos0: i64,
    mate_pos0: i64,
    read1: bool,
    reverse: bool,
    mate_reverse: bool,
    proper: bool,
    cigar: Vec<Cigar>,
    sequence: &[u8],
    qualities: &[u8],
) -> Record {
    let mut record = Record::new();
    record.set(qname, Some(&CigarString(cigar)), sequence, qualities);
    record.set_tid(0);
    record.set_mtid(0);
    record.set_pos(pos0);
    record.set_mpos(mate_pos0);
    record.set_mapq(30);
    let mut flags = 0x1 | if read1 { 0x40 } else { 0x80 };
    if proper {
        flags |= 0x2;
    }
    if reverse {
        flags |= 0x10;
    }
    if mate_reverse {
        flags |= 0x20;
    }
    record.set_flags(flags);
    record
}

fn inward_pair(
    forward_pos0: i64,
    reverse_pos0: i64,
    read_len: usize,
    forward_sequence: Vec<u8>,
    reverse_sequence: Vec<u8>,
    forward_qualities: Vec<u8>,
    reverse_qualities: Vec<u8>,
) -> (Record, Record) {
    (
        record(
            b"pair",
            forward_pos0,
            reverse_pos0,
            true,
            false,
            true,
            true,
            vec![Cigar::Match(read_len as u32)],
            &forward_sequence,
            &forward_qualities,
        ),
        record(
            b"pair",
            reverse_pos0,
            forward_pos0,
            false,
            true,
            false,
            true,
            vec![Cigar::Match(read_len as u32)],
            &reverse_sequence,
            &reverse_qualities,
        ),
    )
}

fn ordinary_pair(forward_pos0: i64, reverse_pos0: i64, read_len: usize) -> (Record, Record) {
    inward_pair(
        forward_pos0,
        reverse_pos0,
        read_len,
        vec![b'A'; read_len],
        vec![b'A'; read_len],
        vec![30; read_len],
        vec![30; read_len],
    )
}

fn evaluate(
    first: &Record,
    second: &Record,
) -> std::result::Result<bam2frag::pair_eligibility::PairCandidate, PairRejection> {
    let cfg = natural_resolved_config();
    let (_temp, mut reference) = reference_cache();
    construct_pair_candidate(first, second, "chr1", CHROM_LEN, &cfg, &mut reference).unwrap()
}

#[test]
fn endpoints_are_zero_based_and_span_is_inclusive() {
    let (read1, read2) = ordinary_pair(500, 530, 30);
    let candidate = evaluate(&read1, &read2).unwrap();
    assert_eq!(candidate.left_5p_genomic, 500);
    assert_eq!(candidate.right_5p_genomic, 559);
    assert_eq!(candidate.fragment_span_bp, 60);
    assert_eq!(
        candidate.window_end0 - candidate.window_start0,
        MOLECULAR_WINDOW_BP as i64
    );
    assert_eq!(candidate.window_start0, 370);
}

#[test]
fn accepts_abutment_overlap_and_full_overlap() {
    for (reverse_start0, expected_span) in [(530, 60), (529, 59), (520, 50), (500, 30)] {
        let (read1, read2) = ordinary_pair(500, reverse_start0, 30);
        assert_eq!(
            evaluate(&read1, &read2).unwrap().fragment_span_bp,
            expected_span
        );
    }
}

#[test]
fn rejects_one_base_internal_gap_and_protruding_geometry() {
    let (read1, read2) = ordinary_pair(500, 531, 30);
    assert_eq!(
        evaluate(&read1, &read2).unwrap_err(),
        PairRejection::IncompleteFragmentCoverage
    );

    let (read1, read2) = ordinary_pair(500, 490, 40);
    assert_eq!(
        evaluate(&read1, &read2).unwrap_err(),
        PairRejection::Protruding3p
    );
}

#[test]
fn rejects_clipping_insertions_and_deletions() {
    for (cigar, sequence_len, expected) in [
        (
            vec![Cigar::SoftClip(1), Cigar::Match(29)],
            30,
            PairRejection::Clipping,
        ),
        (
            vec![Cigar::HardClip(1), Cigar::Match(30)],
            30,
            PairRejection::Clipping,
        ),
        (
            vec![Cigar::Match(10), Cigar::Ins(1), Cigar::Match(19)],
            30,
            PairRejection::IndelCigar,
        ),
        (
            vec![Cigar::Match(10), Cigar::Del(1), Cigar::Match(20)],
            30,
            PairRejection::IndelCigar,
        ),
    ] {
        let read1 = record(
            b"pair",
            500,
            500,
            true,
            false,
            true,
            true,
            cigar,
            &vec![b'A'; sequence_len],
            &vec![30; sequence_len],
        );
        let (_, read2) = ordinary_pair(500, 500, 30);
        assert_eq!(evaluate(&read1, &read2).unwrap_err(), expected);
        assert_eq!(evaluate(&read2, &read1).unwrap_err(), expected);
    }
}

#[test]
fn enforces_orientation_pair_flags_and_mapq_independently() {
    let (mut read1, mut read2) = ordinary_pair(500, 500, 30);
    read1.set_flags(read1.flags() & !0x20);
    read2.set_flags(read2.flags() & !0x10);
    assert_eq!(
        evaluate(&read1, &read2).unwrap_err(),
        PairRejection::NonInwardOrientation
    );

    for flag in [0x100, 0x400, 0x800, 0x200] {
        let (mut read1, read2) = ordinary_pair(500, 500, 30);
        read1.set_flags(read1.flags() | flag);
        assert_eq!(
            evaluate(&read1, &read2).unwrap_err(),
            PairRejection::InvalidFlags
        );
    }

    let (mut read1, read2) = ordinary_pair(500, 500, 30);
    read1.set_flags(read1.flags() & !0x2);
    assert_eq!(
        evaluate(&read1, &read2).unwrap_err(),
        PairRejection::ImproperPair
    );

    let (mut read1, read2) = ordinary_pair(500, 500, 30);
    read1.set_flags(read1.flags() | 0x8);
    assert_eq!(
        evaluate(&read1, &read2).unwrap_err(),
        PairRejection::UnpairedOrMateUnmapped
    );

    let (mut read1, read2) = ordinary_pair(500, 500, 30);
    read1.set_mapq(29);
    assert_eq!(evaluate(&read1, &read2).unwrap_err(), PairRejection::Mapq);
    read1.set_mapq(30);
    assert!(evaluate(&read1, &read2).is_ok());

    let (read1, mut read2) = ordinary_pair(500, 500, 30);
    read2.set_tid(1);
    assert_eq!(
        evaluate(&read1, &read2).unwrap_err(),
        PairRejection::DifferentChromosomes
    );
}

#[test]
fn enforces_inclusive_minimum_and_maximum_span() {
    for (forward_start0, reverse_start0, read_len, expected_span) in
        [(500, 500, 30, 30), (500, 649, 150, 299)]
    {
        let (read1, read2) = ordinary_pair(forward_start0, reverse_start0, read_len);
        assert_eq!(
            evaluate(&read1, &read2).unwrap().fragment_span_bp,
            expected_span
        );
    }
    let (read1, read2) = ordinary_pair(500, 500, 29);
    assert_eq!(
        evaluate(&read1, &read2).unwrap_err(),
        PairRejection::FragmentSpan
    );
    let (read1, read2) = ordinary_pair(500, 650, 150);
    assert_eq!(
        evaluate(&read1, &read2).unwrap_err(),
        PairRejection::FragmentSpan
    );
}

#[test]
fn records_every_kind_of_mate_disagreement_without_rejecting_the_pair() {
    let mut forward = vec![b'A'; 30];
    let mut reverse = vec![b'A'; 30];

    // Concordant: both mates agree on the substitution.
    forward[5] = b'C';
    reverse[5] = b'C';
    let (read1, read2) = inward_pair(
        500,
        500,
        30,
        forward.clone(),
        reverse.clone(),
        vec![30; 30],
        vec![30; 30],
    );
    let pair = evaluate(&read1, &read2).unwrap();
    assert_eq!(pair.observation.concordant_mismatches, 1);
    assert_eq!(pair.observation.unilateral_disagreements, 0);
    assert_eq!(pair.observation.discordant_disagreements, 0);

    // Unilateral: the reverse mate corroborates the reference instead.
    reverse[5] = b'A';
    let (read1, read2) = inward_pair(
        500,
        500,
        30,
        forward.clone(),
        reverse.clone(),
        vec![30; 30],
        vec![30; 30],
    );
    let pair = evaluate(&read1, &read2).unwrap();
    assert_eq!(pair.observation.unilateral_disagreements, 1);
    assert_eq!(pair.observation.concordant_mismatches, 0);

    // Discordant: both mates disagree with the reference and with each other.
    reverse[5] = b'G';
    let (read1, read2) = inward_pair(
        500,
        500,
        30,
        forward.clone(),
        reverse.clone(),
        vec![30; 30],
        vec![30; 30],
    );
    let pair = evaluate(&read1, &read2).unwrap();
    assert_eq!(pair.observation.discordant_disagreements, 1);
    assert_eq!(pair.observation.concordant_mismatches, 0);

    // Concordant at a quality that used to reject the pair. Extraction applies
    // no base-quality threshold at all, so this is admitted and binned.
    //
    // Both mates disagree here, so the forward mate's quality is the one binned,
    // and Q30 lands in the middle bin rather than the high one: the edges follow
    // four-bin Illumina quality binning, where the high bin is Q37. Q30 is a
    // fixture value, not a value these instruments emit.
    reverse[5] = b'C';
    let mut low_quality = vec![30; 30];
    low_quality[5] = 11;
    let (read1, read2) = inward_pair(500, 500, 30, forward, reverse, vec![30; 30], low_quality);
    let pair = evaluate(&read1, &read2).unwrap();
    assert_eq!(pair.observation.concordant_mismatches, 1);
    assert_eq!(pair.observation.disagreements_by_quality_bin, [0, 0, 1, 0]);
}

#[test]
fn singly_covered_bases_keep_both_low_quality_and_mismatching_observations() {
    let mut forward = vec![b'A'; 30];
    let reverse = vec![b'A'; 30];

    // A singly covered reference-matching base at quality zero is retained.
    let mut forward_quality = vec![30; 30];
    forward_quality[5] = 0;
    let (read1, read2) = inward_pair(
        500,
        530,
        30,
        forward.clone(),
        reverse.clone(),
        forward_quality,
        vec![30; 30],
    );
    assert!(evaluate(&read1, &read2).is_ok());

    // So is a singly covered mismatch. Only one mate reaches this position, so
    // there is nothing to corroborate it with; rejecting the fragment would
    // delete every molecule spanning a germline variant outside the overlap,
    // and would do so more often the longer the fragment.
    forward[5] = b'C';
    let (read1, read2) = inward_pair(500, 530, 30, forward, reverse, vec![30; 30], vec![30; 30]);
    let pair = evaluate(&read1, &read2).unwrap();
    assert_eq!(pair.observation.unilateral_disagreements, 1);
    assert_eq!(pair.observation.span_bp, 60);
    assert_eq!(pair.observation.covered_positions, 60);
}

#[test]
fn mate_arrival_order_does_not_change_candidate_or_rejection() {
    let (read1, read2) = ordinary_pair(500, 520, 30);
    let forward = evaluate(&read1, &read2).unwrap();
    let reverse = evaluate(&read2, &read1).unwrap();
    assert_eq!(forward.left_5p_genomic, reverse.left_5p_genomic);
    assert_eq!(forward.right_5p_genomic, reverse.right_5p_genomic);
    assert_eq!(forward.fragment_hash, reverse.fragment_hash);
}

#[test]
fn rejects_incomplete_centered_reference_window() {
    let (read1, read2) = ordinary_pair(10, 10, 30);
    assert_eq!(
        evaluate(&read1, &read2).unwrap_err(),
        PairRejection::IncompleteReferenceWindow
    );
}

#[test]
fn rejection_counters_are_named_and_independent() {
    let mut counters = PairRejectionCounters::default();
    let (read1, read2) = ordinary_pair(500, 500, 30);
    counters.record(&evaluate(&read1, &read2));
    let (read1, read2) = ordinary_pair(500, 531, 30);
    counters.record(&evaluate(&read1, &read2));
    let (mut read1, read2) = ordinary_pair(500, 500, 30);
    read1.set_mapq(29);
    counters.record(&evaluate(&read1, &read2));

    assert_eq!(counters.pairs_accepted, 1);
    assert_eq!(counters.pairs_rejected_incomplete_fragment_coverage, 1);
    assert_eq!(counters.pairs_rejected_mapq, 1);
}
