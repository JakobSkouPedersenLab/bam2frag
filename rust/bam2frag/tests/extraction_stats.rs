//! Persisted extraction QC counters for the natural `extract` command.
//!
//! `extraction_stats.csv` and `extraction_stats.json` must carry every named
//! pair-eligibility rejection reason plus the accepted-but-motif-ineligible
//! count, under one stable set of field names, with a deterministic total.

mod common;

use bam2frag::{
    features::{end_motif_min_quality, record_pair_decision},
    model::{ChromStats, ReservoirKey, ValidatedFragmentCandidate},
    output::{write_stats, EXTRACTION_STATS_FORMAT},
    pair_eligibility::{
        construct_pair_candidate, quality_bin, span_bin, FragmentObservation,
        FragmentObservationStats, PairCandidate, PairRejection, PairRejectionCounters,
        PAIR_COUNTER_FIELDS,
    },
    reference::ReferenceCache,
    reservoir::KeyedReservoir,
};
use common::natural_resolved_config;
use rust_htslib::bam::{
    record::{Cigar, CigarString},
    Record,
};
use serde_json::Value;
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
    pos0: i64,
    mate_pos0: i64,
    read1: bool,
    reverse: bool,
    mate_reverse: bool,
    sequence: &[u8],
    qualities: &[u8],
) -> Record {
    let mut record = Record::new();
    record.set(
        b"pair",
        Some(&CigarString(vec![Cigar::Match(sequence.len() as u32)])),
        sequence,
        qualities,
    );
    record.set_tid(0);
    record.set_mtid(0);
    record.set_pos(pos0);
    record.set_mpos(mate_pos0);
    record.set_mapq(30);
    let mut flags = 0x1 | 0x2 | if read1 { 0x40 } else { 0x80 };
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
    forward_sequence: Vec<u8>,
    reverse_sequence: Vec<u8>,
    forward_qualities: Vec<u8>,
    reverse_qualities: Vec<u8>,
) -> (Record, Record) {
    (
        record(
            forward_pos0,
            reverse_pos0,
            true,
            false,
            true,
            &forward_sequence,
            &forward_qualities,
        ),
        record(
            reverse_pos0,
            forward_pos0,
            false,
            true,
            false,
            &reverse_sequence,
            &reverse_qualities,
        ),
    )
}

fn decide(first: &Record, second: &Record) -> std::result::Result<PairCandidate, PairRejection> {
    let cfg = natural_resolved_config();
    let (_temp, mut reference) = reference_cache();
    construct_pair_candidate(first, second, "chr1", CHROM_LEN, &cfg, &mut reference).unwrap()
}

fn chrom_stats(chrom: &str, reads: u64, counters: PairRejectionCounters) -> ChromStats {
    ChromStats {
        chrom: chrom.to_string(),
        reads,
        skipped: 1,
        evicted: 2,
        invalid_fragments: 3,
        pair_rejections: counters,
        observations: FragmentObservationStats::default(),
    }
}

fn distinct_counters(step: u64) -> PairRejectionCounters {
    // Give every field a different value so a mis-wired merge or column order
    // cannot pass by coincidence.
    let mut counters = PairRejectionCounters::default();
    counters.pairs_accepted = step;
    counters.pairs_rejected_non_inward_orientation = step + 1;
    counters.pairs_rejected_ambiguous_5p = step + 2;
    counters.pairs_rejected_5p_clipping = step + 3;
    counters.pairs_rejected_protruding_3p = step + 4;
    counters.pairs_rejected_unpaired_or_mate_unmapped = step + 5;
    counters.pairs_rejected_improper_pair = step + 6;
    counters.pairs_rejected_invalid_flags = step + 7;
    counters.pairs_rejected_different_chromosomes = step + 8;
    counters.pairs_rejected_mapq = step + 9;
    counters.pairs_rejected_fragment_span = step + 10;
    counters.pairs_rejected_incomplete_fragment_coverage = step + 11;
    counters.pairs_rejected_reference_mismatch_burden = step + 12;
    counters.pairs_rejected_indel_cigar = step + 13;
    counters.pairs_rejected_invalid_5p_span = step + 14;
    counters.pairs_rejected_incomplete_reference_window = step + 15;
    counters.accepted_covered_positions = step + 16;
    counters.accepted_unilateral_disagreements = step + 17;
    counters.accepted_discordant_disagreements = step + 18;
    counters.accepted_concordant_mismatches = step + 19;
    counters.accepted_without_end_motif_bases = step + 20;
    counters
}

fn write_both(rows: &[ChromStats]) -> (TempDir, Vec<csv::StringRecord>, Value) {
    let temp = tempfile::tempdir().unwrap();
    let csv_path = temp.path().join("extraction_stats.csv");
    let json_path = temp.path().join("extraction_stats.json");
    write_stats(&csv_path, &json_path, rows, "S1", "P1").unwrap();

    let mut reader = csv::Reader::from_path(&csv_path).unwrap();
    let header = reader.headers().unwrap().clone();
    let mut records = vec![header];
    for row in reader.records() {
        records.push(row.unwrap());
    }
    let json: Value = serde_json::from_slice(&fs::read(&json_path).unwrap()).unwrap();
    (temp, records, json)
}

#[test]
fn every_rejection_reason_is_persisted_under_a_stable_field_name() {
    let rows = vec![chrom_stats("chr1", 10, distinct_counters(100))];
    let (_temp, records, json) = write_both(&rows);

    // The CSV header carries the counters in the declared order, between the
    // base counts and the identity columns.
    let header: Vec<&str> = records[0].iter().collect();
    assert_eq!(
        header[..5],
        [
            "chrom",
            "representative_reads",
            "skipped",
            "evicted",
            "invalid_fragments"
        ]
    );
    assert_eq!(
        header[5..5 + PAIR_COUNTER_FIELDS.len()],
        PAIR_COUNTER_FIELDS
    );
    assert_eq!(
        header[5 + PAIR_COUNTER_FIELDS.len()..],
        ["patient_id", "sample_id"]
    );

    // The declared order must match the derived serde field names exactly, so
    // renaming a struct field cannot silently desynchronize CSV and JSON.
    let serialized = serde_json::to_value(distinct_counters(0)).unwrap();
    let derived: Vec<String> = serialized.as_object().unwrap().keys().cloned().collect();
    assert_eq!(derived, PAIR_COUNTER_FIELDS.to_vec());

    assert_eq!(json["format"], EXTRACTION_STATS_FORMAT);
    assert_eq!(json["sample_id"], "S1");
    assert_eq!(json["patient_id"], "P1");
    for field in PAIR_COUNTER_FIELDS {
        assert!(
            json["total"].get(field).is_some(),
            "{field} missing from extraction_stats.json total"
        );
        assert!(
            json["chromosomes"][0].get(field).is_some(),
            "{field} missing from extraction_stats.json chromosome row"
        );
    }
}

#[test]
fn chromosome_counters_merge_by_summing_every_field() {
    let mut merged = distinct_counters(100);
    merged.merge(&distinct_counters(200));
    for (index, (name, value)) in merged.as_named_counts().into_iter().enumerate() {
        let expected = (100 + index as u64) + (200 + index as u64);
        assert_eq!(value, expected, "{name} merged incorrectly");
    }

    let mut total = ChromStats::new("total");
    total.merge(&chrom_stats("chr1", 10, distinct_counters(100)));
    total.merge(&chrom_stats("chr2", 20, distinct_counters(200)));
    assert_eq!(total.reads, 30);
    assert_eq!(total.skipped, 2);
    assert_eq!(total.evicted, 4);
    assert_eq!(total.invalid_fragments, 6);
    assert_eq!(total.pair_rejections, merged);
}

#[test]
fn csv_and_json_totals_agree_and_sum_the_chromosome_rows() {
    let rows = vec![
        chrom_stats("chr1", 10, distinct_counters(100)),
        chrom_stats("chr2", 20, distinct_counters(200)),
    ];
    let (_temp, records, json) = write_both(&rows);

    // header + two chromosomes + total
    assert_eq!(records.len(), 4);
    let header: Vec<&str> = records[0].iter().collect();
    let total_row = &records[3];
    assert_eq!(&total_row[0], "total");

    let mut expected = ChromStats::new("total");
    for row in &rows {
        expected.merge(row);
    }

    for (name, value) in expected.pair_rejections.as_named_counts() {
        let column = header.iter().position(|field| *field == name).unwrap();
        assert_eq!(
            total_row[column].parse::<u64>().unwrap(),
            value,
            "{name} disagrees in the CSV total row"
        );
        assert_eq!(
            json["total"][name].as_u64().unwrap(),
            value,
            "{name} disagrees in the JSON total"
        );
    }
    assert_eq!(json["chromosomes"].as_array().unwrap().len(), 2);
    assert_eq!(json["chromosomes"][0]["chrom"], "chr1");
    assert_eq!(json["total"]["representative_reads"].as_u64().unwrap(), 30);
    assert_eq!(
        total_row[1].parse::<u64>().unwrap(),
        json["total"]["representative_reads"].as_u64().unwrap()
    );
}

#[test]
fn a_low_quality_five_prime_motif_base_is_counted_but_never_rejected() {
    // Abutting pair, every base matching the reference. Only the forward read's
    // first sequenced cycle is below Q30, and that position is covered by one
    // mate only, so no mismatch rule applies to it.
    let mut forward_quality = vec![30u8; 30];
    forward_quality[0] = 11;
    let (read1, read2) = inward_pair(
        500,
        530,
        vec![b'A'; 30],
        vec![b'A'; 30],
        forward_quality,
        vec![30; 30],
    );

    let decision = decide(&read1, &read2);
    let candidate = decision
        .as_ref()
        .expect("a low-quality reference-matching motif base must not reject the pair")
        .clone();
    // The motif bases are all A/C/G/T, so the pair contributes its lowest motif
    // quality to the histogram. Q11 lands in the low bin; extraction does not
    // judge it, and a reader picks its own threshold from this
    // distribution.
    assert_eq!(end_motif_min_quality(&read1, &read2), Some(11));

    let mut counters = PairRejectionCounters::default();
    let mut observations = FragmentObservationStats::default();
    record_pair_decision(&mut counters, &mut observations, &decision, &read1, &read2);
    assert_eq!(counters.pairs_accepted, 1);
    assert_eq!(counters.accepted_without_end_motif_bases, 0);
    assert_eq!(observations.end_motif_min_quality_histogram, [0, 1, 0, 0]);
    assert_eq!(counters.pairs_rejected_reference_mismatch_burden, 0);

    // The fragment still reaches the natural reservoir and is retained there.
    let mut reservoir = KeyedReservoir::new(10_000, 17, None);
    reservoir
        .consider(ValidatedFragmentCandidate {
            key: ReservoirKey {
                region_id: 7,
                local_pos: 0,
            },
            anchor_rpos0: candidate.left_5p_genomic,
            pair: candidate,
        })
        .unwrap();
    let (retained, stats) = reservoir.into_parts();
    assert_eq!(retained.len(), 1);
    assert_eq!(stats[0].sampled_count, 1);
}

#[test]
fn a_low_quality_mismatch_is_accepted_and_counted_rather_than_rejected() {
    // Both mates cover position 505 and agree on the substitution, but one
    // observation is below Q30. Extraction no longer thresholds base quality:
    // the pair is accepted, the low quality reaches the model on the fragment bank
    // quality plane, and the disagreement is recorded. The eight motif bases
    // are all Q30 A/C/G/T, so the motif counter must not move.
    let mut substituted = vec![b'A'; 30];
    substituted[5] = b'C';
    let mut low_quality = vec![30u8; 30];
    low_quality[5] = 11;
    let (read1, read2) = inward_pair(
        500,
        500,
        substituted.clone(),
        substituted,
        vec![30; 30],
        low_quality,
    );

    let decision = decide(&read1, &read2);
    assert!(decision.is_ok());
    assert_eq!(end_motif_min_quality(&read1, &read2), Some(30));

    let mut counters = PairRejectionCounters::default();
    let mut observations = FragmentObservationStats::default();
    record_pair_decision(&mut counters, &mut observations, &decision, &read1, &read2);
    assert_eq!(counters.pairs_accepted, 1);
    assert_eq!(counters.accepted_concordant_mismatches, 1);
    assert_eq!(counters.accepted_unilateral_disagreements, 0);
    assert_eq!(counters.accepted_covered_positions, 30);
    assert_eq!(counters.accepted_without_end_motif_bases, 0);

    // Span 30 lands in the first bin, and one mismatch in thirty positions is
    // 3.3%, which lands in histogram bucket 3.
    assert_eq!(observations.accepted_by_span_bin[0], 1);
    assert_eq!(observations.disagreements_by_span_bin[0], 1);
    assert_eq!(observations.mismatch_fraction_histogram[3], 1);
}

#[test]
fn an_ambiguous_end_motif_base_is_counted_separately_from_a_low_quality_one() {
    // An `N` in a motif position is not a quality judgement: there is no base to
    // reconstruct, so the fragment is counted apart from the quality histogram
    // rather than being scored as low quality.
    let mut forward = vec![b'A'; 30];
    forward[0] = b'N';
    let (read1, read2) = inward_pair(
        500,
        530,
        forward,
        vec![b'A'; 30],
        vec![30; 30],
        vec![30; 30],
    );

    let decision = decide(&read1, &read2);
    assert!(
        decision.is_ok(),
        "an ambiguous motif base must not reject the pair"
    );
    assert_eq!(end_motif_min_quality(&read1, &read2), None);

    let mut counters = PairRejectionCounters::default();
    let mut observations = FragmentObservationStats::default();
    record_pair_decision(&mut counters, &mut observations, &decision, &read1, &read2);
    assert_eq!(counters.pairs_accepted, 1);
    assert_eq!(counters.accepted_without_end_motif_bases, 1);
    assert_eq!(observations.end_motif_min_quality_histogram, [0; 4]);
}

#[test]
fn observational_histograms_merge_and_sum_deterministically() {
    let mut left = FragmentObservationStats::default();
    let mut right = FragmentObservationStats::default();

    // A short fully-reference-matching fragment, and a longer one carrying one
    // unilateral disagreement at a low-quality base.
    left.record_accepted(&FragmentObservation {
        span_bp: 60,
        covered_positions: 60,
        ..FragmentObservation::default()
    });
    right.record_accepted(&FragmentObservation {
        span_bp: 200,
        covered_positions: 200,
        unilateral_disagreements: 2,
        disagreements_by_quality_bin: [0, 2, 0, 0],
        ..FragmentObservation::default()
    });

    left.record_end_motif_quality(37);
    right.record_end_motif_quality(11);

    left.merge(&right);
    assert_eq!(left.end_motif_min_quality_histogram, [0, 1, 0, 1]);
    assert_eq!(left.accepted_by_span_bin, [1, 0, 0, 1, 0]);
    assert_eq!(left.covered_positions_by_span_bin, [60, 0, 0, 200, 0]);
    assert_eq!(left.disagreements_by_span_bin, [0, 0, 0, 2, 0]);
    assert_eq!(left.disagreements_by_quality_bin, [0, 2, 0, 0]);
    // 0/60 lands in bucket 0; 2/200 = 1% lands in bucket 1.
    assert_eq!(left.mismatch_fraction_histogram[0], 1);
    assert_eq!(left.mismatch_fraction_histogram[1], 1);
}

#[test]
fn span_and_quality_bins_cover_their_whole_range() {
    // Bin edges are part of the persisted diagnostic, so pin them.
    assert_eq!(span_bin(30), 0);
    assert_eq!(span_bin(99), 0);
    assert_eq!(span_bin(100), 1);
    assert_eq!(span_bin(149), 1);
    assert_eq!(span_bin(150), 2);
    assert_eq!(span_bin(199), 2);
    assert_eq!(span_bin(200), 3);
    assert_eq!(span_bin(249), 3);
    assert_eq!(span_bin(250), 4);
    assert_eq!(span_bin(299), 4);

    // Four-bin Illumina quality binning: no-call, low, middle, high. The edges
    // straddle both observed bin layouts, 2/12/23/37 and 2/11/25/37.
    assert_eq!(quality_bin(2), 0);
    assert_eq!(quality_bin(11), 1);
    assert_eq!(quality_bin(12), 1);
    assert_eq!(quality_bin(23), 2);
    assert_eq!(quality_bin(25), 2);
    assert_eq!(quality_bin(37), 3);

    // Q30 is the middle bin, not the high one: the high bin is the Q37 bin, and
    // a binned instrument never emits Q30 at all. It reads as "high quality" in
    // prose and in test fixtures, which makes it the easy value to assert wrong.
    assert_eq!(quality_bin(30), 2);
    assert_eq!(quality_bin(35), 2);
    assert_eq!(quality_bin(36), 3);
}

#[test]
fn the_mismatch_burden_cap_saturates_its_histogram_bucket() {
    // Everything at or above ten percent collects in the final bucket, so the
    // histogram cannot silently run off its end.
    let observation = FragmentObservation {
        span_bp: 100,
        covered_positions: 100,
        concordant_mismatches: 50,
        ..FragmentObservation::default()
    };
    assert_eq!(observation.mismatch_fraction(), 0.5);
    assert_eq!(observation.mismatch_fraction_bin(), 10);

    let observation = FragmentObservation {
        span_bp: 100,
        covered_positions: 100,
        concordant_mismatches: 9,
        ..FragmentObservation::default()
    };
    assert_eq!(observation.mismatch_fraction_bin(), 9);
}
