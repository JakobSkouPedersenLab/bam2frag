use bam2frag::{
    model::{ReservoirKey, ValidatedFragmentCandidate},
    pair_eligibility::{FragmentObservation, PairCandidate, PairOrientation},
    reservoir::KeyedReservoir,
};
use rust_htslib::bam::Record;

fn candidate(id: i64) -> ValidatedFragmentCandidate {
    let left = 1_000_000 + id * 2;
    let right = left + 159;
    ValidatedFragmentCandidate {
        key: ReservoirKey {
            region_id: 7,
            local_pos: 0,
        },
        pair: PairCandidate {
            read1: Record::new(),
            read2: Record::new(),
            chrom: "chr1".into(),
            left_5p_genomic: left,
            right_5p_genomic: right,
            fragment_span_bp: 160,
            window_start0: left - 80,
            window_end0: left + 240,
            orientation: PairOrientation::R1ForwardR2Reverse,
            fragment_hash: id as u64 ^ 0x9E37_79B9_7F4A_7C15,
            // Reservoir admission does not read the observation tallies, so a
            // fully reference-matching fragment keeps this fixture minimal.
            observation: FragmentObservation {
                span_bp: 160,
                covered_positions: 160,
                ..FragmentObservation::default()
            },
        },
        anchor_rpos0: left + 79,
    }
}

#[test]
fn natural_reservoir_retains_exact_capacity_without_event_strata() {
    let mut reservoir = KeyedReservoir::new(10_000, 17, None);
    for id in 0..12_000i64 {
        reservoir.consider(candidate(id)).unwrap();
    }
    let (candidates, stats) = reservoir.into_parts();
    assert_eq!(candidates.len(), 10_000);
    assert_eq!(stats.len(), 1);
    assert_eq!(stats[0].available_count, 12_000);
    assert_eq!(stats[0].sampled_count, 10_000);
}

#[test]
fn natural_reservoir_is_deterministic_for_fixed_seed_and_fragment_hash() {
    let draw = || {
        let mut reservoir = KeyedReservoir::new(128, 1234, None);
        for id in 0..500i64 {
            reservoir.consider(candidate(id)).unwrap();
        }
        reservoir
            .into_parts()
            .0
            .into_iter()
            .map(|candidate| candidate.pair.fragment_hash)
            .collect::<Vec<_>>()
    };
    assert_eq!(draw(), draw());
}
