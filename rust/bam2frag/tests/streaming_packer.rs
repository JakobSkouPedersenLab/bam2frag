mod common;

use bam2frag::{
    arrays::{decode_payload, decompress_payload, pack_region, payload_nbytes},
    model::{MaterializedFragment, ReservoirKey},
};
use common::natural_resolved_config;

fn synthetic_fragment(anchor: i64, width: usize) -> MaterializedFragment {
    let orientation_code = (anchor % 2) as u8;
    let mut r1_aln = vec![0; width];
    let mut r2_aln = vec![0; width];
    let mut r1_qual = vec![0; width];
    let mut r2_qual = vec![0; width];
    let mut r1_cycle = vec![0; width];
    let mut r2_cycle = vec![0; width];
    let (r1_positions, r2_positions): (Vec<usize>, Vec<usize>) = if orientation_code == 0 {
        ((85..235).collect(), (85..235).rev().collect())
    } else {
        ((85..235).rev().collect(), (85..235).collect())
    };
    for (cycle0, position) in r1_positions.into_iter().enumerate() {
        r1_aln[position] = 2;
        r1_qual[position] = 37;
        r1_cycle[position] = (cycle0 + 1) as u8;
    }
    for (cycle0, position) in r2_positions.into_iter().enumerate() {
        r2_aln[position] = 3;
        r2_qual[position] = 41;
        r2_cycle[position] = (cycle0 + 1) as u8;
    }
    MaterializedFragment {
        key: ReservoirKey { region_id: 7 },
        anchor_rpos0: anchor,
        ref_aln: vec![1; width],
        r1_aln,
        r2_aln,
        r1_qual,
        r2_qual,
        r1_cycle,
        r2_cycle,
        fragment_span_bp: 150,
        orientation_code,
        atac_overlap_fractions: [0.25, 0.5, 0.75],
    }
}

#[test]
fn dnapack8_payload_has_exact_expected_size_and_no_event_plane() {
    let cfg = natural_resolved_config();
    let width = cfg.typed.feature_sets.fragment_sequence.seq_window_w;
    let fragments = (0..10_000)
        .map(|index| synthetic_fragment(index, width))
        .collect::<Vec<_>>();
    let packed = pack_region(7, &fragments, &cfg).unwrap();
    let total = fragments.len();
    let dna = 3 * (total * width * 4).div_ceil(8);
    let quality = 2 * (total * width * 6).div_ceil(8);
    let cycles = 2 * total * width;
    let numeric = total * 3 * 4;
    let structural = total * 4;
    let expected = 112 + dna + quality + cycles + numeric + structural;
    assert_eq!(packed.payloads[0].payload_nbytes as usize, expected);
    assert_eq!(payload_nbytes(total).unwrap(), expected);
    let payload = decompress_payload(&packed.payloads[0].frame, expected).unwrap();
    assert_eq!(&payload[..8], b"fragment bank");
    assert_eq!(payload[11], 9);
    assert!(packed.payloads[0].frame.len() < expected);
}

#[test]
fn dnapack8_reference_packer_round_trips_every_plane() {
    let cfg = natural_resolved_config();
    let width = cfg.typed.feature_sets.fragment_sequence.seq_window_w;
    let fragment = synthetic_fragment(11, width);
    let packed = pack_region(7, &[fragment.clone()], &cfg).unwrap();
    let payload =
        decompress_payload(&packed.payloads[0].frame, payload_nbytes(1).unwrap()).unwrap();
    let decoded = decode_payload(&payload).unwrap();
    assert_eq!(decoded.ref_aln, fragment.ref_aln);
    assert_eq!(decoded.r1_aln, fragment.r1_aln);
    assert_eq!(decoded.r2_aln, fragment.r2_aln);
    assert_eq!(decoded.r1_qual, fragment.r1_qual);
    assert_eq!(decoded.r2_qual, fragment.r2_qual);
    assert_eq!(decoded.r1_cycle, fragment.r1_cycle);
    assert_eq!(decoded.r2_cycle, fragment.r2_cycle);
    assert_eq!(decoded.fragment_span_bp, vec![fragment.fragment_span_bp]);
    assert_eq!(decoded.orientation_code, vec![fragment.orientation_code]);
    assert_eq!(
        decoded.atac_overlap_fractions,
        vec![fragment.atac_overlap_fractions]
    );
}
