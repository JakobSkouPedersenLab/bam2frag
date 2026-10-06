mod common;

use bam2frag::arrays::{
    pack_fixed_bits, unpack_fixed_bits, CompressedPayload, FragmentBankWriter, PackedRegion,
    PACKED_BANK_FILENAME, PACKED_INDEX_FILENAME,
};
use common::{load_from_value, valid_config_json};
use std::fs;
use tempfile::tempdir;

fn le_u16(bytes: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes(bytes[offset..offset + 2].try_into().unwrap())
}
fn le_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().unwrap())
}
fn le_i64(bytes: &[u8], offset: usize) -> i64 {
    i64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}
fn le_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().unwrap())
}

/// The writer stores whatever frame it is handed; these tests check the layout.
fn frame(bytes: &[u8], payload_nbytes: u64) -> CompressedPayload {
    CompressedPayload {
        frame: bytes.to_vec(),
        payload_nbytes,
    }
}

#[test]
fn fixed_bit_packing_round_trips_four_six_and_one_bit_values() {
    for (values, bits) in [
        ((0u8..16).collect::<Vec<_>>(), 4),
        ((0u8..64).collect::<Vec<_>>(), 6),
        (vec![0, 1, 1, 0, 1, 0, 0, 1, 1], 1),
    ] {
        let packed = pack_fixed_bits(&values, bits).unwrap();
        assert_eq!(unpack_fixed_bits(&packed, values.len(), bits), values);
    }
}

#[test]
fn natural_fragment_bank_writer_emits_one_payload_and_dnaidx8_at_10k() {
    let tmp = tempdir().unwrap();
    let regions = vec![
        PackedRegion {
            region_id: 7,
            local_pos: 0,
            payloads: vec![frame(b"natural-a", 1_000)],
            counts: vec![10_000],
        },
        PackedRegion {
            region_id: 8,
            local_pos: 0,
            payloads: vec![frame(b"natural-b", 500)],
            counts: vec![777],
        },
    ];

    let mut writer = FragmentBankWriter::create(tmp.path(), regions.len()).unwrap();
    for region in regions {
        writer.write_region(region).unwrap();
    }
    assert_eq!(writer.byte_totals(), (1_500, 18));
    let rows = writer.finish(tmp.path()).unwrap();

    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].payloads.len(), 1);
    assert_eq!(rows[0].payloads[0].count, 10_000);
    assert_eq!(rows[1].payloads[0].count, 777);

    let bank = fs::read(tmp.path().join(PACKED_BANK_FILENAME)).unwrap();
    assert_eq!(bank, b"natural-anatural-b");

    let index = fs::read(tmp.path().join(PACKED_INDEX_FILENAME)).unwrap();
    assert_eq!(&index[..8], b"bank index\0");
    assert_eq!(le_u16(&index, 8), 8);
    assert_eq!(le_u16(&index, 10), 1);
    assert_eq!(le_u32(&index, 12), 2);

    let row_size = 4 + 8 + 8 + 20;
    let row0 = 16;
    let row1 = row0 + row_size;
    assert_eq!(le_u32(&index, row0), 0);
    assert_eq!(le_i64(&index, row0 + 4), 7);
    assert_eq!(le_u32(&index, row1), 1);
    assert_eq!(le_i64(&index, row1 + 4), 8);
    // Offsets and lengths are those of the stored frames.
    assert_eq!(le_u64(&index, row0 + 24), 0);
    assert_eq!(le_u64(&index, row0 + 32), 9);
    assert_eq!(le_u64(&index, row1 + 24), 9);
    assert_eq!(le_u64(&index, row1 + 32), 9);
}

#[test]
fn natural_writer_rejects_an_empty_frame() {
    let tmp = tempdir().unwrap();
    let mut writer = FragmentBankWriter::create(tmp.path(), 1).unwrap();
    let err = writer
        .write_region(PackedRegion {
            region_id: 7,
            local_pos: 0,
            payloads: vec![frame(b"", 112)],
            counts: vec![0],
        })
        .unwrap_err();
    assert!(format!("{err:#}").contains("empty compressed payload frame"));
}

#[test]
fn natural_writer_rejects_multiple_payloads() {
    let tmp = tempdir().unwrap();
    let mut writer = FragmentBankWriter::create(tmp.path(), 1).unwrap();
    let err = writer
        .write_region(PackedRegion {
            region_id: 7,
            local_pos: 0,
            payloads: vec![frame(b"a", 1), frame(b"b", 1)],
            counts: vec![1, 1],
        })
        .unwrap_err();
    assert!(format!("{err:#}").contains("expected 1"));
}

#[test]
fn canonical_config_is_natural_dnapack8_with_10k_store() {
    let cfg = load_from_value(valid_config_json()).1.unwrap();
    assert_eq!(
        cfg.typed.reservoir_sampling.stored_fragments_per_region,
        10_000
    );
    assert_eq!(cfg.typed.storage.format, "fragment bank");
    assert_eq!(cfg.typed.storage.index_format, "bank index");
}
