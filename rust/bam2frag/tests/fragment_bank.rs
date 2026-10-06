mod common;

use bam2frag::{
    arrays::{
        compress_payload, decode_payload, decompress_payload, encode_payload, payload_nbytes,
        validate_bank_files, FragmentBankWriter, PackedRegion, PayloadKind, PACKED_BANK_FILENAME,
        PACKED_INDEX_FILENAME,
    },
    features::{materialize_validated_into_scratch, MaterializeScratch},
    model::{ReservoirKey, ValidatedFragmentCandidate},
    pair_eligibility::{construct_pair_candidate, PairOrientation},
    reference::ReferenceCache,
};
use common::{natural_resolved_config, valid_materialized_fragment as materialized_fragment};
use rust_htslib::bam::{
    record::{Cigar, CigarString},
    Record,
};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;
use tempfile::TempDir;

/// Every zstd frame starts with this magic number, little-endian 0xFD2FB528.
const ZSTD_FRAME_MAGIC: [u8; 4] = [0x28, 0xB5, 0x2F, 0xFD];
const INDEX_HEADER_BYTES: usize = 16;
const INDEX_ROW_BYTES: usize = 40;

const CHROM_LEN: usize = 2_000;

fn reference_cache() -> (TempDir, ReferenceCache) {
    let temp = tempfile::tempdir().unwrap();
    let fasta = temp.path().join("ref.fa");
    let sequence = "A".repeat(CHROM_LEN);
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
    let reference = ReferenceCache::open(&fasta).unwrap();
    (temp, reference)
}

fn record_with_sequence(
    read1: bool,
    reverse: bool,
    pos0: i64,
    mate_pos0: i64,
    sequence: &[u8],
) -> Record {
    let mut record = Record::new();
    record.set(
        b"centred",
        Some(&CigarString(vec![Cigar::Match(30)])),
        sequence,
        &[35; 30],
    );
    record.set_tid(0);
    record.set_mtid(0);
    record.set_pos(pos0);
    record.set_mpos(mate_pos0);
    record.set_mapq(30);
    record.set_flags(
        0x1 | 0x2
            | if read1 { 0x40 } else { 0x80 }
            | if reverse { 0x10 } else { 0 }
            | if reverse { 0 } else { 0x20 },
    );
    record
}

fn record(read1: bool, reverse: bool, pos0: i64, mate_pos0: i64) -> Record {
    record_with_sequence(read1, reverse, pos0, mate_pos0, &[b'A'; 30])
}

/// `count` valid fragments of one region, encoded as an uncompressed payload.
fn region_payload(region_id: i64, count: usize) -> Vec<u8> {
    let cfg = natural_resolved_config();
    let fragments = (0..count)
        .map(|index| materialized_fragment([0.25, 0.5, 0.75], (index % 2) as u8))
        .collect::<Vec<_>>();
    let refs = fragments.iter().collect::<Vec<_>>();
    encode_payload(region_id, 0, PayloadKind::Natural, &refs, &cfg).unwrap()
}

/// Write `(region_id, payload, count)` regions through the production writer.
fn write_bank(dir: &Path, regions: &[(i64, Vec<u8>, u32)]) -> (PathBuf, PathBuf) {
    let mut writer = FragmentBankWriter::create(dir, regions.len()).unwrap();
    for (region_id, payload, count) in regions {
        writer
            .write_region(PackedRegion {
                region_id: *region_id,
                local_pos: 0,
                payloads: vec![compress_payload(payload).unwrap()],
                counts: vec![*count],
            })
            .unwrap();
    }
    writer.finish(dir).unwrap();
    (
        dir.join(PACKED_BANK_FILENAME),
        dir.join(PACKED_INDEX_FILENAME),
    )
}

fn index_u64(index: &[u8], row: usize, field_offset: usize) -> u64 {
    let start = INDEX_HEADER_BYTES + row * INDEX_ROW_BYTES + field_offset;
    u64::from_le_bytes(index[start..start + 8].try_into().unwrap())
}

// Field offsets inside one bank index row: region index u32, region id i64,
// local position i64, fragment count u32, frame offset u64, stored bytes u64.
const ROW_COUNT: usize = 20;
const ROW_OFFSET: usize = 24;
const ROW_NBYTES: usize = 32;

#[test]
fn headers_round_trip_and_old_or_unknown_versions_are_rejected() {
    let cfg = natural_resolved_config();
    let fragment = materialized_fragment([0.0, 1.0, 0.5], 0);
    let payload = encode_payload(7, 0, PayloadKind::Natural, &[&fragment], &cfg).unwrap();
    assert_eq!(&payload[..8], b"fragment bank");
    assert_eq!(u16::from_le_bytes(payload[8..10].try_into().unwrap()), 8);
    assert_eq!(payload[11], 9);
    assert_eq!(decode_payload(&payload).unwrap().fragment_span_bp, vec![60]);

    let mut old = payload.clone();
    old[8..10].copy_from_slice(&7u16.to_le_bytes());
    assert!(decode_payload(&old)
        .unwrap_err()
        .to_string()
        .contains("version 7"));

    let temp = tempfile::tempdir().unwrap();
    let mut writer = FragmentBankWriter::create(temp.path(), 1).unwrap();
    writer
        .write_region(PackedRegion {
            region_id: 7,
            local_pos: 0,
            payloads: vec![compress_payload(&payload).unwrap()],
            counts: vec![1],
        })
        .unwrap();
    writer.finish(temp.path()).unwrap();
    let index_path = temp.path().join(PACKED_INDEX_FILENAME);
    let index = fs::read(&index_path).unwrap();
    assert_eq!(&index[..8], b"bank index\0");
    assert_eq!(
        validate_bank_files(&temp.path().join(PACKED_BANK_FILENAME), &index_path).unwrap(),
        1
    );
    let mut old_index = index;
    old_index[8..10].copy_from_slice(&7u16.to_le_bytes());
    fs::write(&index_path, old_index).unwrap();
    assert!(
        validate_bank_files(&temp.path().join(PACKED_BANK_FILENAME), &index_path)
            .unwrap_err()
            .to_string()
            .contains("version 7")
    );
}

#[test]
fn stored_frames_decompress_to_the_exact_encoded_payloads_under_an_unchanged_index() {
    let temp = tempfile::tempdir().unwrap();
    let regions = vec![
        (7, region_payload(7, 200), 200),
        (8, region_payload(8, 3), 3),
    ];
    let (bank_path, index_path) = write_bank(temp.path(), &regions);
    assert_eq!(validate_bank_files(&bank_path, &index_path).unwrap(), 2);

    let bank = fs::read(&bank_path).unwrap();
    let index = fs::read(&index_path).unwrap();
    // The index is bank index version 8 with its 40-byte rows; only what the bytes
    // at each offset hold changed.
    assert_eq!(&index[..8], b"bank index\0");
    assert_eq!(u16::from_le_bytes(index[8..10].try_into().unwrap()), 8);
    assert_eq!(
        index.len(),
        INDEX_HEADER_BYTES + regions.len() * INDEX_ROW_BYTES
    );

    let mut expected_offset = 0u64;
    for (row, (_, payload, count)) in regions.iter().enumerate() {
        assert_eq!(index_u64(&index, row, ROW_OFFSET), expected_offset);
        let nbytes = index_u64(&index, row, ROW_NBYTES);
        let frame = &bank[expected_offset as usize..(expected_offset + nbytes) as usize];
        assert_eq!(frame[..4], ZSTD_FRAME_MAGIC);
        assert_eq!(payload_nbytes(*count as usize).unwrap(), payload.len());
        assert_eq!(
            decompress_payload(frame, payload.len()).unwrap(),
            *payload,
            "row {row} did not decompress to its encoded payload"
        );
        expected_offset += nbytes;
    }
    assert_eq!(expected_offset as usize, bank.len());
    // Zero-padded windows and cycle ramps are what make the bank compressible.
    assert!(index_u64(&index, 0, ROW_NBYTES) * 4 < regions[0].1.len() as u64);
}

#[test]
fn the_same_regions_always_give_byte_identical_bank_and_index_files() {
    let regions = vec![
        (7, region_payload(7, 50), 50),
        (9, region_payload(9, 11), 11),
    ];
    let first = tempfile::tempdir().unwrap();
    let second = tempfile::tempdir().unwrap();
    let (first_bank, first_index) = write_bank(first.path(), &regions);
    let (second_bank, second_index) = write_bank(second.path(), &regions);
    assert_eq!(
        fs::read(first_bank).unwrap(),
        fs::read(second_bank).unwrap()
    );
    assert_eq!(
        fs::read(first_index).unwrap(),
        fs::read(second_index).unwrap()
    );
}

/// Write a one-region bank, apply `mutate` to its bank and index bytes, and
/// return the validation error.
fn validation_error(mutate: impl FnOnce(&mut Vec<u8>, &mut Vec<u8>)) -> String {
    let temp = tempfile::tempdir().unwrap();
    let (bank_path, index_path) = write_bank(temp.path(), &[(7, region_payload(7, 20), 20)]);
    let mut bank = fs::read(&bank_path).unwrap();
    let mut index = fs::read(&index_path).unwrap();
    mutate(&mut bank, &mut index);
    fs::write(&bank_path, bank).unwrap();
    fs::write(&index_path, index).unwrap();
    format!(
        "{:#}",
        validate_bank_files(&bank_path, &index_path).unwrap_err()
    )
}

fn set_index_u64(index: &mut [u8], field_offset: usize, value: u64) {
    let start = INDEX_HEADER_BYTES + field_offset;
    index[start..start + 8].copy_from_slice(&value.to_le_bytes());
}

#[test]
fn a_corrupt_truncated_or_mis_indexed_frame_is_rejected() {
    // The last four bytes of the frame are its checksum.
    let error = validation_error(|bank, _| *bank.last_mut().unwrap() ^= 0xff);
    assert!(error.contains("decompressing"), "{error}");

    // A frame cut short while the index still claims its full length.
    let error = validation_error(|bank, _| {
        bank.pop();
    });
    assert!(error.contains("outside the packed bank"), "{error}");

    // A frame cut short with the index agreeing: no longer one complete frame.
    let error = validation_error(|bank, index| {
        bank.pop();
        set_index_u64(index, ROW_NBYTES, bank.len() as u64);
    });
    assert!(error.contains("not a zstd frame"), "{error}");

    // A second frame inside one row's stored bytes.
    let error = validation_error(|bank, index| {
        let frame = bank.clone();
        bank.extend_from_slice(&frame);
        set_index_u64(index, ROW_NBYTES, bank.len() as u64);
    });
    assert!(error.contains("followed by"), "{error}");

    // Bytes after the last indexed frame.
    let error = validation_error(|bank, _| bank.extend_from_slice(b"junk"));
    assert!(error.contains("trailing bytes"), "{error}");

    // A row whose fragment count does not match the payload its frame holds.
    let error = validation_error(|_, index| {
        let start = INDEX_HEADER_BYTES + ROW_COUNT;
        index[start..start + 4].copy_from_slice(&21u32.to_le_bytes());
    });
    assert!(error.contains("fragment count implies"), "{error}");

    // An uncompressed payload, as every bank written before compression holds.
    let error = validation_error(|bank, index| {
        *bank = region_payload(7, 20);
        set_index_u64(index, ROW_NBYTES, bank.len() as u64);
    });
    assert!(error.contains("uncompressed fragment payload"), "{error}");
}

#[test]
fn decompression_requires_one_frame_of_exactly_the_expected_length() {
    let payload = region_payload(7, 5);
    let compressed = compress_payload(&payload).unwrap();
    assert_eq!(compressed.payload_nbytes, payload.len() as u64);
    assert_eq!(
        decompress_payload(&compressed.frame, payload.len()).unwrap(),
        payload
    );
    assert!(decompress_payload(&compressed.frame, payload.len() + 1)
        .unwrap_err()
        .to_string()
        .contains("declares"));
    assert!(decompress_payload(b"not-a-frame", payload.len())
        .unwrap_err()
        .to_string()
        .contains("not a zstd frame"));
}

#[test]
fn exact_round_trip_preserves_planes_structural_fields_and_atac_order() {
    let cfg = natural_resolved_config();
    let fragments = vec![
        materialized_fragment([0.0, 0.0, 0.0], 0),
        materialized_fragment([1.0, 1.0, 1.0], 1),
        materialized_fragment([0.125, 0.5, 0.875], 0),
    ];
    let refs = fragments.iter().collect::<Vec<_>>();
    let decoded =
        decode_payload(&encode_payload(7, 0, PayloadKind::Natural, &refs, &cfg).unwrap()).unwrap();
    assert_eq!(
        decoded.atac_overlap_fractions,
        vec![[0.0; 3], [1.0; 3], [0.125, 0.5, 0.875]]
    );
    assert_eq!(decoded.fragment_span_bp, vec![60; 3]);
    assert_eq!(decoded.orientation_code, vec![0, 1, 0]);
    assert_eq!(
        decoded.ref_aln,
        fragments
            .iter()
            .flat_map(|row| row.ref_aln.clone())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        decoded.r1_qual,
        fragments
            .iter()
            .flat_map(|row| row.r1_qual.clone())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        decoded.r2_cycle,
        fragments
            .iter()
            .flat_map(|row| row.r2_cycle.clone())
            .collect::<Vec<_>>()
    );
}

#[test]
fn validated_endpoints_are_centred_and_cycles_start_at_both_five_prime_ends() {
    let cfg = natural_resolved_config();
    let (temp, mut reference) = reference_cache();
    let read1 = record(true, false, 500, 530);
    let read2 = record(false, true, 530, 500);
    let pair = construct_pair_candidate(
        &read1,
        &read2,
        "chr1",
        CHROM_LEN as i64,
        &cfg,
        &mut reference,
    )
    .unwrap()
    .unwrap();
    assert_eq!(pair.orientation, PairOrientation::R1ForwardR2Reverse);
    assert_eq!(pair.fragment_span_bp, 60);
    assert_eq!(pair.left_5p_genomic - pair.window_start0, 130);
    assert_eq!(pair.right_5p_genomic - pair.window_start0, 189);

    let mut scratch = MaterializeScratch::with_width(320);
    let candidate = ValidatedFragmentCandidate {
        key: ReservoirKey {
            region_id: 7,
            local_pos: 0,
        },
        anchor_rpos0: 529,
        pair,
    };
    let mut second_reference = ReferenceCache::open(&temp.path().join("ref.fa")).unwrap();
    let scalars = materialize_validated_into_scratch(
        candidate,
        &cfg,
        &mut second_reference,
        &mut scratch,
        None,
    )
    .unwrap();
    assert_eq!(scalars.fragment_span_bp, 60);
    assert_eq!(scalars.orientation_code, 0);
    assert_eq!(scratch.ref_aln, vec![1; 320]);
    assert_eq!(scratch.r1_cycle[130], 1);
    assert_eq!(scratch.r1_cycle[159], 30);
    assert_eq!(scratch.r2_cycle[160], 30);
    assert_eq!(scratch.r2_cycle[189], 1);
    assert!(scratch.r1_cycle[..130].iter().all(|cycle| *cycle == 0));
    assert!(scratch.r2_cycle[190..].iter().all(|cycle| *cycle == 0));
}

#[test]
fn accepted_snv_is_only_a_ref_versus_read_difference() {
    let cfg = natural_resolved_config();
    let (temp, mut reference) = reference_cache();
    let mut forward = [b'A'; 30];
    let mut reverse = [b'A'; 30];
    forward[25] = b'C';
    reverse[5] = b'C';
    let read1 = record_with_sequence(true, false, 500, 520, &forward);
    let read2 = record_with_sequence(false, true, 520, 500, &reverse);
    let pair = construct_pair_candidate(
        &read1,
        &read2,
        "chr1",
        CHROM_LEN as i64,
        &cfg,
        &mut reference,
    )
    .unwrap()
    .unwrap();
    let mut scratch = MaterializeScratch::with_width(320);
    let candidate = ValidatedFragmentCandidate {
        key: ReservoirKey {
            region_id: 7,
            local_pos: 0,
        },
        anchor_rpos0: 524,
        pair,
    };
    let mut second_reference = ReferenceCache::open(&temp.path().join("ref.fa")).unwrap();
    materialize_validated_into_scratch(candidate, &cfg, &mut second_reference, &mut scratch, None)
        .unwrap();
    let snv_idx = (320 - 50) / 2 + 25;
    assert_eq!(scratch.ref_aln[snv_idx], 1);
    assert_eq!(scratch.r1_aln[snv_idx], 2);
    assert_eq!(scratch.r2_aln[snv_idx], 2);
}

#[test]
fn malformed_plane_table_payload_and_index_are_rejected() {
    let cfg = natural_resolved_config();
    let fragment = materialized_fragment([0.0; 3], 0);
    let mut payload = encode_payload(7, 0, PayloadKind::Natural, &[&fragment], &cfg).unwrap();
    let mut malformed_structural = payload.clone();
    let orientation_offset = malformed_structural.len() - 2;
    malformed_structural[orientation_offset] = 2;
    assert!(decode_payload(&malformed_structural)
        .unwrap_err()
        .to_string()
        .contains("structural record"));
    payload[40..48].copy_from_slice(&1u64.to_le_bytes());
    assert!(decode_payload(&payload)
        .unwrap_err()
        .to_string()
        .contains("plane sizes"));

    let temp = tempfile::tempdir().unwrap();
    fs::write(temp.path().join(PACKED_BANK_FILENAME), b"short").unwrap();
    fs::write(temp.path().join(PACKED_INDEX_FILENAME), b"bank index\0").unwrap();
    assert!(validate_bank_files(
        &temp.path().join(PACKED_BANK_FILENAME),
        &temp.path().join(PACKED_INDEX_FILENAME)
    )
    .is_err());
}

#[test]
#[ignore = "measurement-only benchmark; run explicitly on the Linux extraction host"]
fn reports_dnapack8_decode_throughput() {
    let cfg = natural_resolved_config();
    let fragments = (0..10_000)
        .map(|index| materialized_fragment([0.25, 0.5, 0.75], (index % 2) as u8))
        .collect::<Vec<_>>();
    let refs = fragments.iter().collect::<Vec<_>>();
    let payload = encode_payload(7, 0, PayloadKind::Natural, &refs, &cfg).unwrap();

    let iterations = 20usize;
    let started = Instant::now();
    for _ in 0..iterations {
        let decoded = decode_payload(&payload).unwrap();
        assert_eq!(decoded.fragment_span_bp.len(), fragments.len());
    }
    let elapsed = started.elapsed().as_secs_f64();
    let fragments_per_second = (iterations * fragments.len()) as f64 / elapsed;
    eprintln!(
        "fragment payload_bytes={} fragments_per_second={fragments_per_second:.0}",
        payload.len()
    );
}

#[test]
#[ignore = "measurement-only benchmark; run explicitly on the Linux extraction host"]
fn reports_payload_compression_ratio_and_throughput() {
    // Synthetic fragments are far more uniform than real ones, so the ratio here
    // overstates a real bank's; the speeds are what this measures.
    let payload = region_payload(7, 10_000);
    let iterations = 10usize;
    let started = Instant::now();
    let mut compressed = compress_payload(&payload).unwrap();
    for _ in 1..iterations {
        compressed = compress_payload(&payload).unwrap();
    }
    let compress_seconds = started.elapsed().as_secs_f64();
    let started = Instant::now();
    for _ in 0..iterations {
        let restored = decompress_payload(&compressed.frame, payload.len()).unwrap();
        assert_eq!(restored.len(), payload.len());
    }
    let decompress_seconds = started.elapsed().as_secs_f64();
    let megabytes = (iterations * payload.len()) as f64 / 1e6;
    eprintln!(
        "zstd payload_bytes={} stored_bytes={} ratio={:.1} compress_mb_per_s={:.0} \
         decompress_mb_per_s={:.0}",
        payload.len(),
        compressed.frame.len(),
        payload.len() as f64 / compressed.frame.len() as f64,
        megabytes / compress_seconds,
        megabytes / decompress_seconds,
    );
}
