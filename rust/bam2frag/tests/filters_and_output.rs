mod common;

use bam2frag::{
    arrays::{compress_payload, encode_payload, FragmentBankWriter, PackedRegion, PayloadKind},
    cli::InspectArgs,
    filters::{load_patient_variant_blacklists, FilterContext, IntervalSet},
    output::inspect_output,
};
use ndarray::Array1;
use std::{fs, path::Path};
use tempfile::tempdir;

fn write_file(path: &Path, contents: impl AsRef<[u8]>) {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).unwrap();
    }
    fs::write(path, contents).unwrap();
}

fn write_complete_packed_output_layout(output_dir: &Path) {
    write_file(
        &output_dir.join("sample_meta.json"),
        br#"{
          "generator": "bam2frag",
          "sample_id": "S1",
          "patient_id": "P1",
          "n_regions": 1,
          "sampling_mode": "natural",
          "payload_order": ["natural"],
          "max_stored_fragments_per_region": 10000,
          "fragment_bank_file": "fragment_bank.bin",
          "fragment_bank_compression": {
            "codec": "zstd",
            "level": 3,
            "frame_content_size": true,
            "frame_checksum": true
          },
          "fragment_bank_index_file": "fragment_bank_index.bin",
          "seq_window_w": 320,
          "endpoint_definition": "validated_outer_5p",
          "complete_reference_window_required": true,
          "strict_protruding_3p_rejection": true,
          "stored_structural_metadata_fields": ["fragment_span_bp", "fragment_orientation"],
          "fragment_num_feature_names": [
            "atac_blood_overlap_fraction",
            "atac_colon_overlap_fraction",
            "atac_bladder_overlap_fraction"
          ],
          "fragment_atac_mask_bin_size_bp": 10,
          "fragment_atac_source_hash": "fixture-hash",
          "fragment_atac_usage": "optional add-on; zero-filled when fragment ATAC covariates are disabled"
        }"#,
    );

    ndarray_npy::write_npy(
        output_dir.join("region_available_count.npy"),
        &Array1::from(vec![1u64]),
    )
    .unwrap();

    for name in ["extraction_stats.csv", "extraction_stats.json"] {
        write_file(&output_dir.join(name), b"dummy");
    }
    let cfg = common::natural_resolved_config();
    let payload = encode_payload(7, PayloadKind::Natural, &[], &cfg).unwrap();
    write_one_region_bank(output_dir, &payload);
}

/// Replace the layout's bank with one region holding `payload`, compressed.
fn write_one_region_bank(output_dir: &Path, payload: &[u8]) {
    let mut writer = FragmentBankWriter::create(output_dir, 1).unwrap();
    writer
        .write_region(PackedRegion {
            region_id: 7,
            payloads: vec![compress_payload(payload).unwrap()],
            counts: vec![0],
        })
        .unwrap();
    writer.finish(output_dir).unwrap();
}

fn make_layout_natural(output_dir: &Path) {
    write_complete_packed_output_layout(output_dir);
}

#[test]
fn interval_set_uses_half_open_bed_coordinates() {
    let tmp = tempdir().unwrap();
    let bed = tmp.path().join("intervals.bed");
    write_file(&bed, b"chr1\t10\t20\nchr1\t30\t40\nchr2\t5\t8\n");
    let intervals = IntervalSet::load(&bed).unwrap();
    assert!(!intervals.contains("chr1", 9));
    assert!(intervals.contains("chr1", 10));
    assert!(intervals.contains("chr1", 19));
    assert!(!intervals.contains("chr1", 20));
    assert!(intervals.contains("chr2", 5));
    assert!(!intervals.contains("chr2", 8));
}

#[test]
fn filter_context_combines_all_interval_sources() {
    let tmp = tempdir().unwrap();
    let genomic = tmp.path().join("genomic.bed");
    let donor = tmp.path().join("donor.bed");
    write_file(&genomic, b"chr1\t10\t20\n");
    write_file(&donor, b"chr1\t30\t40\n");
    let mut ctx = FilterContext::default();
    ctx.genomic = Some(IntervalSet::load(&genomic).unwrap());
    ctx.donor = Some(IntervalSet::load(&donor).unwrap());
    assert!(ctx.position_filtered("chr1", 10));
    assert!(ctx.position_filtered("chr1", 35));
    assert!(!ctx.position_filtered("chr1", 20));
}

#[test]
fn patient_variant_blacklists_load_bed_coordinates_and_count_exact_hits() {
    let tmp = tempdir().unwrap();
    let base = tmp.path().join("variant_blacklists");
    let patient_dir = base.join("P001");
    fs::create_dir_all(&patient_dir).unwrap();
    write_file(
        &patient_dir.join("P001_snvs_calls.tsv"),
        b"#CHR\tSTART\tREF\tALT\n1\t10\tA\tC\n",
    );
    write_file(
        &patient_dir.join("P001_indels_calls.tsv"),
        b"#CHR\tSTART\tEND\tREF\tALT\n3\t100\t103\tATG\tA\n",
    );
    let (snv, indel) = load_patient_variant_blacklists(&base, "P001").unwrap();
    let mut ctx = FilterContext::default();
    ctx.snv = snv;
    ctx.indel = indel;
    assert!(ctx.snv_blacklisted("chr1", 10, "a", "c"));
    assert!(ctx.indel_blacklisted("chr3", 100, 103, "atg", "a"));
    assert_eq!(ctx.blacklist_hit_counts(), (1, 1));
}

#[test]
fn inspect_output_accepts_complete_natural_v6_10k_layout() {
    let tmp = tempdir().unwrap();
    write_complete_packed_output_layout(tmp.path());
    inspect_output(InspectArgs {
        output_dir: tmp.path().to_path_buf(),
    })
    .unwrap();
}

#[test]
fn inspect_output_rejects_metadata_from_another_generator() {
    let tmp = tempdir().unwrap();
    make_layout_natural(tmp.path());
    let meta_path = tmp.path().join("sample_meta.json");
    let mut meta: serde_json::Value =
        serde_json::from_slice(&fs::read(&meta_path).unwrap()).unwrap();
    meta["generator"] = serde_json::json!("something-else");
    fs::write(&meta_path, serde_json::to_vec_pretty(&meta).unwrap()).unwrap();

    let err = inspect_output(InspectArgs {
        output_dir: tmp.path().to_path_buf(),
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("bam2frag"));
}

#[test]
fn inspect_output_rejects_wrong_endpoint_metadata() {
    let tmp = tempdir().unwrap();
    make_layout_natural(tmp.path());
    let meta_path = tmp.path().join("sample_meta.json");
    let mut meta: serde_json::Value =
        serde_json::from_slice(&fs::read(&meta_path).unwrap()).unwrap();
    meta["endpoint_definition"] = serde_json::json!("outer_cigar");
    fs::write(&meta_path, serde_json::to_vec_pretty(&meta).unwrap()).unwrap();

    let err = inspect_output(InspectArgs {
        output_dir: tmp.path().to_path_buf(),
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("endpoint/reference contract"));
}

#[test]
fn inspect_output_rejects_malformed_payload() {
    // A well-formed frame around a payload that is not a fragment payload: decompression
    // succeeds and the payload validation is what refuses it.
    let tmp = tempdir().unwrap();
    make_layout_natural(tmp.path());
    let cfg = common::natural_resolved_config();
    let mut payload = encode_payload(7, PayloadKind::Natural, &[], &cfg).unwrap();
    payload[0] = b'X';
    write_one_region_bank(tmp.path(), &payload);

    let err = inspect_output(InspectArgs {
        output_dir: tmp.path().to_path_buf(),
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("payload magic"));
}

#[test]
fn inspect_output_rejects_a_corrupted_frame() {
    let tmp = tempdir().unwrap();
    make_layout_natural(tmp.path());
    let bank_path = tmp.path().join("fragment_bank.bin");
    let mut bank = fs::read(&bank_path).unwrap();
    *bank.last_mut().unwrap() ^= 0xff;
    fs::write(&bank_path, bank).unwrap();

    let err = inspect_output(InspectArgs {
        output_dir: tmp.path().to_path_buf(),
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("decompressing"));
}

#[test]
fn inspect_output_rejects_metadata_without_payload_compression() {
    let tmp = tempdir().unwrap();
    make_layout_natural(tmp.path());
    let meta_path = tmp.path().join("sample_meta.json");
    let mut meta: serde_json::Value =
        serde_json::from_slice(&fs::read(&meta_path).unwrap()).unwrap();
    meta.as_object_mut()
        .unwrap()
        .remove("fragment_bank_compression");
    fs::write(&meta_path, serde_json::to_vec_pretty(&meta).unwrap()).unwrap();

    let err = inspect_output(InspectArgs {
        output_dir: tmp.path().to_path_buf(),
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("payload compression"));
}

#[test]
fn inspect_output_rejects_missing_available_count() {
    let tmp = tempdir().unwrap();
    write_complete_packed_output_layout(tmp.path());
    fs::remove_file(tmp.path().join("region_available_count.npy")).unwrap();
    let err = inspect_output(InspectArgs {
        output_dir: tmp.path().to_path_buf(),
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("region_available_count.npy"));
}

#[test]
fn inspect_output_rejects_missing_packed_bank() {
    let tmp = tempdir().unwrap();
    write_complete_packed_output_layout(tmp.path());
    fs::remove_file(tmp.path().join("fragment_bank.bin")).unwrap();
    let err = inspect_output(InspectArgs {
        output_dir: tmp.path().to_path_buf(),
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("fragment_bank.bin"));
}

#[test]
fn inspect_output_rejects_available_count_not_aligned_to_dnaidx8() {
    let tmp = tempdir().unwrap();
    write_complete_packed_output_layout(tmp.path());
    ndarray_npy::write_npy(
        tmp.path().join("region_available_count.npy"),
        &Array1::<u64>::zeros(0),
    )
    .unwrap();
    let err = inspect_output(InspectArgs {
        output_dir: tmp.path().to_path_buf(),
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("does not match bank index row count"));
}

#[test]
fn patient_blacklist_accepts_one_based_pos_headers_and_normalizes_chromosomes() {
    let tmp = tempdir().unwrap();
    let base = tmp.path().join("variant_blacklists");
    let patient_id = "P002";
    let patient_dir = base.join(patient_id);
    fs::create_dir_all(&patient_dir).unwrap();
    write_file(
        &patient_dir.join(format!("{patient_id}_snvs_calls.tsv")),
        b"CHR\tPOS\tREF\tALT\n2\t11\tC\tT\n",
    );
    write_file(
        &patient_dir.join(format!("{patient_id}_indels_calls.tsv")),
        b"CHR\tPOS\tREF\tALT\n2\t21\tGA\tG\n",
    );

    let (snv, indel) = load_patient_variant_blacklists(&base, patient_id).unwrap();
    assert!(snv.contains(&("chr2".to_string(), 10, "C".to_string(), "T".to_string())));
    assert!(indel.contains(&(
        "chr2".to_string(),
        20,
        22,
        "GA".to_string(),
        "G".to_string(),
    )));
}
