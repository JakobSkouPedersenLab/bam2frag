use crate::{
    arrays::{
        bank_compression_meta, ArrayMeta, BankCompressionMeta, BANK_COMPRESSION_CODEC,
        BANK_COMPRESSION_LEVEL, CYCLE_NORMALIZATION, CYCLE_PLANE_ENCODING, ENDPOINT_DEFINITION,
        STORED_STRUCTURAL_METADATA_FIELDS,
    },
    atac::ATAC_TISSUE_FEATURE_NAMES,
    cli::InspectArgs,
    config::{AtacTissueGroupConfig, ResolvedConfig},
    model::ChromStats,
    pair_eligibility::{FragmentObservationStats, PairRejectionCounters, PAIR_COUNTER_FIELDS},
};
use anyhow::{bail, Result};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs::{self, File},
    path::Path,
};

pub const EXTRACTION_STATS_FORMAT: &str = "bam2frag_extraction_stats";

/// One `extraction_stats.json` row: the base counts plus every named
/// pair-eligibility counter flattened alongside them.
#[derive(Serialize)]
struct StatsRow<'a> {
    chrom: &'a str,
    representative_reads: u64,
    skipped: u64,
    evicted: u64,
    invalid_fragments: u64,
    #[serde(flatten)]
    pair_rejections: &'a PairRejectionCounters,
    /// Span- and quality-stratified histograms. JSON only: these are arrays,
    /// and `extraction_stats.csv` carries one scalar per column.
    observations: &'a FragmentObservationStats,
}

impl<'a> StatsRow<'a> {
    fn new(row: &'a ChromStats) -> Self {
        Self {
            chrom: &row.chrom,
            representative_reads: row.reads,
            skipped: row.skipped,
            evicted: row.evicted,
            invalid_fragments: row.invalid_fragments,
            pair_rejections: &row.pair_rejections,
            observations: &row.observations,
        }
    }
}

#[derive(Serialize)]
struct ExtractionStats<'a> {
    format: &'static str,
    sample_id: &'a str,
    patient_id: &'a str,
    total: StatsRow<'a>,
    chromosomes: Vec<StatsRow<'a>>,
}

/// Write per-chromosome and total extraction statistics as CSV and JSON.
///
/// Both files carry the same columns under the same names, and the trailing
/// `total` row/object is the deterministic sum of the per-chromosome rows.
pub fn write_stats(
    csv_path: &Path,
    json_path: &Path,
    rows: &[ChromStats],
    sample_id: &str,
    patient_id: &str,
) -> Result<()> {
    let mut total = ChromStats::new("total");
    for row in rows {
        total.merge(row);
    }

    let mut writer = csv::Writer::from_path(csv_path)?;
    let mut header: Vec<&str> = vec![
        "chrom",
        "representative_reads",
        "skipped",
        "evicted",
        "invalid_fragments",
    ];
    header.extend(PAIR_COUNTER_FIELDS);
    header.extend(["patient_id", "sample_id"]);
    writer.write_record(&header)?;
    for row in rows.iter().chain(std::iter::once(&total)) {
        write_stat(&mut writer, row, sample_id, patient_id)?;
    }
    writer.flush()?;

    let stats = ExtractionStats {
        format: EXTRACTION_STATS_FORMAT,
        sample_id,
        patient_id,
        total: StatsRow::new(&total),
        chromosomes: rows.iter().map(StatsRow::new).collect(),
    };
    fs::write(json_path, serde_json::to_vec_pretty(&stats)?)?;
    Ok(())
}

fn write_stat(
    writer: &mut csv::Writer<File>,
    row: &ChromStats,
    sample_id: &str,
    patient_id: &str,
) -> Result<()> {
    let mut record: Vec<String> = vec![
        row.chrom.clone(),
        row.reads.to_string(),
        row.skipped.to_string(),
        row.evicted.to_string(),
        row.invalid_fragments.to_string(),
    ];
    record.extend(
        row.pair_rejections
            .as_named_counts()
            .iter()
            .map(|(_, value)| value.to_string()),
    );
    record.push(patient_id.to_string());
    record.push(sample_id.to_string());
    writer.write_record(&record)?;
    Ok(())
}

#[derive(Serialize)]
struct Meta<'a> {
    generator: &'static str,
    sample_id: &'a str,
    patient_id: &'a str,
    n_regions: usize,
    sampling_mode: &'static str,
    payload_order: [&'static str; 1],
    max_stored_fragments_per_region: usize,
    fragment_bank_file: String,
    fragment_bank_compression: BankCompressionMeta,
    fragment_bank_index_file: String,
    region_available_count_file: &'static str,
    region_available_count_alignment: &'static str,
    seq_window_w: usize,
    centred_window_formula: &'static str,
    centred_window_formula_version: &'static str,
    complete_reference_window_required: bool,
    endpoint_definition: &'static str,
    strict_protruding_3p_rejection: bool,
    cycle_plane_encoding: &'static str,
    cycle_normalization: &'static str,
    min_len: i64,
    max_len: i64,
    max_reference_mismatch_fraction: f64,
    min_mapq: u8,
    strict_trustworthy_ends: bool,
    reject_clipping: bool,
    reject_protruding_3p: bool,
    require_full_fragment_coverage: bool,
    require_complete_reference_window: bool,
    patient_snp_blacklist_enabled: bool,
    donor_blacklist_enabled: bool,
    low_power_blacklist_enabled: bool,
    problematic_blacklist_enabled: bool,
    stored_structural_metadata_fields: [&'static str; 2],
    removed_molecular_planes: [&'static str; 3],
    region_stride: u64,
    seq_feature_names: Vec<String>,
    fragment_num_feature_names: Vec<String>,
    fragment_atac_covariates_enabled: bool,
    fragment_atac_mask_bin_size_bp: Option<u64>,
    fragment_atac_groups: Option<&'a [AtacTissueGroupConfig]>,
    fragment_atac_atlas_provenance: &'static str,
    fragment_atac_source_hash: String,
    fragment_atac_usage: &'static str,
    array_files: BTreeMap<String, String>,
    summary_files: BTreeMap<String, String>,
    n_sampled_rows_input: usize,
    n_sampled_rows_kept: usize,
    n_sampled_rows_dropped_overflow: usize,
    config_sha256: String,
    implementation: &'static str,
    version: &'static str,
}

pub fn write_meta(
    path: &Path,
    sample: &str,
    patient: &str,
    cfg: &ResolvedConfig,
    arrays: &ArrayMeta,
) -> Result<()> {
    let config_sha256 = format!("{:x}", Sha256::digest(serde_json::to_vec(&cfg.raw)?));
    let array_files = [("region_available_count", "region_available_count.npy")]
        .into_iter()
        .map(|(key, value)| (key.to_string(), value.to_string()))
        .collect();
    let summary_files = [
        ("extraction_stats", "extraction_stats.csv"),
        ("extraction_stats_json", "extraction_stats.json"),
        ("fragment_bank_index", "fragment_bank_index.bin"),
        ("fragment_bank", "fragment_bank.bin"),
    ]
    .into_iter()
    .map(|(key, value)| (key.to_string(), value.to_string()))
    .collect();
    let fragment_atac_source_hash = format!(
        "{:x}",
        Sha256::digest(serde_json::to_vec(&(
            cfg.typed.fragment_atac_covariates.mask_bin_size_bp,
            &cfg.typed.fragment_atac_covariates.groups,
        ))?)
    );
    let meta = Meta {
        generator: "bam2frag",
        sample_id: sample, patient_id: patient, n_regions: arrays.n_regions,
        sampling_mode: "natural", payload_order: ["natural"],
        max_stored_fragments_per_region: cfg.typed.reservoir_sampling.stored_fragments_per_region,
        fragment_bank_file: arrays.fragment_bank_file.clone(),
        fragment_bank_compression: bank_compression_meta(),
        fragment_bank_index_file: arrays.fragment_bank_index_file.clone(),
        region_available_count_file: "region_available_count.npy",
        region_available_count_alignment: "row i corresponds exactly to bank index row i",
        seq_window_w: cfg.typed.feature_sets.fragment_sequence.seq_window_w,
        centred_window_formula: "left_cut_idx=floor((320-fragment_span_bp)/2); right_cut_idx=left_cut_idx+fragment_span_bp-1",
        centred_window_formula_version: "centred_validated_outer_5p",
        complete_reference_window_required: true, endpoint_definition: ENDPOINT_DEFINITION,
        strict_protruding_3p_rejection: true,
        cycle_plane_encoding: CYCLE_PLANE_ENCODING,
        cycle_normalization: CYCLE_NORMALIZATION,
        min_len: cfg.min_len,
        max_len: cfg.max_len,
        max_reference_mismatch_fraction: cfg.max_reference_mismatch_fraction,
        min_mapq: cfg.min_mapq,
        strict_trustworthy_ends: cfg.typed.filters.pair_eligibility.strict_trustworthy_ends,
        reject_clipping: cfg.typed.filters.pair_eligibility.reject_clipping,
        reject_protruding_3p: cfg.typed.filters.pair_eligibility.reject_protruding_3p,
        require_full_fragment_coverage: cfg.typed.filters.pair_eligibility.require_full_fragment_coverage,
        require_complete_reference_window: cfg.typed.filters.pair_eligibility.require_complete_reference_window,
        patient_snp_blacklist_enabled: cfg.typed.filters.patient_snp_blacklist.enabled,
        donor_blacklist_enabled: cfg.typed.filters.donor_blacklist.enabled,
        low_power_blacklist_enabled: cfg.typed.filters.low_power_blacklist.enabled,
        problematic_blacklist_enabled: cfg.typed.filters.problematic_blacklist.enabled,
        stored_structural_metadata_fields: STORED_STRUCTURAL_METADATA_FIELDS,
        removed_molecular_planes: ["event_class", "genomic_start_offset", "genomic_end_offset"],
        region_stride: cfg.typed.coordinate_system.region_stride,
        seq_feature_names: cfg.typed.feature_sets.fragment_sequence.array_names.clone(),
        fragment_num_feature_names: ATAC_TISSUE_FEATURE_NAMES.iter().map(|name| (*name).to_string()).collect(),
        fragment_atac_covariates_enabled: cfg.typed.fragment_atac_covariates.enabled,
        fragment_atac_mask_bin_size_bp: cfg.typed.fragment_atac_covariates.enabled.then_some(cfg.typed.fragment_atac_covariates.mask_bin_size_bp),
        fragment_atac_groups: cfg.typed.fragment_atac_covariates.enabled.then_some(cfg.typed.fragment_atac_covariates.groups.as_slice()),
        fragment_atac_atlas_provenance: "configured aggregate tissue ATAC group definition",
        fragment_atac_source_hash,
        fragment_atac_usage: "optional add-on; zero-filled when fragment ATAC covariates are disabled",
        array_files, summary_files,
        n_sampled_rows_input: arrays.kept + arrays.dropped,
        n_sampled_rows_kept: arrays.kept,
        n_sampled_rows_dropped_overflow: arrays.dropped,
        config_sha256, implementation: "rust", version: env!("CARGO_PKG_VERSION"),
    };
    fs::write(path, serde_json::to_vec_pretty(&meta)?)?;
    Ok(())
}

pub fn inspect_output(args: InspectArgs) -> Result<()> {
    for name in [
        "sample_meta.json",
        "region_available_count.npy",
        "fragment_bank.bin",
        "fragment_bank_index.bin",
        "extraction_stats.csv",
        "extraction_stats.json",
    ] {
        if !args.output_dir.join(name).is_file() {
            bail!("missing required output {name}");
        }
    }
    let available: ndarray::Array1<u64> =
        ndarray_npy::read_npy(args.output_dir.join("region_available_count.npy"))?;
    let index_rows = crate::arrays::validate_bank_files(
        &args.output_dir.join("fragment_bank.bin"),
        &args.output_dir.join("fragment_bank_index.bin"),
    )?;
    if index_rows != available.len() {
        bail!(
            "region_available_count length {} does not match bank index row count {index_rows}",
            available.len()
        );
    }
    if available.iter().any(|count| *count == 0) {
        bail!(
            "region_available_count contains an active bank index row with zero eligible fragments"
        );
    }
    let meta = fs::read_to_string(args.output_dir.join("sample_meta.json"))?;
    let value: serde_json::Value = serde_json::from_str(&meta)?;
    if value.get("generator").and_then(|v| v.as_str()) != Some("bam2frag")
        || value.get("sampling_mode").and_then(|v| v.as_str()) != Some("natural")
        || value
            .get("max_stored_fragments_per_region")
            .and_then(|v| v.as_u64())
            .is_none()
    {
        bail!("sample_meta.json was not written by a natural bam2frag extraction");
    }
    if value.get("seq_window_w").and_then(|v| v.as_u64()) != Some(320)
        || value.get("endpoint_definition").and_then(|v| v.as_str()) != Some(ENDPOINT_DEFINITION)
        || value
            .get("complete_reference_window_required")
            .and_then(|v| v.as_bool())
            != Some(true)
        || value
            .get("strict_protruding_3p_rejection")
            .and_then(|v| v.as_bool())
            != Some(true)
    {
        bail!("sample metadata does not declare the centred endpoint/reference contract");
    }
    let compression = value.get("fragment_bank_compression");
    if compression
        .and_then(|c| c.get("codec"))
        .and_then(|v| v.as_str())
        != Some(BANK_COMPRESSION_CODEC)
        || compression
            .and_then(|c| c.get("level"))
            .and_then(|v| v.as_i64())
            != Some(i64::from(BANK_COMPRESSION_LEVEL))
    {
        bail!(
            "sample metadata does not declare zstd level-{BANK_COMPRESSION_LEVEL} payload \
             compression"
        );
    }
    println!("{meta}");
    Ok(())
}
