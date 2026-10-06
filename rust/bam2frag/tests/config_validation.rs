mod common;

use bam2frag::cli::ConfigArgs;
use bam2frag::config::{load_resolved, ALL_AUTOSOMES, DEFAULT_STORED_FRAGMENTS_PER_REGION};
use common::{load_from_value, valid_config_json};
use serde_json::json;
use std::{fs, path::Path};

#[test]
fn canonical_fixture_stores_10k_fragments_per_region() {
    let cfg = load_from_value(valid_config_json()).1.unwrap();
    assert_eq!(
        cfg.typed.reservoir_sampling.stored_fragments_per_region,
        DEFAULT_STORED_FRAGMENTS_PER_REGION
    );
    assert_eq!(DEFAULT_STORED_FRAGMENTS_PER_REGION, 10_000);
}

#[test]
fn rejects_any_store_size_other_than_10k() {
    for invalid in [1usize, 2048, 4096, 8192, 9999, 10_001] {
        let mut data = valid_config_json();
        data["reservoir_sampling"]["stored_fragments_per_region"] = json!(invalid);
        let error = load_from_value(data).1.unwrap_err();
        let message = format!("{error:#}");
        assert!(
            message.contains("10000") || message.contains("10_000"),
            "invalid={invalid} unexpected error: {message}"
        );
    }
}

#[test]
fn config_applies_documented_defaults() {
    let cfg = load_from_value(valid_config_json()).1.unwrap();
    assert_eq!(cfg.min_mapq, 30);
    assert_eq!(cfg.max_reference_mismatch_fraction, 0.10);
    assert_eq!(cfg.min_len, 30);
    assert_eq!(cfg.max_len, 299);
}

#[test]
fn a_configuration_still_carrying_min_bq_is_rejected() {
    // Extraction applies no base-quality threshold. A configuration written for
    // the old contract must fail loudly rather than have the key ignored, so an
    // extraction cannot quietly run under different rules than its config states.
    let mut data = valid_config_json();
    data["min_bq"] = json!(30);
    let message = format!("{:#}", load_from_value(data).1.unwrap_err());
    assert!(
        message.contains("min_bq"),
        "unexpected error for a config carrying min_bq: {message}"
    );
}

#[test]
fn the_reference_mismatch_burden_cap_must_stay_inside_its_band() {
    // The cap is a structural mismapping guard with no single canonical value,
    // so it is validated as a band: far above the disagreement a real molecule
    // shows, far below a mismapped read's.
    for invalid in [json!(0.0), json!(0.009), json!(0.51), json!(1.0)] {
        let mut data = valid_config_json();
        data["max_reference_mismatch_fraction"] = invalid.clone();
        let message = format!("{:#}", load_from_value(data).1.unwrap_err());
        assert!(
            message.contains("max_reference_mismatch_fraction"),
            "invalid={invalid} unexpected error: {message}"
        );
    }

    for valid in [json!(0.01), json!(0.10), json!(0.50)] {
        let mut data = valid_config_json();
        data["max_reference_mismatch_fraction"] = valid.clone();
        assert!(
            load_from_value(data).1.is_ok(),
            "valid={valid} was rejected"
        );
    }
}

#[test]
fn rejects_noncanonical_pair_eligibility_thresholds_and_switches() {
    for (key, value) in [
        ("min_mapq", json!(29)),
        ("min_len", json!(29)),
        ("max_len", json!(300)),
    ] {
        let mut data = valid_config_json();
        data[key] = value;
        assert!(load_from_value(data)
            .1
            .unwrap_err()
            .to_string()
            .contains("strict pair eligibility"));
    }

    let mut data = valid_config_json();
    data["filters"]["pair_eligibility"]["require_full_fragment_coverage"] = json!(false);
    assert!(load_from_value(data)
        .1
        .unwrap_err()
        .to_string()
        .contains("pair_eligibility"));

    let mut data = valid_config_json();
    data["feature_sets"]["fragment_sequence"]["seq_window_w"] = json!(321);
    assert!(load_from_value(data)
        .1
        .unwrap_err()
        .to_string()
        .contains("seq_window_w must be 320"));
}

#[test]
fn a_config_still_declaring_block_window_scanning_is_rejected() {
    // Regions are always scanned whole. A config carrying the old subsampler's
    // keys must fail rather than have them quietly ignored, which would leave
    // the author believing they had restricted the scan.
    for key in [
        "strategy",
        "blocks_per_region",
        "block_bp",
        "windows_per_block",
        "window_bp",
        "sampling",
        "merge_overlapping_windows",
        "collapse_blocks_to_region_token",
        "enabled",
    ] {
        let mut data = valid_config_json();
        data["reservoir_sampling"][key] = json!("whatever");
        let error = load_from_value(data).1.unwrap_err();
        assert!(
            format!("{error:#}").contains(key),
            "stale scan key {key:?} was accepted: {error:#}"
        );
    }
}

#[test]
fn rejects_removed_molecular_contract_fields() {
    for key in [
        "fragment_event_class",
        "event_class_cardinality",
        "fragment_start_offset_bp",
        "fragment_end_offset_bp",
        "region_width_bp",
        "coarse_offset_bins",
        "fine_offset_bins",
        "span_bins",
        "orientation_cardinality",
        "cut",
        "cut_mask_fraction",
        "ref_context_bp",
        "ref_context_left_bp",
        "ref_context_right_bp",
        "cut_ref_inward_mask_min",
        "cut_ref_inward_mask_max",
        "cut_ref_min_visible_bp",
    ] {
        let mut data = valid_config_json();
        data["removed_probe"] = json!({(key): true});
        assert!(load_from_value(data)
            .1
            .unwrap_err()
            .to_string()
            .contains("removed contract field"));
    }
}

#[test]
fn active_source_and_configs_do_not_contain_removed_event_fields() {
    // The files below are the rejection lists themselves: they must name the
    // removed fields in order to refuse them. Excluding them by path keeps the
    // needles readable as plain literals, so an ordinary repository-wide grep
    // for these names stays meaningful instead of silently reporting "clean".
    const REJECTION_LIST_FILES: [&str; 3] =
        ["config.rs", "config_validation.rs", "config_contract.py"];

    fn inspect_tree(path: &Path, needles: &[&str]) {
        for entry in fs::read_dir(path).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_dir() {
                if path.file_name().and_then(|name| name.to_str()) != Some("target") {
                    inspect_tree(&path, needles);
                }
            } else if matches!(
                path.extension().and_then(|ext| ext.to_str()),
                Some("rs" | "py" | "json")
            ) {
                let name = path
                    .file_name()
                    .and_then(|name| name.to_str())
                    .unwrap_or("");
                if REJECTION_LIST_FILES.contains(&name) {
                    continue;
                }
                let text = fs::read_to_string(&path).unwrap();
                for needle in needles {
                    assert!(
                        !text.contains(needle),
                        "{} contains removed field {needle}",
                        path.display()
                    );
                }
            }
        }
    }

    let repo = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let needles = ["fragment_event_class", "event_class_cardinality"];
    inspect_tree(&repo.join("src"), &needles);
    inspect_tree(&repo.join("rust"), &needles);
}

#[test]
fn cli_exposes_only_the_natural_extraction_endpoint() {
    use bam2frag::cli::Cli;
    use clap::CommandFactory;

    let command = Cli::command();
    let names = command
        .get_subcommands()
        .map(|sub| sub.get_name().to_string())
        .collect::<Vec<_>>();

    assert!(names.contains(&"extract".to_string()));
    assert!(names.contains(&"validate-config".to_string()));
    assert!(names.contains(&"inspect".to_string()));
    // bam2frag extracts natural cfDNA fragments only. The lesion-anchored
    // ctDNA endpoint is deliberately absent.
    assert!(!names.contains(&"extract-ctdna-fragments".to_string()));

    let extract = command
        .get_subcommands()
        .find(|sub| sub.get_name() == "extract")
        .unwrap();
    let args = extract
        .get_arguments()
        .map(|arg| arg.get_id().as_str().to_string())
        .collect::<Vec<_>>();

    assert!(args.contains(&"bam".to_string()));
    assert!(args.contains(&"reference".to_string()));
    assert!(args.contains(&"output_dir".to_string()));
    // Patient and sample covariate tables are not part of bam2frag.
    assert!(!args.contains(&"clinical_info".to_string()));
    assert!(!args.contains(&"qcstats".to_string()));
    assert!(!args.contains(&"sample_metadata".to_string()));
}

/// The `config.json` shipped at the repository root is the contract every user
/// runs with, so it must satisfy the same validator as any other config.
#[test]
fn shipped_repository_config_validates() {
    let config = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../config.json");
    assert!(
        config.is_file(),
        "missing repository config: {}",
        config.display()
    );

    let resolved = load_resolved(&ConfigArgs {
        config: config.clone(),
    })
    .unwrap_or_else(|error| panic!("repository config.json must validate: {error:#}"));

    assert_eq!(resolved.min_mapq, 30);
    assert_eq!(resolved.min_len, 30);
    assert_eq!(resolved.max_len, 299);
    assert_eq!(
        resolved
            .typed
            .reservoir_sampling
            .stored_fragments_per_region,
        DEFAULT_STORED_FRAGMENTS_PER_REGION
    );
    // ATAC covariates are an optional add-on, off in the shipped config so that
    // a plain run needs no chromatin atlas.
    assert!(!resolved.typed.fragment_atac_covariates.enabled);

    // The shipped config uses the "all" shortcut; it must resolve to the autosomes.
    assert_eq!(resolved.typed.coordinate_system.chromosomes, ALL_AUTOSOMES);

    // bam2frag ships with every cohort blacklist switched off, so a plain run
    // needs no patient-specific resources.
    assert!(!resolved.requires_snp_blacklist());
    assert!(!resolved.typed.filters.donor_blacklist.enabled);
    assert!(!resolved.typed.filters.low_power_blacklist.enabled);
    assert!(!resolved.typed.filters.problematic_blacklist.enabled);
}

#[test]
fn chromosomes_all_shortcut_expands_to_the_autosomes() {
    let mut config = valid_config_json();
    config["coordinate_system"]["chromosomes"] = json!("all");

    let cfg = load_from_value(config).1.unwrap();

    assert_eq!(cfg.typed.coordinate_system.chromosomes.len(), 22);
    assert_eq!(cfg.typed.coordinate_system.chromosomes, ALL_AUTOSOMES);
    // The sex chromosomes are never included implicitly.
    assert!(!cfg
        .typed
        .coordinate_system
        .chromosomes
        .iter()
        .any(|chrom| chrom == "chrX" || chrom == "chrY"));
}

#[test]
fn chromosomes_all_shortcut_is_case_insensitive() {
    let mut config = valid_config_json();
    config["coordinate_system"]["chromosomes"] = json!("ALL");

    let cfg = load_from_value(config).1.unwrap();

    assert_eq!(cfg.typed.coordinate_system.chromosomes, ALL_AUTOSOMES);
}

#[test]
fn explicit_chromosome_lists_are_still_accepted() {
    let mut config = valid_config_json();
    config["coordinate_system"]["chromosomes"] = json!(["chr7", "chrX"]);

    let cfg = load_from_value(config).1.unwrap();

    assert_eq!(cfg.typed.coordinate_system.chromosomes, ["chr7", "chrX"]);
}

#[test]
fn unknown_chromosome_shortcut_is_rejected() {
    let mut config = valid_config_json();
    config["coordinate_system"]["chromosomes"] = json!("autosomes");

    let error = load_from_value(config).1.unwrap_err();

    // A silently empty scan would write a valid but empty sample directory.
    assert!(
        format!("{error:#}").contains("unknown chromosome shortcut"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn empty_chromosome_list_is_rejected() {
    let mut config = valid_config_json();
    config["coordinate_system"]["chromosomes"] = json!([]);

    let error = load_from_value(config).1.unwrap_err();

    assert!(
        format!("{error:#}").contains("at least one chromosome"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn duplicate_chromosomes_are_rejected() {
    let mut config = valid_config_json();
    config["coordinate_system"]["chromosomes"] = json!(["chr1", "chr2", "chr1"]);

    let error = load_from_value(config).1.unwrap_err();

    assert!(
        format!("{error:#}").contains("more than once"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn fragment_atac_covariates_may_be_disabled() {
    // The add-on is optional: with it off, no chromatin atlas is needed and the
    // atlas declarations are not checked, so a config may leave them out.
    let mut config = valid_config_json();
    config["fragment_atac_covariates"] = json!({"enabled": false});

    let cfg = load_from_value(config).1.unwrap();

    assert!(!cfg.typed.fragment_atac_covariates.enabled);
}

#[test]
fn disabling_atac_does_not_require_a_matching_bin_size() {
    // mask_bin_size_bp must equal chromatin_bin_size_bp only while enabled.
    let mut config = valid_config_json();
    config["fragment_atac_covariates"] = json!({"enabled": false, "mask_bin_size_bp": 999});

    assert!(load_from_value(config).1.is_ok());
}

#[test]
fn enabled_atac_still_validates_its_atlas_declarations() {
    let mut config = valid_config_json();
    config["fragment_atac_covariates"]["groups"][0]["name"] = json!("plasma");

    let error = load_from_value(config).1.unwrap_err();

    assert!(
        format!("{error:#}").contains("blood, colon, bladder"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn chromosome_batch_size_is_a_cli_flag_with_the_historical_default() {
    use bam2frag::cli::Cli;
    use clap::CommandFactory;

    let extract = Cli::command();
    let extract = extract
        .get_subcommands()
        .find(|sub| sub.get_name() == "extract")
        .unwrap();
    let arg = extract
        .get_arguments()
        .find(|arg| arg.get_id() == "chromosome_batch_size")
        .expect("extract must expose --chromosome-batch-size");

    assert_eq!(arg.get_long(), Some("chromosome-batch-size"));
    // 12 was the hardcoded value before the flag existed; keeping it as the
    // default means adding the flag changed no existing run's behaviour.
    assert_eq!(
        arg.get_default_values(),
        ["12"],
        "default must stay 12 so behaviour is unchanged when the flag is omitted"
    );
}

#[test]
fn a_config_still_declaring_the_fixed_encoding_is_rejected() {
    // These keys used to be validated against one literal and then ignored.
    // Now that they are gone, a config carrying them must fail loudly rather
    // than have them silently dropped, which is what made them worth removing.
    for key in [
        "dna_code",
        "quality_scale",
        "seq_in_channels",
        "window_placement",
        "centred_window_formula",
        "cycle_encoding",
        "cycle_normalization",
        "source_columns",
        "complete_reference_window_required",
    ] {
        let mut config = valid_config_json();
        config["feature_sets"]["fragment_sequence"][key] = json!("whatever");
        let error = load_from_value(config).1.unwrap_err();
        assert!(
            format!("{error:#}").contains(key),
            "stale key {key:?} was accepted: {error:#}"
        );
    }

    for key in [
        "dtype",
        "values",
        "member_aggregation",
        "fragment_aggregation",
        "feature_names",
    ] {
        let mut config = valid_config_json();
        config["fragment_atac_covariates"][key] = json!("whatever");
        let error = load_from_value(config).1.unwrap_err();
        assert!(
            format!("{error:#}").contains(key),
            "stale ATAC key {key:?} was accepted: {error:#}"
        );
    }
}
