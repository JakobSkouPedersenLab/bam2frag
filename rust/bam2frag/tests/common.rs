use bam2frag::{
    cli::ConfigArgs,
    config::{load_resolved, ResolvedConfig},
    model::{MaterializedFragment, ReservoirKey},
};
use serde_json::{json, Value};
use std::{fs, path::Path};
use tempfile::TempDir;

pub fn valid_config_json() -> Value {
    json!({
        "coordinate_system": {
            "genome": "hg38",
            "chromosomes": ["chr1", "chr2"],
            "region_window_bp": 1000,
            "region_stride": 1000000
        },
        "reservoir_sampling": {
            "fetch_padding_bp": 500,
            "stored_fragments_per_region": 10000,
            "seed": 7
        },
        "fragment_atac_covariates": {
            "enabled": true,
            "mask_bin_size_bp": 10,
            "groups": [
                {"name": "blood", "members": [
                    "neutrophil__pc", "monocyte__pc", "b_cell__pc", "t_cell__pc",
                    "natural_killer_cell__pc", "erythroblasts__pc", "mega__pc"
                ]},
                {"name": "colon", "members": ["colonic_mucosa__tis"]},
                {"name": "bladder", "members": ["mucosa_of_urinary_bladder__tis"]}
            ]
        },
        "feature_sets": {
            "fragment_sequence": {
                "array_names": [
                    "fragment_ref", "fragment_r1", "fragment_r2",
                    "fragment_r1_qual", "fragment_r2_qual",
                    "fragment_r1_cycle", "fragment_r2_cycle"
                ],
                "seq_window_w": 320
            }
        },
        "filters": {
            "patient_snp_blacklist": {
                "enabled": true
            },
            "donor_blacklist": {
                "enabled": false
            },
            "low_power_blacklist": {
                "enabled": false
            },
            "problematic_blacklist": {
                "enabled": false
            },
            "pair_eligibility": {
                "strict_trustworthy_ends": true,
                "reject_clipping": true,
                "reject_protruding_3p": true,
                "require_full_fragment_coverage": true,
                "require_complete_reference_window": true
            }
        },
        "min_mapq": 30,
        "max_reference_mismatch_fraction": 0.10,
        "min_len": 30,
        "max_len": 299
    })
}

pub fn natural_config_json() -> Value {
    valid_config_json()
}

pub fn natural_resolved_config() -> ResolvedConfig {
    load_from_value(natural_config_json()).1.unwrap()
}

pub fn write_json(path: &Path, value: &Value) {
    fs::write(path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
}

pub fn load_from_value(config: Value) -> (TempDir, anyhow::Result<ResolvedConfig>) {
    let dir = tempfile::tempdir().unwrap();
    let config_path = dir.path().join("fragment_data.json");
    write_json(&config_path, &config);
    let result = load_resolved(&ConfigArgs {
        config: config_path,
    });
    (dir, result)
}

pub fn resolved_config() -> ResolvedConfig {
    load_from_value(valid_config_json()).1.unwrap()
}

/// One 60 bp fragment that passes the full fragment payload validation: two
/// 30-cycle reads meeting in the middle of the centred window, with cycle 1 at
/// both validated 5' endpoints.
pub fn valid_materialized_fragment(atac: [f32; 3], orientation_code: u8) -> MaterializedFragment {
    let width = 320;
    let mut r1_aln = vec![0; width];
    let mut r2_aln = vec![0; width];
    let mut r1_qual = vec![0; width];
    let mut r2_qual = vec![0; width];
    let mut r1_cycle = vec![0; width];
    let mut r2_cycle = vec![0; width];
    let (r1_positions, r2_positions): (Vec<usize>, Vec<usize>) = if orientation_code == 0 {
        ((130..160).collect(), (160..190).rev().collect())
    } else {
        ((160..190).rev().collect(), (130..160).collect())
    };
    for (cycle0, position) in r1_positions.into_iter().enumerate() {
        r1_aln[position] = 2;
        r1_qual[position] = 31;
        r1_cycle[position] = (cycle0 + 1) as u8;
    }
    for (cycle0, position) in r2_positions.into_iter().enumerate() {
        r2_aln[position] = 3;
        r2_qual[position] = 47;
        r2_cycle[position] = (cycle0 + 1) as u8;
    }
    MaterializedFragment {
        key: ReservoirKey { region_id: 7 },
        anchor_rpos0: 10,
        ref_aln: vec![1; width],
        r1_aln,
        r2_aln,
        r1_qual,
        r2_qual,
        r1_cycle,
        r2_cycle,
        fragment_span_bp: 60,
        orientation_code,
        atac_overlap_fractions: atac,
    }
}
