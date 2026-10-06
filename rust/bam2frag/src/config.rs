use crate::cli::ConfigArgs;
use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::{fs, path::Path};

pub const DEFAULT_STORED_FRAGMENTS_PER_REGION: usize = 10_000;
pub const MOLECULAR_WINDOW_BP: usize = 320;
const ATAC_TISSUE_GROUP_COUNT: usize = 3;

/// The autosomes, in karyotype order, that the `"all"` shortcut expands to.
///
/// The sex chromosomes are deliberately absent. Their copy number depends on
/// donor sex, so including them implicitly would change what a region means
/// between samples; a run that wants them must name them explicitly.
pub const ALL_AUTOSOMES: [&str; 22] = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
    "chr22",
];

fn all_autosomes() -> Vec<String> {
    ALL_AUTOSOMES
        .iter()
        .map(|chrom| (*chrom).to_string())
        .collect()
}

/// Accept either an explicit chromosome list or the `"all"` shortcut.
///
/// Expansion happens here, at the edge, so every consumer downstream still sees
/// a concrete list and no caller has to know the shortcut exists.
fn deserialize_chromosomes<'de, D>(deserializer: D) -> std::result::Result<Vec<String>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    struct ChromosomesVisitor;

    impl<'de> serde::de::Visitor<'de> for ChromosomesVisitor {
        type Value = Vec<String>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a list of chromosome names, or the shortcut \"all\"")
        }

        fn visit_str<E>(self, value: &str) -> std::result::Result<Self::Value, E>
        where
            E: serde::de::Error,
        {
            if value.eq_ignore_ascii_case("all") {
                return Ok(all_autosomes());
            }
            Err(E::custom(format!(
                "unknown chromosome shortcut {value:?}; use \"all\" for the 22 autosomes, \
                 or list the chromosomes explicitly"
            )))
        }

        fn visit_seq<A>(self, mut seq: A) -> std::result::Result<Self::Value, A::Error>
        where
            A: serde::de::SeqAccess<'de>,
        {
            let mut chromosomes = Vec::with_capacity(seq.size_hint().unwrap_or(0));
            while let Some(chrom) = seq.next_element::<String>()? {
                chromosomes.push(chrom);
            }
            Ok(chromosomes)
        }
    }

    deserializer.deserialize_any(ChromosomesVisitor)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoordinateSystem {
    pub genome: String,
    /// Chromosomes to extract. Written either as an explicit list or as the
    /// string `"all"`, which expands to [`ALL_AUTOSOMES`] while loading.
    #[serde(deserialize_with = "deserialize_chromosomes")]
    pub chromosomes: Vec<String>,
    pub region_window_bp: u64,
    pub region_stride: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlacklistConfig {
    #[serde(default)]
    pub enabled: bool,
}

impl Default for BlacklistConfig {
    fn default() -> Self {
        Self { enabled: false }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FiltersConfig {
    #[serde(default)]
    pub patient_snp_blacklist: BlacklistConfig,
    #[serde(default)]
    pub donor_blacklist: BlacklistConfig,
    #[serde(default)]
    pub low_power_blacklist: BlacklistConfig,
    #[serde(default)]
    pub problematic_blacklist: BlacklistConfig,
    pub pair_eligibility: PairEligibilityConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PairEligibilityConfig {
    pub strict_trustworthy_ends: bool,
    pub reject_clipping: bool,
    pub reject_protruding_3p: bool,
    pub require_full_fragment_coverage: bool,
    pub require_complete_reference_window: bool,
}

/// Uniform reservoir sampling over every eligible fragment in a region.
///
/// Always on: it is what the tool does. The block/window subsampler this block
/// used to configure is gone, so the only settings left are how much padding to
/// fetch, how many fragments to keep, and the seed.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReservoirSamplingConfig {
    #[serde(default = "default_fetch_padding")]
    pub fetch_padding_bp: u64,
    pub stored_fragments_per_region: usize,
    #[serde(default)]
    pub seed: u64,
    #[serde(default)]
    pub max_region_keys: Option<usize>,
}

/// The two sequence settings the extractor actually reads.
///
/// The encoding this config used to restate — 4-bit DNA codes, 6-bit
/// qualities, the centred-window formula, the cycle-plane encoding — has
/// exactly one form, fixed in `arrays.rs`. Declaring it here meant a value the
/// validator compared against a literal and then ignored, so it is described
/// in the README instead and no longer configurable-looking.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FragmentSequenceConfig {
    pub array_names: Vec<String>,
    pub seq_window_w: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AtacTissueGroupConfig {
    pub name: String,
    pub members: Vec<String>,
}

/// The optional ATAC add-on: whether it runs, at what bin size, over which
/// atlas groups.
///
/// The dtype, value range and aggregation rules it used to declare are fixed
/// properties of the covariate plane, not choices, so they live in the README
/// and in `atac.rs` rather than in a config field that could only hold one
/// value. Feature names come from [`ATAC_TISSUE_FEATURE_NAMES`].
///
/// [`ATAC_TISSUE_FEATURE_NAMES`]: crate::atac::ATAC_TISSUE_FEATURE_NAMES
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FragmentAtacCovariatesConfig {
    #[serde(default)]
    pub enabled: bool,
    #[serde(default = "default_chromatin_bin_size_bp")]
    pub mask_bin_size_bp: u64,
    #[serde(default)]
    pub groups: Vec<AtacTissueGroupConfig>,
}

impl Default for FragmentAtacCovariatesConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            mask_bin_size_bp: default_chromatin_bin_size_bp(),
            groups: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureSetsConfig {
    pub fragment_sequence: FragmentSequenceConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RootConfig {
    pub coordinate_system: CoordinateSystem,
    pub reservoir_sampling: ReservoirSamplingConfig,
    pub feature_sets: FeatureSetsConfig,
    #[serde(default)]
    pub fragment_atac_covariates: FragmentAtacCovariatesConfig,
    pub filters: FiltersConfig,
    pub min_mapq: u8,
    pub min_len: i64,
    pub max_len: i64,
    /// Largest fraction of the validated span that may disagree with hg38.
    ///
    /// A structural mismapping guard, not a base-quality threshold. Sequence
    /// quality is carried into the fragment bank on the two quality planes
    /// instead of being thresholded away at extraction.
    pub max_reference_mismatch_fraction: f64,
    #[serde(default)]
    pub max_mate_buffer: Option<usize>,
    #[serde(default)]
    pub chromatin_bin_size_bp: Option<u64>,
}

#[derive(Debug, Clone)]
pub struct ResolvedConfig {
    pub raw: Value,
    pub typed: RootConfig,
    pub min_mapq: u8,
    pub min_len: i64,
    pub max_len: i64,
    pub max_reference_mismatch_fraction: f64,
    pub max_mate_buffer: usize,
    pub chromatin_bin_size_bp: u64,
}

fn default_fetch_padding() -> u64 {
    500
}
fn default_chromatin_bin_size_bp() -> u64 {
    10
}
impl ResolvedConfig {
    pub fn requires_snp_blacklist(&self) -> bool {
        self.typed.filters.patient_snp_blacklist.enabled
    }
}

fn load_json(path: &Path) -> Result<Value> {
    let bytes = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
    serde_json::from_slice(&bytes).with_context(|| format!("parsing {}", path.display()))
}

fn reject_removed_contract_fields(value: &Value) -> Result<()> {
    let removed = [
        // `min_bq` used to reject fragments over a single disagreeing base.
        // Base quality now reaches the reader on the fragment bank quality planes,
        // and only the end-motif task thresholds it, so a configuration still
        // carrying the old key is rejected rather than silently ignored.
        "min_bq",
        "fragment_event_class",
        "event_class_cardinality",
        "fragment_start_offset_bp",
        "fragment_end_offset_bp",
        "region_width_bp",
        "coarse_offset_bins",
        "fine_offset_bins",
        "span_bins",
        "orientation_cardinality",
        "span_embedding",
        "orientation_embedding",
        "cut",
        "cut_mask_fraction",
        "structural_mask_fraction",
        "ref_context_bp",
        "ref_context_left_bp",
        "ref_context_right_bp",
        "cut_ref_inward_mask_min",
        "cut_ref_inward_mask_max",
        "cut_ref_min_visible_bp",
        "fixed_flank_bp",
        "left_flank_bp",
        "right_flank_bp",
    ];
    match value {
        Value::Object(values) => {
            for (key, child) in values {
                if removed.contains(&key.as_str()) {
                    bail!("removed contract field {key:?} is not accepted");
                }
                reject_removed_contract_fields(child)?;
            }
        }
        Value::Array(values) => {
            for child in values {
                reject_removed_contract_fields(child)?;
            }
        }
        Value::String(value) if removed.contains(&value.as_str()) => {
            bail!("removed contract field {value:?} is not accepted");
        }
        _ => {}
    }
    Ok(())
}

pub fn load_resolved(args: &ConfigArgs) -> Result<ResolvedConfig> {
    let raw = load_json(&args.config)?;
    reject_removed_contract_fields(&raw)?;
    let typed: RootConfig = serde_json::from_value(raw.clone())
        .context("deserializing fragment extraction configuration")?;
    validate(&typed)?;
    Ok(ResolvedConfig {
        min_mapq: typed.min_mapq,
        min_len: typed.min_len,
        max_len: typed.max_len,
        max_reference_mismatch_fraction: typed.max_reference_mismatch_fraction,
        max_mate_buffer: typed.max_mate_buffer.unwrap_or(1_000_000),
        chromatin_bin_size_bp: typed.chromatin_bin_size_bp.unwrap_or(10),
        raw,
        typed,
    })
}

fn validate(cfg: &RootConfig) -> Result<()> {
    let c = &cfg.coordinate_system;
    let r = &cfg.reservoir_sampling;

    let pair = &cfg.filters.pair_eligibility;
    if !pair.strict_trustworthy_ends
        || !pair.reject_clipping
        || !pair.reject_protruding_3p
        || !pair.require_full_fragment_coverage
        || !pair.require_complete_reference_window
    {
        bail!("all strict filters.pair_eligibility switches must be enabled");
    }
    if cfg.min_mapq != 30 || cfg.min_len != 30 || cfg.max_len != 299 {
        bail!("strict pair eligibility requires min_len=30, max_len=299, and min_mapq=30");
    }
    // The guard is a band rather than an equality because this is a structural
    // mismapping threshold with no single canonical value: it must stay far
    // above the per-fragment disagreement a real molecule shows (well under 1%)
    // and far below a mismapped read's.
    if !(0.01..=0.50).contains(&cfg.max_reference_mismatch_fraction) {
        bail!(
            "max_reference_mismatch_fraction must lie in [0.01, 0.50], got {}",
            cfg.max_reference_mismatch_fraction
        );
    }

    // An empty list would otherwise scan nothing and still write a valid,
    // empty sample directory, and a repeated chromosome would be scanned twice.
    if c.chromosomes.is_empty() {
        bail!("coordinate_system.chromosomes must name at least one chromosome");
    }
    let mut seen = std::collections::BTreeSet::new();
    for chrom in &c.chromosomes {
        if !seen.insert(chrom.as_str()) {
            bail!("coordinate_system.chromosomes lists {chrom:?} more than once");
        }
    }

    if c.region_window_bp == 0 || c.region_stride == 0 {
        bail!("region width and stride must be > 0");
    }
    if r.stored_fragments_per_region != DEFAULT_STORED_FRAGMENTS_PER_REGION {
        bail!(
            "reservoir_sampling.stored_fragments_per_region must be {} exactly, got {}",
            DEFAULT_STORED_FRAGMENTS_PER_REGION,
            r.stored_fragments_per_region
        );
    }

    // Fragment ATAC covariates are an optional add-on. With them disabled the
    // covariate plane is written as zeros and no chromatin atlas is needed, so
    // the declarations below are only checked when the add-on is switched on.
    let atac = &cfg.fragment_atac_covariates;
    if atac.enabled {
        if atac.mask_bin_size_bp == 0 {
            bail!("fragment_atac_covariates.mask_bin_size_bp must be > 0");
        }
        let names = atac
            .groups
            .iter()
            .map(|group| group.name.as_str())
            .collect::<Vec<_>>();
        if names.as_slice() != ["blood", "colon", "bladder"] {
            bail!(
                "fragment_atac_covariates.groups must be ordered as blood, colon, bladder; got {names:?}"
            );
        }
        if atac.groups.iter().any(|group| group.members.is_empty()) {
            bail!("every enabled fragment_atac_covariates group needs at least one atlas member");
        }
        let expected_members = [
            vec![
                "neutrophil__pc",
                "monocyte__pc",
                "b_cell__pc",
                "t_cell__pc",
                "natural_killer_cell__pc",
                "erythroblasts__pc",
                "mega__pc",
            ],
            vec!["colonic_mucosa__tis"],
            vec!["mucosa_of_urinary_bladder__tis"],
        ];
        for (group, expected) in atac.groups.iter().zip(expected_members) {
            if group.members.iter().map(String::as_str).collect::<Vec<_>>() != expected {
                bail!(
                    "fragment_atac_covariates group {:?} does not match the authoritative atlas members",
                    group.name
                );
            }
        }
        if atac.mask_bin_size_bp != cfg.chromatin_bin_size_bp.unwrap_or(10) {
            bail!("fragment_atac_covariates.mask_bin_size_bp must equal chromatin_bin_size_bp");
        }
        if atac.groups.len() != ATAC_TISSUE_GROUP_COUNT {
            bail!(
                "fragment_atac_covariates requires exactly {ATAC_TISSUE_GROUP_COUNT} tissue groups"
            );
        }
    }

    let seq = &cfg.feature_sets.fragment_sequence;
    if seq.seq_window_w != MOLECULAR_WINDOW_BP {
        bail!("feature_sets.fragment_sequence.seq_window_w must be {MOLECULAR_WINDOW_BP}");
    }
    let expected_arrays = [
        "fragment_ref",
        "fragment_r1",
        "fragment_r2",
        "fragment_r1_qual",
        "fragment_r2_qual",
        "fragment_r1_cycle",
        "fragment_r2_cycle",
    ];
    if seq
        .array_names
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>()
        .as_slice()
        != expected_arrays.as_slice()
    {
        bail!("fragment sequence arrays are missing or out of order");
    }

    Ok(())
}

pub fn run_validate(args: ConfigArgs) -> Result<()> {
    let cfg = load_resolved(&args)?;
    let out = serde_json::json!({
        "chromosomes": cfg.typed.coordinate_system.chromosomes,
        "region_window_bp": cfg.typed.coordinate_system.region_window_bp,
        "region_stride": cfg.typed.coordinate_system.region_stride,
        "sampling_mode": "natural",
        "stored_fragments_per_region": cfg.typed.reservoir_sampling.stored_fragments_per_region,
        "sequence_width": cfg.typed.feature_sets.fragment_sequence.seq_window_w,
        "filters": {
            "min_mapq": cfg.min_mapq,
            "max_reference_mismatch_fraction": cfg.max_reference_mismatch_fraction,
            "min_len": cfg.min_len,
            "max_len": cfg.max_len,
            "pair_eligibility": cfg.typed.filters.pair_eligibility,
            "patient_snp_blacklist_required": cfg.requires_snp_blacklist()
        }
    });
    println!("{}", serde_json::to_string_pretty(&out)?);
    Ok(())
}
