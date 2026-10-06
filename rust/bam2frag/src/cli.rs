use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

/// Extract representative natural cfDNA fragments from BAM files into the
/// compressed fragment-bank dataset format.
#[derive(Debug, Parser)]
#[command(
    name = "bam2frag",
    version,
    about = "cfDNA fragment-reservoir BAM extractor",
    propagate_version = true
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Extract one BAM file into a model-ready sample directory.
    Extract(ExtractArgs),

    /// Validate the extraction configuration, then print its resolved
    /// extraction contract.
    ValidateConfig(ConfigArgs),

    /// Validate and print metadata for an existing extracted sample.
    Inspect(InspectArgs),

    /// Print the package version.
    Version,
}

/// Canonical fragment extraction configuration.
///
/// All extraction, feature, storage and filtering settings are defined in one
/// config.json file.
#[derive(Debug, Clone, Args)]
pub struct ConfigArgs {
    /// Main extraction/data configuration JSON.
    #[arg(long = "config", value_name = "FILE")]
    pub config: PathBuf,
}

#[derive(Debug, Clone, Args)]
pub struct ExtractArgs {
    /// Coordinate-sorted, indexed input BAM file.
    #[arg(long, value_name = "FILE")]
    pub bam: PathBuf,

    /// Indexed reference FASTA used for sequence reconstruction and event
    /// normalization.
    #[arg(long, value_name = "FILE")]
    pub reference: PathBuf,

    /// Biological sample identifier written to the output metadata.
    #[arg(long = "sample-id", value_name = "ID")]
    pub sample_id: String,

    /// Patient identifier used to resolve patient-specific filters. When
    /// omitted, --sample-id is used.
    #[arg(long = "patient-id", value_name = "ID")]
    pub patient_id: Option<String>,

    /// Destination directory for arrays, summaries, statistics, and metadata.
    #[arg(long = "output-dir", value_name = "DIR")]
    pub output_dir: PathBuf,

    #[command(flatten)]
    pub configs: ConfigArgs,

    /// Static chromatin atlas root, required when fragment ATAC covariates are enabled.
    #[arg(long = "chromatin-dir", value_name = "DIR")]
    pub chromatin_dir: Option<PathBuf>,

    /// Genome-wide positional blacklist in BED-like interval format.
    #[arg(long = "genomic-blacklist", value_name = "FILE")]
    pub genomic_blacklist: Option<PathBuf>,

    /// Healthy-donor positional blacklist in BED-like interval format.
    #[arg(long = "donor-blacklist", value_name = "FILE")]
    pub donor_blacklist: Option<PathBuf>,

    /// Base directory containing exact patient-specific variant blacklists:
    /// {patient_id}/{patient_id}_snvs_calls.tsv and
    /// {patient_id}/{patient_id}_indels_calls.tsv.
    #[arg(long = "snp-blacklist", value_name = "DIR")]
    pub snp_blacklist: Option<PathBuf>,

    /// Patient-specific low-power positional blacklist.
    #[arg(long = "lowpower-blacklist", value_name = "FILE")]
    pub lowpower_blacklist: Option<PathBuf>,

    /// Patient-specific problematic-region positional blacklist.
    #[arg(long = "problematic-blacklist", value_name = "FILE")]
    pub problematic_blacklist: Option<PathBuf>,

    /// Number of chromosome workers used by Rayon.
    #[arg(long, default_value_t = 1, value_name = "N")]
    pub threads: usize,

    /// How many chromosomes are scanned before their fragments are packed.
    ///
    /// This is the dominant memory control, and it is not `--threads`. A
    /// batch's surviving candidates all hold both mates' BAM records and stay
    /// resident until the whole batch is scanned, so peak memory scales with
    /// the total region count of the largest batch. Lowering `--threads` does
    /// not change that; lowering this does, at the cost of packing in smaller
    /// bursts. Raising it above `--threads` buys nothing but memory use.
    #[arg(long = "chromosome-batch-size", default_value_t = 12, value_name = "N")]
    pub chromosome_batch_size: usize,

    /// Random seed used for reservoir sampling.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,

    /// Emit progress after this many processed paired fragments. Set to zero
    /// to disable periodic progress messages.
    #[arg(
        long = "progress-every-fragments",
        default_value_t = 5_000_000,
        value_name = "N"
    )]
    pub progress_every_fragments: u64,

    /// Stop after this many processed paired fragments per chromosome. Zero
    /// disables the debug limit.
    #[arg(long = "debug-max-fragments", default_value_t = 0, value_name = "N")]
    pub debug_max_fragments: u64,

    /// Replace an existing output directory.
    #[arg(long)]
    pub overwrite: bool,
}

#[derive(Debug, Clone, Args)]
pub struct InspectArgs {
    /// Extracted sample directory to validate and inspect.
    #[arg(long = "output-dir", value_name = "DIR")]
    pub output_dir: PathBuf,
}
