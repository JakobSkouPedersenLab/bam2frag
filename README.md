# bam2frag

A Rust tool to extract cfDNA fragments from BAM files very fast.

`bam2frag` scans a coordinate-sorted BAM, reconstructs every eligible paired-end cfDNA fragment against the reference, and stores a uniform random sample of them per genomic region in a compressed binary fragment bank. A thin Python CLI drives the Rust extractor on a per-sample basis.

The tool extracts **natural raw cfDNA fragments only**. It has no mutation input, no tumor-catalog or ctDNA-anchored extraction path, and no patient or sample covariate tables: the output describes the molecules in the BAM and nothing else.

For each sample the tool outputs:

- `fragment_bank.bin` — the fragments themselves, as zstd-compressed payloads
- `fragment_bank_index.bin` — the row index into the bank
- `region_available_count.npy` — how many eligible fragments each region actually had
- `sample_meta.json` — the resolved extraction contract and output manifest
- `extraction_stats.csv` / `extraction_stats.json` — per-chromosome acceptance and rejection counts

The pipeline is designed for whole-genome BAM files. Output is written to an `.incomplete.<pid>` directory and atomically renamed once every array and sidecar is in place, so a half-written sample directory never looks finished.

---

## Setting up your environment

### 1. Create a new Conda environment

```bash
conda create -n bam2frag python=3.12
conda activate bam2frag
```

### 2. Install Poetry

```bash
pip install poetry
poetry install
```

### 3. Build the Rust extractor

The Python CLI is a wrapper; the extraction itself is Rust and must be compiled once.

```bash
rustup toolchain install stable
cargo build --release --manifest-path rust/bam2frag/Cargo.toml
```

This produces `rust/bam2frag/target/release/bam2frag-core`, which the Python CLI finds automatically. The crate links `htslib`, so it builds on Linux and macOS but not on native Windows.

### 4. Install gwf, if needed

```bash
conda config --add channels gwforg
conda install gwf
```

---

## How to use

### Command line

```bash
bam2frag \
    /path/to/sample.bam \
    /path/to/reference.fasta \
    ./config.json \
    --output_dir /path/to/output/SAMPLE001 \
    --genomic_blacklist_file /path/to/genomic_blacklist.bed \
    --chromatin_dir /path/to/chromatin \
    --sample_id SAMPLE001 \
    --patient_id PATIENT001 \
    --threads 12 \
    --seed 0 \
    --overwrite \
    --verbose
```

The three positional arguments and `--sample_id` are required. Everything else, including the chromatin atlas, is optional, so the shortest useful run is:

```bash
bam2frag \
    /path/to/sample.bam \
    /path/to/reference.fasta \
    ./config.json \
    --sample_id SAMPLE001
```

Both the BAM and the reference FASTA must be indexed: the extractor fetches by region and reconstructs reference sequence through `faidx`.

To check a configuration without reading a BAM:

```bash
bam2frag sample.bam reference.fasta ./config.json --sample_id SAMPLE001 --validate_only
```

### Calling the Rust extractor directly

The Python CLI adds input checking, defaults and a pre-flight config validation. The underlying executable can also be called on its own:

```bash
rust/bam2frag/target/release/bam2frag-core extract \
    --bam sample.bam \
    --reference hg38.fa \
    --config config.json \
    --sample-id SAMPLE001 \
    --output-dir out/SAMPLE001 \
    --chromatin-dir /path/to/chromatin \
    --threads 12 --seed 0 --overwrite

rust/bam2frag/target/release/bam2frag-core validate-config --config config.json
rust/bam2frag/target/release/bam2frag-core inspect --output-dir out/SAMPLE001
```

---

## Arguments

### Required positional arguments

| Argument | Description |
|---|---|
| `bam_path` | Input BAM file. Must be coordinate-sorted and indexed. |
| `reference_path` | Reference FASTA (hg38). Must have a `.fai` index. |
| `config_path` | JSON extraction contract. See [Configuration](#configuration). |

### Output and logging

| Argument | Description |
|---|---|
| `--output_dir` | Destination directory for the extracted sample. Defaults to `./data/<sample_id>`. Parent directories are created automatically. |
| `--overwrite` | Replace an existing output directory. Without it, an existing directory is an error. |
| `--validate_only` | Validate the contract, print it, and exit without reading the BAM. |
| `--verbose` | Raise the Rust log level from `info` to `debug`. |

### Sample identity

| Argument | Description |
|---|---|
| `--sample_id` | **Required.** Biological sample ID, written to `sample_meta.json` and the statistics files. |
| `--patient_id` | Patient ID used to resolve patient-specific blacklists. Defaults to the sample ID, treating one sample as one patient. |

`--sample_id` is required and is never inferred from the BAM path. A BAM file name usually carries an assay or material suffix that is **not** part of the sample's identity — in

```text
.../C05638/C72A05638D_cfdna_N295-284/C72A05638D_cfdna.aligned.sorted.markdup.bam
```

the sample is `C72A05638D`; `_cfdna` says which material was sequenced. Guessing would silently write `C72A05638D_cfdna` into every row of output, and no downstream analysis would notice.

Naming schemes are cohort-specific, so they belong in the workflow rather than the tool. `src/gwf/workflow.py` has `get_sample_id` and `get_patient_id` for exactly this, implemented for the layout above and meant to be replaced:

```python
MATERIAL_SUFFIXES = ("_cfdna", "_buffycoat", "_plasma")
# .../<patient_id>/<sample_id>_cfdna_<library>/<sample_id>_cfdna.*.bam
#   -> sample_id=C72A05638D, patient_id=C05638
```

### ATAC covariates (optional add-on)

| Argument | Description |
|---|---|
| `--chromatin_dir` | Root of the aggregate tissue ATAC atlas. |

Fragment ATAC covariates are an **optional add-on, disabled by default**. With them off, bam2frag needs no chromatin atlas, the `--chromatin_dir` argument is unused, and the covariate plane is written as zeros for every fragment. The plane is always present, so the bank layout does not change when the add-on is switched on or off.

To enable it, set `fragment_atac_covariates.enabled` to `true` in `config.json` and pass `--chromatin_dir`. The extractor then expects one sparse `.npz` coverage mask per atlas member and chromosome:

```text
<chromatin_dir>/chromatin_masks/ATAC/sparse_overlaps_by_chromosome/<member>/<chrom>.npz
```

Each `.npz` holds `data.npy`, `indices.npy`, `indptr.npy` and `shape.npy` for a CSR or CSC matrix whose values are per-bin covered fractions in `[0, 1]`, binned at `mask_bin_size_bp`. The `<member>` names are the ones listed under `fragment_atac_covariates.groups` in the config, and the group names and ordering are validated against the atlas the shipped config describes.

`sample_meta.json` records `fragment_atac_covariates_enabled`, so a reader can always tell whether the plane holds real values or zeros.

### Blacklist filters

| Argument | Description |
|---|---|
| `--genomic_blacklist_file` | BED-like interval file with genomic regions to exclude. |
| `--donor_blacklist_file` | Healthy-donor positional blacklist, BED-like intervals. |
| `--snp_blacklist_path` | Directory of patient-specific variant blacklists. |
| `--lowpower_blacklist_path` | Patient-specific low-power region blacklist. |
| `--problematic_blacklist_path` | Patient-specific problematic-region blacklist. |

All are optional. The last four are only consulted when the matching switch under `filters` in the config is enabled; they are disabled in the shipped contract. The SNP blacklist argument is a **directory**, laid out per patient:

```text
<snp_blacklist_path>/<patient_id>/<patient_id>_snvs_calls.tsv
<snp_blacklist_path>/<patient_id>/<patient_id>_indels_calls.tsv
```

### Runtime

| Argument | Description |
|---|---|
| `--threads` | Number of chromosome workers. Defaults to `$SLURM_CPUS_PER_TASK`, or 1. |
| `--chromosome_batch_size` | How many chromosomes are scanned before their fragments are packed. Default 12. The dominant memory control, see [Memory](#memory). |
| `--seed` | Seed for reservoir sampling. Extraction is deterministic for a given seed, BAM and contract. |
| `--progress_every_fragments` | Emit a progress line every N processed paired fragments. `0` disables. Default 5,000,000. |
| `--debug_max_fragments` | Stop after N paired fragments per chromosome. `0` disables. For smoke tests. |

---

## Output

For an output directory such as:

```bash
--output_dir /path/to/output/SAMPLE001
```

the pipeline writes:

```text
/path/to/output/SAMPLE001/
  fragment_bank.bin
  fragment_bank_index.bin
  region_available_count.npy
  sample_meta.json
  extraction_stats.csv
  extraction_stats.json
```

| File | Description |
|---|---|
| `fragment_bank.bin` | One zstd frame per region, each holding that region's fragment payload. |
| `fragment_bank_index.bin` | The bank index. One row per region, giving its offset and stored frame size. |
| `region_available_count.npy` | `uint64` per region: how many eligible fragments the region had before reservoir capping. Row `i` corresponds exactly to bank index row `i`. |
| `sample_meta.json` | Resolved contract, array manifest, `config_sha256`, and the counts of fragments seen, kept and dropped to overflow. |
| `extraction_stats.csv` | Per-chromosome acceptance and rejection counters. |
| `extraction_stats.json` | The same statistics in JSON. |

### How fragments are selected

The genome is divided into fixed regions (1 Mb by default). Within each region, every paired-end fragment passing the eligibility filters is a candidate, and the extractor keeps a **uniform random sample** of up to `stored_fragments_per_region` (10,000) of them via reservoir sampling. Regions with fewer eligible fragments keep all of them; `region_available_count.npy` records the true count, so downstream work can tell a sparse region from a capped one.

This is a natural sample: no mutation, locus or allele influences which fragments are kept.

### Compression

Each region's payload is stored as a single zstd frame at level 3. The bytes inside the frame are a byte-identical uncompressed fragment payload, and the index is unchanged apart from each row's `nbytes` recording the stored frame size — compression is how the bank is written, not a different representation.

Most of a payload is zero padding outside the reads plus monotone cycle ramps, so a whole-sample bank compresses roughly **12x** (44 GB down to about 3.6 GB on one colon sample). Level 3 already reaches that ratio; higher levels cost write time for little further gain.

---

## Fragment representation

Every stored fragment is a fixed-width 320 bp window, centred on the fragment:

```text
left_cut_idx  = floor((320 - fragment_span_bp) / 2)
right_cut_idx = left_cut_idx + fragment_span_bp - 1
```

A fragment is only stored when the full 320 bp reference window is available, so no window is truncated at a contig edge.

Each fragment is written across nine planes:

| Plane | Description |
|---|---|
| `fragment_ref` | Reference base at each window position. |
| `fragment_r1` | Read 1 base at each window position. |
| `fragment_r2` | Read 2 base at each window position. |
| `fragment_r1_qual` | Read 1 base quality, 6-bit, scaled to a maximum of 63. |
| `fragment_r2_qual` | Read 2 base quality, 6-bit. |
| `fragment_r1_cycle` | Read 1 sequencing cycle: `0` = absent, otherwise the original cycle plus one. |
| `fragment_r2_cycle` | Read 2 sequencing cycle, same encoding. |
| `fragment_atac_covariates` | Three `float32` base-weighted ATAC overlap fractions: blood, colon, bladder. All zero unless the optional ATAC add-on is enabled. |
| `fragment_structural` | `fragment_span_bp` (`u16`), `fragment_orientation` (`u8`), one reserved zero byte. |

Bases are 4-bit coded `N=0, A=1, C=2, G=3, T=4`. Qualities are 6-bit. Fragment endpoints are the **validated outer 5′** positions of the two mates.

Note that base quality is carried into the bank on the two quality planes rather than being thresholded away at extraction: a low-quality base is stored with its quality, not discarded.

### Fragment eligibility

A read pair must clear every one of the following to be stored. `extraction_stats.csv` reports the counter for each rejection reason, which makes it straightforward to see why a sample yielded fewer fragments than expected:

| Counter | Meaning |
|---|---|
| `pairs_accepted` | Pairs stored as candidate fragments. |
| `pairs_rejected_non_inward_orientation` | Not an inward-facing pair. |
| `pairs_rejected_ambiguous_5p` | The 5′ end could not be resolved unambiguously. |
| `pairs_rejected_5p_clipping` | Soft or hard clipping at a 5′ end. |
| `pairs_rejected_protruding_3p` | A 3′ end protrudes past its mate's 5′ end. |
| `pairs_rejected_unpaired_or_mate_unmapped` | Unpaired, or the mate is unmapped. |
| `pairs_rejected_improper_pair` | Not flagged as a proper pair. |
| `pairs_rejected_invalid_flags` | Secondary, supplementary, duplicate or QC-fail. |
| `pairs_rejected_different_chromosomes` | Mates on different chromosomes. |
| `pairs_rejected_mapq` | Below `min_mapq`. |
| `pairs_rejected_fragment_span` | Span outside `[min_len, max_len]`. |
| `pairs_rejected_incomplete_fragment_coverage` | The pair does not cover the full fragment span. |
| `pairs_rejected_reference_mismatch_burden` | Disagreement with the reference above `max_reference_mismatch_fraction`, a structural mismapping guard. |
| `pairs_rejected_indel_cigar` | Indel in the CIGAR. |
| `pairs_rejected_invalid_5p_span` | The 5′ span is not valid. |
| `pairs_rejected_incomplete_reference_window` | The full 320 bp reference window was unavailable. |

The statistics also report `accepted_covered_positions`, `accepted_unilateral_disagreements`, `accepted_discordant_disagreements`, `accepted_concordant_mismatches` and `accepted_without_end_motif_bases` for the accepted pairs.

---

## Configuration

Extraction settings are defined in `config.json`. The contract is strict and self-validating: `validate-config` rejects unknown model-facing features and removed fields rather than silently ignoring them.

### Key fields

| Field | Description |
|---|---|
| `coordinate_system.genome` | Reference genome. `hg38`. |
| `coordinate_system.chromosomes` | Chromosomes to process. Either `"all"` or an explicit list. See [Selecting chromosomes](#selecting-chromosomes). |
| `coordinate_system.region_window_bp` | Region width in bp. Regions tile each chromosome contiguously. Default 1,000,000. |
| `coordinate_system.region_stride` | Width of each chromosome's slice of the **region-id space**, not a genomic distance. Region ids are `chromosome_index * region_stride + region_bin`. Default 1,000,000. |
| `reservoir_sampling.stored_fragments_per_region` | Fragments kept per region. Default 10,000. |
| `reservoir_sampling.fetch_padding_bp` | Padding added when fetching a region, so fragments overlapping the edge are seen. |
| `reservoir_sampling.seed` | Default seed, overridden by `--seed`. |
| `min_len` / `max_len` | Fragment span bounds in bp. Default 30 and 299. |
| `min_mapq` | Minimum mapping quality. Default 30. |
| `max_reference_mismatch_fraction` | Largest fraction of the validated span that may disagree with the reference. A mismapping guard, not a quality threshold. Default 0.10. |
| `max_mate_buffer` | Maximum number of reads buffered while waiting for a mate. |
| `filters.pair_eligibility.*` | The five strict eligibility switches. All must be enabled. |
| `filters.*.enabled` | Whether each blacklist is consulted. All disabled by default. |
| `feature_sets.fragment_sequence` | `seq_window_w` (the 320 bp window) and `array_names` (the seven sequence planes, in order). |
| `fragment_atac_covariates` | Optional ATAC add-on: `enabled`, plus bin size, tissue groups and their atlas members. Off by default. |

The strict profile is coupled: enabling all `pair_eligibility` switches requires `min_len=30`, `max_len=299` and `min_mapq=30`, and `validate-config` enforces that.

### What is not configurable

Two things that look like settings are not.

**The fragment encoding** has exactly one form, so it is a property of the output: 4-bit bases `N=0, A=1, C=2, G=3, T=4`, 6-bit qualities scaled to 63, the centred-window formula, the `0 = absent` cycle-plane encoding, and the ATAC plane's `float32` fractions aggregated by element-wise max over atlas members. These are fixed in `arrays.rs` and `atac.rs` and documented under [Fragment representation](#fragment-representation).

**Scanning is always whole-region.** Regions tile each chromosome contiguously and every region is read in full, so there is no scan plan to configure. Earlier revisions carried a block/window subsampler — `blocks_per_region`, `block_bp`, `windows_per_block`, `window_bp`, `sampling`, `merge_overlapping_windows`, `collapse_blocks_to_region_token` — that every configuration ran in whole-region mode, where it built one window per region and merged them straight back into one span per chromosome. It is gone, along with `reservoir_sampling.enabled`: the reservoir is what the tool does.

A config that still declares any of these — whether a `dna_code` block or a `scan` block carried over from an older contract — is **rejected**, not ignored, so an extraction can never quietly run under rules its config appears to set.

### Selecting chromosomes

`coordinate_system.chromosomes` takes either the shortcut `"all"`:

```json
"coordinate_system": {
  "genome": "hg38",
  "chromosomes": "all",
  "region_window_bp": 1000000,
  "region_stride": 1000000
}
```

or an explicit list, when only part of the genome is wanted:

```json
"chromosomes": ["chr1", "chr2", "chr3"]
```

`"all"` expands to the 22 autosomes, `chr1` through `chr22`, while the config is loaded. The match is case-insensitive, and `validate-config` prints the expanded list.

**The sex chromosomes are not included in `"all"`.** Their copy number depends on donor sex, so including them implicitly would change what a region means between samples. Name them explicitly to extract them:

```json
"chromosomes": ["chr1", "chr2", "chrX", "chrY"]
```

An empty list, a repeated chromosome, or an unrecognised shortcut is rejected at load time rather than producing an empty or double-scanned run. Every configured chromosome must also be present in the BAM header, which the extractor checks before scanning starts.

---

## Workflow

`src/gwf/` holds a worked example for running one extraction target per BAM on a cluster with [gwf](https://gwf.app):

- `src/gwf/template.py` — the per-sample target: stages the reference to node-local scratch, clears any partial directory, and runs `bam2frag`
- `src/gwf/workflow.py` — reads a list of BAM paths and builds one target per sample

Both are templates to copy and adapt. The site-specific settings at the top of `workflow.py` are placeholders, so no real paths or cluster account names are in version control. Fill in `slurm_account`, `bam_file_list_path`, `reference_file` and `genomic_blacklist_file` — and `chromatin_dir` only if you enable the ATAC add-on — then:

```bash
gwf -f src/gwf/workflow.py status
gwf -f src/gwf/workflow.py run
```

The template declares `sample_meta.json` among its outputs. Because that file is written only after the atomic rename, gwf treats a sample as done only once extraction fully finished, and reruns anything that was interrupted.

---

## Parallel execution

The extractor runs one worker per chromosome through Rayon. Set the worker count with `--threads`, or leave it unset and let it follow the environment variable:

```bash
SLURM_CPUS_PER_TASK=12 bam2frag \
    /path/to/sample.bam \
    /path/to/reference.fasta \
    ./config.json \
    --sample_id SAMPLE001
```

If `SLURM_CPUS_PER_TASK` is not set and `--threads` is not given, the pipeline defaults to one worker.

Each worker opens its own BAM reader and reference cache, and every thread reads the reference through `faidx` for the whole run, which is why the gwf template copies the FASTA to node-local scratch before starting.

### Memory

Peak memory is driven by the scan phase, and the control for it is `--chromosome_batch_size`, not `--threads`.

Chromosomes are scanned in batches. A batch's surviving fragment candidates — each holding both mates' BAM records — all stay resident until the whole batch has been scanned, and only then are they packed and released. So peak memory scales with the **total region count of the largest batch**, roughly:

```text
peak_bytes  ≈  regions_in_largest_batch × stored_fragments_per_region × ~1 KB
```

At the default batch size of 12, the first batch is `chr1`–`chr12`: about 2,100 regions at up to 10,000 fragments each, so a whole-genome run peaks near **20 GB**.

Lowering `--threads` does **not** lower this peak — it only scans fewer chromosomes concurrently *within* the same batch. Lowering `--chromosome_batch_size` does, at the cost of packing in smaller bursts:

| `--chromosome_batch_size` | Largest batch | Approximate peak |
|---|---|---|
| 12 (default) | `chr1`–`chr12`, ~2,100 regions | ~20 GB |
| 6 | `chr1`–`chr6`, ~1,230 regions | ~12 GB |
| 4 | `chr1`–`chr4`, ~880 regions | ~9 GB |
| 2 | `chr1`–`chr2`, ~490 regions | ~5 GB |
| 1 | `chr1`, ~249 regions | ~3 GB |

Add headroom for the mate buffers and per-worker reference caches. Budget **50 GB** for a default whole-genome run, or lower the batch size to fit a smaller allocation. Raising it above `--threads` costs memory without buying scan parallelism, and the extractor logs a note when you do.

The other lever is `reservoir_sampling.stored_fragments_per_region`, which scales peak memory and output size together.

---

## Testing on a cluster

`scripts/smoke_test.sh` runs the cheapest thing that exercises the whole pipeline on a real BAM: it validates the contract, extracts one small chromosome, validates the output with `inspect`, and re-runs at the same seed to confirm the bank is byte-identical.

```bash
srun --account <account> --mem 25g --cpus-per-task 4 --time 2:00:00 --pty bash

cargo build --release --manifest-path rust/bam2frag/Cargo.toml
cargo test --manifest-path rust/bam2frag/Cargo.toml

scripts/smoke_test.sh \
    --bam /path/to/sample.bam \
    --reference /path/to/hg38.fna \
    --sample-id SAMPLE001 \
    --patient-id PATIENT001 \
    --chrom chr21 \
    --threads 4
```

`chr21` is about 47 Mb, so the smoke test covers roughly 47 regions and produces a bank of a few hundred MB — comfortable in 25 GB. Add `--max-fragments 2000000` to cut it to a couple of minutes when you only want to see the pipeline run end to end.

Once that passes, the same sample genome-wide is the stock `config.json` with no `--chrom` restriction, in a larger allocation:

```bash
sbatch --account <account> --mem 60g --cpus-per-task 12 --time 24:00:00 \
    --wrap "bam2frag /path/to/sample.bam /path/to/hg38.fna ./config.json \
            --sample_id SAMPLE001 --patient_id PATIENT001 \
            --output_dir /path/to/out/SAMPLE001 \
            --threads 12 --chromosome_batch_size 12 --overwrite --verbose"
```

In a tighter allocation, trade packing bursts for memory by lowering the batch size:

```bash
sbatch --account <account> --mem 25g --cpus-per-task 4 --time 24:00:00 \
    --wrap "bam2frag /path/to/sample.bam /path/to/hg38.fna ./config.json \
            --sample_id SAMPLE001 --patient_id PATIENT001 \
            --output_dir /path/to/out/SAMPLE001 \
            --threads 4 --chromosome_batch_size 4 --overwrite --verbose"
```

---

## Tests

Python:

```bash
poetry install --with dev
poetry run pytest
```

Rust:

```bash
cargo test --manifest-path rust/bam2frag/Cargo.toml
cargo fmt --all --manifest-path rust/bam2frag/Cargo.toml -- --check
```

---

## Contributions

Contributions are welcome. If you find a bug or have an idea for a feature, please open an issue or submit a pull request.

---

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
