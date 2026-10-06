#!/usr/bin/env bash
#
# End-to-end smoke test for bam2frag on one real sample.
#
# Runs the cheapest thing that exercises the whole pipeline: validate the
# contract, extract a single small chromosome, then validate the output. A
# whole-genome run is the same command with the stock config.json; see the
# memory note in the README before submitting one.
#
# Usage:
#   scripts/smoke_test.sh --bam BAM --reference FASTA [options]
#
# Options:
#   --bam PATH          Coordinate-sorted, indexed BAM            (required)
#   --reference PATH    Indexed reference FASTA                   (required)
#   --sample-id ID      Biological sample ID                      (required)
#   --patient-id ID     Patient ID                                (default: sample ID)
#   --chrom NAME        Chromosome to restrict the test to        (default: chr21)
#   --output-dir PATH   Where to write the sample directory       (default: ./smoke_out/<sample>)
#   --threads N         Chromosome workers                        (default: 1)
#   --max-fragments N   Stop after N paired fragments; 0 = whole chromosome (default: 0)
#   --batch-size N      Chromosomes scanned before packing; the memory lever   (default: 1)
#
set -euo pipefail

BAM=""
REFERENCE=""
SAMPLE_ID=""
PATIENT_ID=""
CHROM="chr21"
OUTPUT_DIR=""
THREADS=1
MAX_FRAGMENTS=0
# One chromosome is all this test scans, so a batch of one keeps the smoke run
# inside a small allocation no matter which chromosome is picked.
BATCH_SIZE=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)           BAM="$2"; shift 2 ;;
        --reference)     REFERENCE="$2"; shift 2 ;;
        --sample-id)     SAMPLE_ID="$2"; shift 2 ;;
        --patient-id)    PATIENT_ID="$2"; shift 2 ;;
        --chrom)         CHROM="$2"; shift 2 ;;
        --output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --max-fragments) MAX_FRAGMENTS="$2"; shift 2 ;;
        --batch-size)    BATCH_SIZE="$2"; shift 2 ;;
        -h|--help)       sed -n '2,23p' "$0"; exit 0 ;;
        *) echo "[ERROR] Unknown option: $1" >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

fail() { echo "[ERROR] $*" >&2; exit 1; }
step() { echo; echo "=== $* ==="; }

[[ -n "$BAM" ]] || fail "--bam is required"
[[ -n "$REFERENCE" ]] || fail "--reference is required"
# Deliberately not derived from the BAM name: that usually carries an assay or
# material suffix such as _cfdna, which is not part of the sample's identity.
[[ -n "$SAMPLE_ID" ]] || fail "--sample-id is required"

step "Checking inputs"
[[ -s "$BAM" ]] || fail "BAM not found: $BAM"
# The extractor fetches by region, so both indexes must exist up front rather
# than failing an hour into a scan.
if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" && ! -s "${BAM}.csi" ]]; then
    fail "No BAM index beside $BAM (expected .bam.bai, .bai or .bam.csi)"
fi
[[ -s "$REFERENCE" ]] || fail "Reference not found: $REFERENCE"
[[ -s "${REFERENCE}.fai" ]] || fail "No FASTA index: ${REFERENCE}.fai (run: samtools faidx $REFERENCE)"
echo "bam        $BAM"
echo "reference  $REFERENCE"

PATIENT_ID="${PATIENT_ID:-$SAMPLE_ID}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO_ROOT/smoke_out/$SAMPLE_ID}"
echo "sample_id  $SAMPLE_ID"
echo "patient_id $PATIENT_ID"
echo "chrom      $CHROM"
echo "batch_size $BATCH_SIZE"
echo "output_dir $OUTPUT_DIR"

step "Locating the extractor"
if [[ -n "${BAM2FRAG_BIN:-}" ]]; then
    EXE="$BAM2FRAG_BIN"
elif command -v bam2frag-core >/dev/null 2>&1; then
    EXE="$(command -v bam2frag-core)"
elif [[ -x "$REPO_ROOT/rust/bam2frag/target/release/bam2frag-core" ]]; then
    EXE="$REPO_ROOT/rust/bam2frag/target/release/bam2frag-core"
else
    fail "bam2frag-core not found. Build it:
    cargo build --release --manifest-path $REPO_ROOT/rust/bam2frag/Cargo.toml"
fi
echo "$EXE"
"$EXE" version

step "Building a single-chromosome config"
# The stock config scans all 22 autosomes. Restricting it to one small
# chromosome is what makes this a smoke test rather than a production run.
SMOKE_CONFIG="$(mktemp -t bam2frag-smoke-XXXXXX.json)"
trap 'rm -f "$SMOKE_CONFIG"' EXIT
python3 - "$REPO_ROOT/config.json" "$SMOKE_CONFIG" "$CHROM" <<'PY'
import json, sys

source, destination, chrom = sys.argv[1], sys.argv[2], sys.argv[3]
with open(source) as handle:
    config = json.load(handle)
config["coordinate_system"]["chromosomes"] = [chrom]
with open(destination, "w") as handle:
    json.dump(config, handle, indent=2)
print(f"restricted to {chrom}")
PY

step "validate-config (stock)"
"$EXE" validate-config --config "$REPO_ROOT/config.json" | head -40

step "validate-config (smoke)"
"$EXE" validate-config --config "$SMOKE_CONFIG" >/dev/null
echo "ok"

step "extract"
rm -rf "$OUTPUT_DIR" "${OUTPUT_DIR}".incomplete.* 2>/dev/null || true
mkdir -p "$(dirname "$OUTPUT_DIR")"
export RUST_LOG="${RUST_LOG:-info}"
export RUST_BACKTRACE=1
/usr/bin/time -v "$EXE" extract \
    --bam "$BAM" \
    --reference "$REFERENCE" \
    --config "$SMOKE_CONFIG" \
    --sample-id "$SAMPLE_ID" \
    --patient-id "$PATIENT_ID" \
    --output-dir "$OUTPUT_DIR" \
    --threads "$THREADS" \
    --chromosome-batch-size "$BATCH_SIZE" \
    --seed 0 \
    --debug-max-fragments "$MAX_FRAGMENTS" \
    --overwrite 2>&1 | tail -40

step "inspect"
"$EXE" inspect --output-dir "$OUTPUT_DIR"

step "Output"
ls -lh "$OUTPUT_DIR"
echo
echo "--- extraction_stats.csv (first rows) ---"
head -3 "$OUTPUT_DIR/extraction_stats.csv"
echo
echo "--- sample_meta.json (selected keys) ---"
python3 - "$OUTPUT_DIR/sample_meta.json" <<'PY'
import json, sys

with open(sys.argv[1]) as handle:
    meta = json.load(handle)
for key in [
    "generator",
    "sample_id",
    "patient_id",
    "n_regions",
    "sampling_mode",
    "max_stored_fragments_per_region",
    "n_sampled_rows_kept",
    "fragment_atac_covariates_enabled",
    "seq_window_w",
    "config_sha256",
]:
    print(f"  {key}: {meta.get(key)}")
PY

step "Determinism check"
# Same seed, same inputs must give a byte-identical bank. This is the cheapest
# check that reservoir sampling and scan planning are not picking up ambient
# state such as thread scheduling order.
SECOND_DIR="${OUTPUT_DIR}.rerun"
rm -rf "$SECOND_DIR"
"$EXE" extract \
    --bam "$BAM" \
    --reference "$REFERENCE" \
    --config "$SMOKE_CONFIG" \
    --sample-id "$SAMPLE_ID" \
    --patient-id "$PATIENT_ID" \
    --output-dir "$SECOND_DIR" \
    --threads "$THREADS" \
    --chromosome-batch-size "$BATCH_SIZE" \
    --seed 0 \
    --debug-max-fragments "$MAX_FRAGMENTS" \
    --overwrite >/dev/null 2>&1
if cmp -s "$OUTPUT_DIR/fragment_bank.bin" "$SECOND_DIR/fragment_bank.bin" &&
   cmp -s "$OUTPUT_DIR/fragment_bank_index.bin" "$SECOND_DIR/fragment_bank_index.bin"; then
    echo "ok: identical bank and index across two runs at seed 0"
    rm -rf "$SECOND_DIR"
else
    fail "bank differs between two runs at the same seed: $OUTPUT_DIR vs $SECOND_DIR"
fi

step "PASSED"
echo "Sample directory: $OUTPUT_DIR"
