#!/usr/bin/env python3
"""Per-sample command line entry point for the bam2frag cfDNA fragment extractor.

This is a thin orchestration layer. All extraction work is done by the Rust
executable in ``rust/bam2frag``; this script resolves that executable, validates
the extraction contract, maps per-sample inputs onto its ``extract`` subcommand,
and reports the resulting sample directory.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Name of the compiled Rust executable. It differs from this script's console
# entry point ("bam2frag") so that a PATH lookup can never resolve back here.
RUST_EXECUTABLE = "bam2frag-core"

ENV_OVERRIDE = "BAM2FRAG_BIN"

REPO_ROOT = Path(__file__).resolve().parents[1]

# Written by the Rust extractor only after every array and sidecar is in place,
# so its presence marks a complete sample directory.
COMPLETION_MARKER = "sample_meta.json"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="bam2frag",
        description="Extract natural cfDNA fragments from a single BAM file into a "
                    "compressed fragment bank.",
    )

    parser.add_argument("bam_path", type=str, help="Coordinate-sorted, indexed BAM file")
    parser.add_argument("reference_path", type=str, help="Indexed reference FASTA (hg38)")
    parser.add_argument("config_path", type=str, help="Extraction contract JSON")

    parser.add_argument("--output_dir", type=str, default=None,
                        help="Destination directory for the extracted sample. Defaults to "
                             "./data/<sample_id>.")
    parser.add_argument("--chromatin_dir", type=str, default=None,
                        help="Chromatin atlas root. Only needed when the config enables "
                             "fragment_atac_covariates, which is off by default.")

    parser.add_argument("--sample_id", type=str, required=True,
                        help="Biological sample ID written to the output metadata. Required: "
                             "it cannot be derived from the BAM file name, which commonly "
                             "carries an assay or material suffix that is not part of the "
                             "sample identity.")
    parser.add_argument("--patient_id", type=str, default=None,
                        help="Patient ID used to resolve patient-specific blacklists. "
                             "Defaults to the sample ID, treating one sample as one patient.")

    parser.add_argument("--genomic_blacklist_file", type=str, default=None)
    parser.add_argument("--donor_blacklist_file", type=str, default=None)
    parser.add_argument("--snp_blacklist_path", type=str, default=None)
    parser.add_argument("--lowpower_blacklist_path", type=str, default=None)
    parser.add_argument("--problematic_blacklist_path", type=str, default=None)

    parser.add_argument("--threads", type=int, default=None,
                        help="Chromosome workers. Defaults to $SLURM_CPUS_PER_TASK, or 1.")
    parser.add_argument("--chromosome_batch_size", type=int, default=12,
                        help="How many chromosomes are scanned before their fragments are "
                             "packed. This, not --threads, is the dominant memory control: "
                             "lower it to fit a smaller allocation.")
    parser.add_argument("--seed", type=int, default=0,
                        help="Seed for reservoir sampling.")
    parser.add_argument("--progress_every_fragments", type=int, default=5_000_000,
                        help="Emit progress every N paired fragments. 0 disables progress.")
    parser.add_argument("--debug_max_fragments", type=int, default=0,
                        help="Stop after N paired fragments per chromosome. 0 disables.")

    parser.add_argument("--overwrite", action="store_true",
                        help="Replace an existing output directory.")
    parser.add_argument("--validate_only", action="store_true",
                        help="Validate the extraction contract, print it, and exit without "
                             "reading the BAM.")
    parser.add_argument("--verbose", action="store_true",
                        help="Raise the Rust log level from info to debug.")

    return parser.parse_args(argv)


def resolve_executable() -> str:
    """
    Locate the compiled Rust extractor.

    The search order is explicit override, then PATH, then the in-tree cargo
    build directories, so a cluster deployment can pin one executable while a
    working copy falls back to its own `cargo build`.
    """
    suffix = ".exe" if os.name == "nt" else ""

    override = os.environ.get(ENV_OVERRIDE)
    if override:
        candidate = Path(override).expanduser()
        if not candidate.is_file():
            raise FileNotFoundError(
                f"${ENV_OVERRIDE} is set to {candidate}, which is not a file"
            )
        return str(candidate)

    on_path = shutil.which(RUST_EXECUTABLE)
    if on_path:
        return on_path

    crate_dir = REPO_ROOT / "rust" / "bam2frag" / "target"
    candidates = [
        crate_dir / "release" / f"{RUST_EXECUTABLE}{suffix}",
        crate_dir / "debug" / f"{RUST_EXECUTABLE}{suffix}",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)

    raise FileNotFoundError(
        f"Could not find the {RUST_EXECUTABLE} executable. Build it with "
        f"`cargo build --release --manifest-path rust/bam2frag/Cargo.toml`, put it on PATH, "
        f"or set ${ENV_OVERRIDE}. Checked PATH and: "
        + ", ".join(str(path) for path in candidates)
    )


def require_file(path: str, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} does not exist or is not a file: {resolved}")
    return resolved


def require_dir(path: str, label: str) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_dir():
        raise NotADirectoryError(f"{label} does not exist or is not a directory: {resolved}")
    return resolved


def require_bam_index(bam: Path) -> Path:
    """
    Return the BAM index. The extractor fetches by region, so an index is required.
    """
    candidates = [Path(f"{bam}.bai"), bam.with_suffix(".bai"), Path(f"{bam}.csi")]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"No BAM index found for {bam}. Looked for: "
        + ", ".join(str(path) for path in candidates)
    )


def default_threads() -> int:
    return int(os.environ.get("SLURM_CPUS_PER_TASK", 1))


def run(command: list[str], env: dict[str, str], label: str) -> None:
    print(f"[bam2frag] {label}: {' '.join(command)}", flush=True)
    result = subprocess.run(command, env=env)
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    executable = resolve_executable()
    config = require_file(args.config_path, "config")

    env = dict(os.environ)
    env["RUST_LOG"] = "debug" if args.verbose else "info"
    env.setdefault("RUST_BACKTRACE", "1")

    # The contract is validated on its own first, so a malformed config fails in
    # seconds rather than after the BAM scan has started.
    validate_command = [executable, "validate-config", "--config", str(config)]
    if args.validate_only:
        run(validate_command, env, "validate-config")
        return 0
    validated = subprocess.run(validate_command, env=env, capture_output=True, text=True)
    if validated.returncode != 0:
        # The resolved contract is only interesting when it fails to resolve.
        sys.stderr.write(validated.stdout)
        sys.stderr.write(validated.stderr)
        raise SystemExit(validated.returncode)

    bam = require_file(args.bam_path, "BAM")
    require_bam_index(bam)
    reference = require_file(args.reference_path, "reference FASTA")
    require_file(f"{reference}.fai", "reference FASTA index")

    sample_id = args.sample_id
    # One sample is treated as one patient unless told otherwise. That is only a
    # labelling default; it decides which patient-specific blacklists are read.
    patient_id = args.patient_id or sample_id

    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir \
        else (REPO_ROOT / "data" / sample_id)
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    threads = args.threads if args.threads is not None else default_threads()

    command = [
        executable, "extract",
        "--bam", str(bam),
        "--reference", str(reference),
        "--config", str(config),
        "--sample-id", sample_id,
        "--patient-id", patient_id,
        "--output-dir", str(output_dir),
        "--threads", str(threads),
        "--chromosome-batch-size", str(args.chromosome_batch_size),
        "--seed", str(args.seed),
        "--progress-every-fragments", str(args.progress_every_fragments),
        "--debug-max-fragments", str(args.debug_max_fragments),
    ]

    optional_paths = {
        "--chromatin-dir": (args.chromatin_dir, require_dir, "chromatin atlas"),
        "--genomic-blacklist": (args.genomic_blacklist_file, require_file, "genomic blacklist"),
        "--donor-blacklist": (args.donor_blacklist_file, require_file, "donor blacklist"),
        "--snp-blacklist": (args.snp_blacklist_path, require_dir, "SNP blacklist root"),
        "--lowpower-blacklist": (args.lowpower_blacklist_path, require_file, "low-power blacklist"),
        "--problematic-blacklist": (args.problematic_blacklist_path, require_file,
                                    "problematic-region blacklist"),
    }
    for flag, (value, check, label) in optional_paths.items():
        if value:
            command += [flag, str(check(value, label))]

    if args.overwrite:
        command.append("--overwrite")

    print(f"[bam2frag] sample_id={sample_id} patient_id={patient_id} "
          f"threads={threads} output_dir={output_dir}", flush=True)
    run(command, env, "extract")

    marker = output_dir / COMPLETION_MARKER
    if not marker.is_file():
        raise SystemExit(f"[bam2frag] extraction reported success but {marker} is missing")

    print(f"[bam2frag] wrote {output_dir}", flush=True)
    return 0


def cli() -> int:
    """
    Console entry point. Expected input problems are reported as one line rather
    than a traceback; anything unexpected still raises so it can be debugged.
    """
    try:
        return main()
    except (FileNotFoundError, NotADirectoryError, ValueError) as error:
        sys.stderr.write(f"[bam2frag] error: {error}\n")
        return 2


if __name__ == "__main__":
    sys.exit(cli())
