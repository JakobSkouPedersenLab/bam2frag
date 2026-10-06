#!/usr/bin/env python3
"""Example gwf workflow: one bam2frag extraction target per BAM file.

This is a template to copy and adapt, not a production workflow. Fill in the
site-specific settings below; they are placeholders so that no real paths or
cluster account names end up in version control.

Run it with:

    gwf -f src/gwf/workflow.py status
    gwf -f src/gwf/workflow.py run
"""
from __future__ import annotations

from pathlib import Path
from typing import List

from gwf import Workflow

from src.gwf.template import make_fragments


def read_text_file(path: str) -> List[str]:
    """
    Reads a text file
    """
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def resolve_patient_dir(base: Path, patient_id: str) -> Path | None:
    """
    Returns the actual on-disk patient directory under `base`,
    matching case-insensitively.
    """
    if not base.exists():
        return None

    pid_lower = patient_id.lower()

    for child in base.iterdir():
        if child.is_dir() and child.name.lower() == pid_lower:
            return child

    return None


# The two functions below are where a cohort's naming scheme lives. bam2frag
# itself never guesses either ID: a BAM file name usually carries an assay or
# material suffix that is not part of the sample's identity, and a tool that
# stripped it by guessing would mislabel output silently.
#
# The implementations here are a worked example for paths shaped like:
#
#   .../main_data/<patient_id>/<sample_id>_cfdna_<library>/<sample_id>_cfdna.aligned.sorted.markdup.bam
#
# which yields sample_id=C72A05638D, patient_id=C05638 for
#
#   .../main_data/C05638/C72A05638D_cfdna_N295-284/C72A05638D_cfdna.aligned.sorted.markdup.bam
#
# Replace both for any other layout.

# Material and assay suffixes that appear in file names but are not part of the
# sample ID. Extend for other assay types rather than widening the split.
MATERIAL_SUFFIXES = ("_cfdna", "_buffycoat", "_plasma")


def get_sample_id(bam_path: str) -> str:
    """
    Returns the biological sample ID for a BAM path.

    The file name's leading token is the sample ID plus a material suffix, so
    the suffix is stripped explicitly rather than by splitting on the first
    underscore, which would also truncate IDs that legitimately contain one.
    """
    leading_token = Path(bam_path).name.split(".")[0]
    lowered = leading_token.lower()
    for suffix in MATERIAL_SUFFIXES:
        if lowered.endswith(suffix):
            return leading_token[: -len(suffix)]
    return leading_token


def get_patient_id(bam_path: str) -> str:
    """
    Returns the patient ID, read from the directory above the sample directory.

    Falls back to the sample ID, treating one sample as one patient, when the
    path is not nested that way.
    """
    parents = Path(bam_path).resolve().parents
    if len(parents) >= 2:
        return parents[1].name
    return get_sample_id(bam_path)


gwf = Workflow()

# Site-specific settings. Fill these in before running.
slurm_account = ""
bam_file_list_path = "/path/to/bam_paths.txt"  # text file listing one BAM path per line
reference_file = "/path/to/hg38.fasta"
# Optional ATAC add-on. Leave as None unless config.json enables
# fragment_atac_covariates, in which case point this at the chromatin atlas root.
chromatin_dir = None
genomic_blacklist_file = "/path/to/genomic_blacklist.bed"

# Cohort-specific blacklists. These stay disabled unless the matching
# `filters.*.enabled` switch in config.json is turned on.
donor_blacklist_file = None
snp_blacklist_path = None
lowpower_blacklist_path = None
problematic_blacklist_path = None

config_path = str(Path(__file__).resolve().parents[2] / "config.json")

# Prepare output directory once
output_root = Path("../../data")
output_root.mkdir(parents=True, exist_ok=True)

bam_paths = read_text_file(bam_file_list_path)

processed_samples: set[str] = set()

# Process BAM files
for bam_path in bam_paths:

    patient_id = get_patient_id(bam_path)
    sample_id = get_sample_id(bam_path)
    if sample_id in processed_samples:
        continue
    processed_samples.add(sample_id)

    # One sample directory per BAM
    sample_dir = output_root / sample_id

    # Resolve patient-specific SNP blacklist. The extractor expects a directory
    # holding {patient_id}/{patient_id}_snvs_calls.tsv and _indels_calls.tsv.
    patient_snp = None
    if snp_blacklist_path:
        base = Path(snp_blacklist_path)
        patient_dir = resolve_patient_dir(base, patient_id)

        if patient_dir:
            patient_snp = str(base)
        else:
            print(f"[WARN] No SNP directory for {patient_id} under {base}", flush=True)

    # Resolve patient-specific lowpower blacklist
    patient_lowpower = None
    if lowpower_blacklist_path:
        base = Path(lowpower_blacklist_path)
        patient_dir = resolve_patient_dir(base, patient_id)

        if patient_dir:
            lp_candidate = patient_dir / f"{patient_dir.name}.lowpower.lt100X.bg"
            if lp_candidate.exists():
                patient_lowpower = str(lp_candidate)
            else:
                print(f"[WARN] No lowpower blacklist for {patient_id}: {lp_candidate}", flush=True)
        else:
            print(f"[WARN] No lowpower directory for {patient_id} under {base}", flush=True)

    # Resolve patient-specific problematic blacklist
    patient_problematic = None
    if problematic_blacklist_path:
        base = Path(problematic_blacklist_path)
        patient_dir = resolve_patient_dir(base, patient_id)

        if patient_dir:
            matches = sorted(patient_dir.glob(f"{patient_dir.name}.problematic.merged.gt*X.bg"))

            if len(matches) == 1:
                patient_problematic = str(matches[0])
            elif len(matches) > 1:
                patient_problematic = str(matches[0])
                print(
                    f"[WARN] Multiple problematic blacklist files for {patient_id}, "
                    f"using first: {matches}",
                    flush=True,
                )
            else:
                print(
                    f"[WARN] No problematic blacklist for {patient_id} in {patient_dir}",
                    flush=True,
                )
        else:
            print(
                f"[WARN] No problematic blacklist directory for {patient_id} under {base}",
                flush=True,
            )

    # Pass everything to the template
    template_args = {
        "bam_file_path": bam_path,
        "reference_path": reference_file,
        "config_path": config_path,
        "chromatin_dir": chromatin_dir,
        "genomic_blacklist_file": genomic_blacklist_file,
        "donor_blacklist_file": donor_blacklist_file,
        "snp_blacklist_path": patient_snp,
        "lowpower_blacklist_path": patient_lowpower,
        "problematic_blacklist_path": patient_problematic,
        "sample_id": sample_id,
        "patient_id": patient_id,
        "account": slurm_account,
        "output_dir": str(sample_dir),
        "cores": 12,
        "memory": "50g",
        "walltime": "24:00:00",
        "seed": 0,
    }

    gwf.target_from_template(
        name=f"make_fragments_{sample_id}",
        template=make_fragments(**template_args),
    )
