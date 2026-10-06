"""Tests for the bam2frag Python CLI.

These cover the orchestration layer only: executable resolution, input checking
and the mapping from this tool's snake_case flags onto the Rust extractor's
kebab-case ones. Extraction behaviour itself is covered by the Rust test suite.
"""
from __future__ import annotations

import json
import os
import stat
from pathlib import Path

import pytest

from src import bam2frag


@pytest.fixture
def fake_binary(tmp_path: Path) -> Path:
    """An executable file standing in for the compiled Rust extractor."""
    suffix = ".exe" if os.name == "nt" else ""
    binary = tmp_path / f"bam2frag-core{suffix}"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(binary.stat().st_mode | stat.S_IEXEC)
    return binary


@pytest.fixture
def sample_inputs(tmp_path: Path) -> dict[str, Path]:
    """A minimal set of input files that pass the CLI's existence checks."""
    bam = tmp_path / "SAMPLE001.bam"
    bam.write_bytes(b"")
    (tmp_path / "SAMPLE001.bam.bai").write_bytes(b"")

    reference = tmp_path / "hg38.fasta"
    reference.write_bytes(b"")
    (tmp_path / "hg38.fasta.fai").write_bytes(b"")

    config = tmp_path / "config.json"
    config.write_text(json.dumps({"min_mapq": 30}))

    return {"bam": bam, "reference": reference, "config": config}


def test_env_override_wins(fake_binary: Path, monkeypatch) -> None:
    monkeypatch.setenv(bam2frag.ENV_OVERRIDE, str(fake_binary))
    assert bam2frag.resolve_executable() == str(fake_binary)


def test_env_override_must_exist(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv(bam2frag.ENV_OVERRIDE, str(tmp_path / "nope"))
    with pytest.raises(FileNotFoundError, match=bam2frag.ENV_OVERRIDE):
        bam2frag.resolve_executable()


def test_missing_binary_names_the_build_command(monkeypatch) -> None:
    monkeypatch.delenv(bam2frag.ENV_OVERRIDE, raising=False)
    monkeypatch.setattr(bam2frag.shutil, "which", lambda _: None)
    monkeypatch.setattr(bam2frag, "REPO_ROOT", Path("/nonexistent"))
    with pytest.raises(FileNotFoundError, match="cargo build --release"):
        bam2frag.resolve_executable()


def test_sample_id_is_required(sample_inputs) -> None:
    """A BAM name carries an assay suffix (`_cfdna`) that is not the sample ID,
    so the tool asks rather than guessing and mislabelling the output."""
    with pytest.raises(SystemExit):
        bam2frag.parse_args([
            str(sample_inputs["bam"]),
            str(sample_inputs["reference"]),
            str(sample_inputs["config"]),
        ])


def test_resolved_binary_name_cannot_collide_with_this_cli() -> None:
    """A PATH lookup must never resolve back to this script's console entry point."""
    assert bam2frag.RUST_EXECUTABLE != "bam2frag"


def test_missing_bam_index_is_rejected(tmp_path: Path) -> None:
    bam = tmp_path / "unindexed.bam"
    bam.write_bytes(b"")
    with pytest.raises(FileNotFoundError, match="No BAM index"):
        bam2frag.require_bam_index(bam)


@pytest.mark.parametrize("index_name", ["indexed.bam.bai", "indexed.bai", "indexed.bam.csi"])
def test_accepted_bam_index_names(tmp_path: Path, index_name: str) -> None:
    bam = tmp_path / "indexed.bam"
    bam.write_bytes(b"")
    (tmp_path / index_name).write_bytes(b"")
    assert bam2frag.require_bam_index(bam).name == index_name


def test_default_threads_follows_slurm(monkeypatch) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    assert bam2frag.default_threads() == 8
    monkeypatch.delenv("SLURM_CPUS_PER_TASK")
    assert bam2frag.default_threads() == 1


def _capture_extract_command(monkeypatch, sample_inputs, fake_binary, extra_argv):
    """Run main() with the subprocess layer stubbed, returning the extract argv."""
    calls: list[list[str]] = []

    monkeypatch.setenv(bam2frag.ENV_OVERRIDE, str(fake_binary))

    class Completed:
        returncode = 0
        stdout = ""
        stderr = ""

    monkeypatch.setattr(bam2frag.subprocess, "run", lambda cmd, **kw: calls.append(cmd) or Completed())
    # main() checks for the completion marker the real extractor would write.
    monkeypatch.setattr(bam2frag.Path, "is_file", lambda self: True)

    argv = [
        str(sample_inputs["bam"]),
        str(sample_inputs["reference"]),
        str(sample_inputs["config"]),
        *([] if "--sample_id" in extra_argv else ["--sample_id", "SAMPLE001"]),
        *extra_argv,
    ]
    bam2frag.main(argv)
    return [call for call in calls if "extract" in call][0]


def test_extract_command_maps_flags(monkeypatch, sample_inputs, fake_binary, tmp_path) -> None:
    command = _capture_extract_command(
        monkeypatch, sample_inputs, fake_binary,
        ["--output_dir", str(tmp_path / "out"), "--threads", "4", "--seed", "7", "--overwrite"],
    )

    assert command[1] == "extract"
    assert "--sample-id" in command and "SAMPLE001" in command
    # patient_id defaults to sample_id, treating one sample as one patient.
    assert command[command.index("--patient-id") + 1] == "SAMPLE001"
    assert command[command.index("--threads") + 1] == "4"
    assert command[command.index("--seed") + 1] == "7"
    assert "--overwrite" in command

    # Covariate tables are not part of bam2frag and must never be passed.
    for absent in ["--clinical-info", "--qcstats", "--sample-metadata"]:
        assert absent not in command


def test_optional_resources_are_omitted_when_unset(
    monkeypatch, sample_inputs, fake_binary, tmp_path
) -> None:
    command = _capture_extract_command(
        monkeypatch, sample_inputs, fake_binary, ["--output_dir", str(tmp_path / "out")]
    )
    # An unset optional path must not render as a bare flag, which would make the
    # extractor read the next flag as its value.
    for absent in [
        "--chromatin-dir",
        "--genomic-blacklist",
        "--donor-blacklist",
        "--snp-blacklist",
        "--lowpower-blacklist",
        "--problematic-blacklist",
    ]:
        assert absent not in command


def test_explicit_sample_and_patient_ids(
    monkeypatch, sample_inputs, fake_binary, tmp_path
) -> None:
    command = _capture_extract_command(
        monkeypatch, sample_inputs, fake_binary,
        ["--output_dir", str(tmp_path / "out"), "--sample_id", "S9", "--patient_id", "P9"],
    )
    assert command[command.index("--sample-id") + 1] == "S9"
    assert command[command.index("--patient-id") + 1] == "P9"


def test_chromosome_batch_size_defaults_to_twelve(
    monkeypatch, sample_inputs, fake_binary, tmp_path
) -> None:
    command = _capture_extract_command(
        monkeypatch, sample_inputs, fake_binary, ["--output_dir", str(tmp_path / "out")]
    )
    assert command[command.index("--chromosome-batch-size") + 1] == "12"


def test_chromosome_batch_size_is_passed_through(
    monkeypatch, sample_inputs, fake_binary, tmp_path
) -> None:
    # The memory lever: lowering it is how a whole-genome run fits a smaller
    # allocation, so it must reach the extractor rather than be silently dropped.
    command = _capture_extract_command(
        monkeypatch, sample_inputs, fake_binary,
        ["--output_dir", str(tmp_path / "out"), "--chromosome_batch_size", "3"],
    )
    assert command[command.index("--chromosome-batch-size") + 1] == "3"
