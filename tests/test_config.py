"""The one check on config.json that the Rust test suite cannot make.

`shipped_repository_config_validates` in the Rust suite loads this exact file
through the real validator, so every rule the validator enforces is already
covered there. Restating those rules here would just be a second copy that
drifts out of step with the first.

What the validator cannot catch is a *re-introduced* block: `RootConfig` is not
`deny_unknown_fields`, so an unknown top-level key is silently ignored rather
than rejected. Pasting a `storage` or `patient_covariates` block back in would
look like it worked while doing nothing at all. That is what this file guards,
and it is the only reason it exists.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

CONFIG_PATH = Path(__file__).resolve().parents[1] / "config.json"

# Blocks bam2frag deliberately does not have. Each is either a frag2vec concern
# (patient/sample covariates, model-side preprocessing) or an internal
# format/version declaration with exactly one possible value.
REMOVED_BLOCKS = [
    "patient_covariates",
    "sample_covariates",
    "preprocessing",
    "quality_statistics",
    "eda",
    "model_input",
    "storage",
    "fragment_bank",
    "array_contract",
    "version",
]


@pytest.fixture(scope="module")
def config() -> dict:
    with CONFIG_PATH.open() as handle:
        return json.load(handle)


@pytest.mark.parametrize("block", REMOVED_BLOCKS)
def test_removed_block_is_not_silently_reintroduced(config: dict, block: str) -> None:
    assert block not in config, (
        f"{block!r} is back in config.json. The Rust validator ignores unknown "
        f"top-level keys, so it would be read by nothing and silently do nothing."
    )
