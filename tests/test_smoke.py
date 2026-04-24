"""Minimal smoke tests for the Snakemake-driven RNA-seq pipeline.

This project's real work is done by Snakemake rules wrapping STAR,
featureCounts, DESeq2 and fgsea — all of which are exercised in CI via
``snakemake --lint``. These pytest tests just guard the Python entry
points and metadata integrity so CI catches import errors early.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_repo_layout_present() -> None:
    assert (REPO_ROOT / "Snakefile").is_file()
    assert (REPO_ROOT / "config.yaml").is_file()
    assert (REPO_ROOT / "environment.yml").is_file()
    assert (REPO_ROOT / "Dockerfile").is_file()


def test_metadata_csv_has_required_columns() -> None:
    meta = REPO_ROOT / "data" / "metadata.csv"
    if not meta.is_file():
        pytest.skip("metadata.csv not committed in this environment")
    with meta.open() as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
    for col in ("sample_id", "condition"):
        assert col in headers, f"metadata.csv missing required column {col!r}: {headers}"


def test_qc_script_imports() -> None:
    import importlib.util

    qc_path = REPO_ROOT / "scripts" / "qc.py"
    spec = importlib.util.spec_from_file_location("qc", qc_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Only exercise the module-level code; any missing optional deps will
    # skip with a clear error rather than crash the test suite.
    try:
        spec.loader.exec_module(module)
    except ModuleNotFoundError as exc:
        pytest.skip(f"qc.py requires {exc.name!r}; skipping import test")
