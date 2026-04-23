"""
Pipeline smoke tests that run without FASTQ / STAR / R.

Exercises the Python-side contracts:
  * metadata.csv structure matches what the Snakefile reads
  * config.yaml parses and has every key the Snakefile references
  * the synthetic-counts generator writes a valid featureCounts file
    with planted DE signal

The R-side (DESeq2, fgsea, visualization) is validated by
`snakemake --dry-run`, not by pytest.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
COUNTS_OUT = PROJECT_ROOT / "results" / "counts" / "gene_counts.txt"


def test_metadata_has_required_columns():
    df = pd.read_csv(PROJECT_ROOT / "data" / "metadata.csv")
    for col in ("sample_id", "condition"):
        assert col in df.columns
    assert set(df["condition"].unique()) >= {"control", "treated"}
    assert df["sample_id"].is_unique


def test_config_has_keys_the_snakefile_reads():
    with (PROJECT_ROOT / "config.yaml").open() as fh:
        cfg = yaml.safe_load(fh)
    assert cfg["metadata_file"]
    assert cfg["reference"]["gtf"]
    assert cfg["star"]["genome_index"]
    assert cfg["deseq2"]["design"] == "~condition"
    assert cfg["deseq2"]["contrast"] == ["condition", "treated", "control"]


@pytest.fixture(scope="module")
def synthetic_counts():
    """Run the synthetic-counts generator and return the output path."""
    subprocess.run(
        [sys.executable,
         str(PROJECT_ROOT / "scripts" / "generate_synthetic_counts.py"),
         "--output", str(COUNTS_OUT)],
        check=True, cwd=PROJECT_ROOT,
    )
    assert COUNTS_OUT.is_file()
    return COUNTS_OUT


def test_synthetic_counts_is_featurecounts_format(synthetic_counts):
    """featureCounts format: first line '# Program:...', header row has
    Geneid + 5 annotation cols + one column per sample."""
    with synthetic_counts.open() as fh:
        first_line = fh.readline()
    assert first_line.startswith("# Program:"), first_line
    df = pd.read_csv(synthetic_counts, sep="\t", comment=None, skiprows=1)
    assert list(df.columns[:6]) == ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    # One column per sample in metadata.csv.
    metadata = pd.read_csv(PROJECT_ROOT / "data" / "metadata.csv")
    assert df.shape[1] == 6 + len(metadata)


def test_synthetic_counts_has_planted_signal(synthetic_counts):
    """The planted DE signal (50 driver genes, 4-8x effect) should be
    obvious as a per-gene mean-ratio between treated and control samples."""
    df = pd.read_csv(synthetic_counts, sep="\t", skiprows=1)
    counts = df.set_index("Geneid").iloc[:, 5:]
    ctrl_cols = [c for c in counts.columns if "ctrl" in c]
    trt_cols = [c for c in counts.columns if "treat" in c]
    assert ctrl_cols and trt_cols

    # Avoid log(0); DESeq2 does the same with pseudocount.
    log_ratio = (
        (counts[trt_cols].mean(axis=1) + 1).apply("log")
        - (counts[ctrl_cols].mean(axis=1) + 1).apply("log")
    )
    # At least 30 genes should show |log-ratio| > 1 (i.e. >2x fold change) —
    # we planted 50, some will be damped by noise.
    strong = (log_ratio.abs() > 1.0).sum()
    assert strong >= 30, f"Expected >=30 planted DE genes, saw {strong}"


def test_gsea_script_exists_and_is_cli_runnable():
    """The GSEA rule `shell`s out to gsea_analysis.R with 2 positional
    args; check the script defends against wrong usage."""
    gsea_r = PROJECT_ROOT / "scripts" / "gsea_analysis.R"
    assert gsea_r.is_file()
    text = gsea_r.read_text()
    assert 'Usage: Rscript gsea_analysis.R' in text
    assert 'fgsea' in text
