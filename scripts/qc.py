#!/usr/bin/env python3
"""QC Summary Script for RNA-seq Pipeline.

Parses FastQC output files and alignment logs to generate a summary table
with pass/fail flagging for key quality metrics.

Usage:
    python scripts/qc.py <fastqc_dir> <output_csv>

Arguments:
    fastqc_dir  Directory containing FastQC output folders (*_fastqc/)
    output_csv  Path for the summary CSV output
"""

import csv
import os
import sys
import zipfile
from pathlib import Path

# ---- QC thresholds ---------------------------------------------------------
THRESHOLDS = {
    "total_sequences_min": 1_000_000,
    "percent_gc_max": 60,
    "percent_gc_min": 35,
    "percent_duplicates_max": 60,
    "mean_quality_min": 28,
}


def parse_fastqc_data(fastqc_path: str) -> dict:
    """Parse a fastqc_data.txt file and extract key metrics."""
    metrics = {
        "sample": "",
        "total_sequences": 0,
        "percent_gc": 0,
        "percent_duplicates": 0.0,
        "mean_quality": 0.0,
        "sequence_length": "",
    }

    with open(fastqc_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("Filename"):
                metrics["sample"] = line.split("\t")[1]
            elif line.startswith("Total Sequences"):
                metrics["total_sequences"] = int(line.split("\t")[1])
            elif line.startswith("%GC"):
                metrics["percent_gc"] = int(line.split("\t")[1])
            elif line.startswith("#Total Deduplicated Percentage"):
                dedup_pct = float(line.split("\t")[1])
                metrics["percent_duplicates"] = round(100.0 - dedup_pct, 2)
            elif line.startswith("Sequence length"):
                metrics["sequence_length"] = line.split("\t")[1]

    return metrics


def extract_fastqc_data(fastqc_dir: str) -> str:
    """If a zip archive exists, extract fastqc_data.txt from it."""
    data_file = os.path.join(fastqc_dir, "fastqc_data.txt")
    if os.path.isfile(data_file):
        return data_file

    zip_path = fastqc_dir + ".zip"
    if os.path.isfile(zip_path):
        dirname = os.path.basename(fastqc_dir)
        with zipfile.ZipFile(zip_path, "r") as zf:
            target = f"{dirname}/fastqc_data.txt"
            zf.extract(target, os.path.dirname(fastqc_dir))
        return data_file

    return ""


def flag_qc(metrics: dict) -> str:
    """Return PASS or FAIL based on thresholds."""
    reasons = []

    if metrics["total_sequences"] < THRESHOLDS["total_sequences_min"]:
        reasons.append("low_reads")
    if metrics["percent_gc"] > THRESHOLDS["percent_gc_max"]:
        reasons.append("high_gc")
    if metrics["percent_gc"] < THRESHOLDS["percent_gc_min"]:
        reasons.append("low_gc")
    if metrics["percent_duplicates"] > THRESHOLDS["percent_duplicates_max"]:
        reasons.append("high_dup")

    if reasons:
        return "FAIL:" + ";".join(reasons)
    return "PASS"


def main():
    if len(sys.argv) < 3:
        print("Usage: python scripts/qc.py <fastqc_dir> <output_csv>", file=sys.stderr)
        sys.exit(1)

    fastqc_dir = sys.argv[1]
    output_csv = sys.argv[2]

    if not os.path.isdir(fastqc_dir):
        print(f"Error: directory not found: {fastqc_dir}", file=sys.stderr)
        sys.exit(1)

    # Find all FastQC output folders
    qc_dirs = sorted(
        p for p in Path(fastqc_dir).iterdir() if p.is_dir() and p.name.endswith("_fastqc")
    )

    if not qc_dirs:
        print("Warning: no *_fastqc/ directories found, checking for zip files")
        qc_dirs = sorted(
            Path(str(p).replace(".zip", "")) for p in Path(fastqc_dir).glob("*_fastqc.zip")
        )

    results = []
    for qc_path in qc_dirs:
        data_file = extract_fastqc_data(str(qc_path))
        if not data_file:
            print(f"Warning: no data found for {qc_path.name}, skipping")
            continue
        metrics = parse_fastqc_data(data_file)
        metrics["qc_flag"] = flag_qc(metrics)
        results.append(metrics)

    # Write summary CSV
    fieldnames = [
        "sample",
        "total_sequences",
        "sequence_length",
        "percent_gc",
        "percent_duplicates",
        "qc_flag",
    ]

    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    print(f"QC summary written to {output_csv} ({len(results)} samples)")

    # Print summary to stdout
    passed = sum(1 for r in results if r["qc_flag"] == "PASS")
    failed = len(results) - passed
    print(f"  PASS: {passed}  |  FAIL: {failed}")


if __name__ == "__main__":
    main()
