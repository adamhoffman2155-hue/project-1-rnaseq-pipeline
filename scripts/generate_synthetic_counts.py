#!/usr/bin/env python3
"""
Generate a synthetic featureCounts-format counts file so the downstream
stages (DESeq2, visualization, GSEA) of the pipeline run end-to-end
without needing real FASTQ / STAR / reference data.

Output
------
results/counts/gene_counts.txt : featureCounts 2.x format
    # Program:featureCounts v2.x.x; ...
    Geneid \t Chr \t Start \t End \t Strand \t Length \t <sample>.bam columns...

Signal
------
6 samples (3 control, 3 treated), 1000 genes. 50 "driver" genes are
up- or down-regulated 4-8x between conditions, the rest are noise.
Poisson/negative-binomial-ish dispersion so DESeq2 produces a
non-degenerate set of results.

Usage
-----
    python scripts/generate_synthetic_counts.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

N_GENES = 1000
N_DE = 50
BASE_LAMBDA = 100
RANDOM_STATE = 42

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT = PROJECT_ROOT / "results" / "counts" / "gene_counts.txt"


def load_sample_order(metadata_csv: Path) -> list[str]:
    df = pd.read_csv(metadata_csv)
    return df["sample_id"].tolist()


def generate_counts(samples: list[str], rng: np.random.Generator) -> pd.DataFrame:
    gene_ids = [f"ENSG{i:011d}" for i in range(N_GENES)]
    # Per-gene base expression, log-normal
    base_expr = np.exp(rng.normal(loc=np.log(BASE_LAMBDA), scale=1.0, size=N_GENES))

    # Pick N_DE driver genes; half upregulated, half downregulated in treated.
    driver_idx = rng.choice(N_GENES, size=N_DE, replace=False)
    up_idx = driver_idx[: N_DE // 2]
    down_idx = driver_idx[N_DE // 2 :]

    counts = np.zeros((N_GENES, len(samples)), dtype=int)
    for j, sid in enumerate(samples):
        is_treated = sid.startswith("SRR_treat")
        mu = base_expr.copy()
        if is_treated:
            mu[up_idx] *= rng.uniform(4, 8, size=up_idx.size)
            mu[down_idx] /= rng.uniform(4, 8, size=down_idx.size)
        # Negative-binomial-ish: draw Poisson with per-gene overdispersion
        # simplified as multiplicative lognormal noise.
        sample_mu = mu * np.exp(rng.normal(0, 0.2, size=N_GENES))
        counts[:, j] = rng.poisson(sample_mu)

    df = pd.DataFrame(counts, columns=[f"results/alignment/{s}.bam" for s in samples])
    df.insert(0, "Length", 1500)
    df.insert(0, "Strand", "+")
    df.insert(0, "End", np.arange(N_GENES) * 1500 + 1500)
    df.insert(0, "Start", np.arange(N_GENES) * 1500 + 1)
    df.insert(0, "Chr", "chr1")
    df.insert(0, "Geneid", gene_ids)
    return df


def write_counts(df: pd.DataFrame, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("# Program:featureCounts v2.0.x; synthetic fixture (not real data)\n")
        df.to_csv(fh, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", default=str(PROJECT_ROOT / "data" / "metadata.csv"))
    parser.add_argument("--output", default=str(OUT))
    args = parser.parse_args()

    rng = np.random.default_rng(RANDOM_STATE)
    samples = load_sample_order(Path(args.metadata))
    df = generate_counts(samples, rng)
    write_counts(df, Path(args.output))
    print(f"Wrote {args.output}  ({N_GENES} genes x {len(samples)} samples, "
          f"{N_DE} planted DE genes)")


if __name__ == "__main__":
    main()
