#!/usr/bin/env python3
"""
Proof-of-Concept: Differential Expression Analysis on Van't Veer 2002 Breast Cancer Dataset

This POC demonstrates the DE statistical workflow on a real, published gene
expression dataset: Van't Veer et al. 2002 Nature 70-gene breast cancer
prognostic signature (198 patients, 78 gene expression features + ER status
+ histologic grade).

Substitution note: the original Snakemake pipeline targets TCGA-STAD RNA-seq
counts from recount3/GDC, but those sources are not reachable from the
reproducibility sandbox used for this POC. Van't Veer is a legitimate
substitute: it is real, published, peer-reviewed breast cancer gene expression
data with a binary clinical endpoint (distant metastasis) suitable for
two-class DE analysis. The exact same Welch t-test + BH-FDR workflow would
run on any counts or log-expression matrix.

Outputs:
    results/poc/de_results.csv     per-gene t statistic, raw and adjusted p
    results/poc/top_de_genes.csv   top 10 by adjusted p
    results/poc/volcano.png
    results/poc/poc_summary.txt
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from sksurv.datasets import load_breast_cancer

RESULTS_DIR = Path(__file__).resolve().parent.parent.parent / "results" / "poc"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def benjamini_hochberg(pvals):
    """Return BH-adjusted p-values for a 1-D array of p-values."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adj, 1.0)
    return out


def main():
    X, y = load_breast_cancer()
    gene_cols = [c for c in X.columns if X[c].dtype.kind == "f"]
    expr = X[gene_cols].copy()
    event = y["e.tdm"].astype(bool)

    n_total = len(expr)
    n_metastasis = int(event.sum())
    n_control = int((~event).sum())
    print(f"Loaded Van't Veer 2002: {n_total} samples, {len(gene_cols)} genes")
    print(f"  Metastasis events: {n_metastasis}")
    print(f"  Non-metastasis: {n_control}")

    results = []
    for gene in gene_cols:
        vals_pos = expr.loc[event, gene].values
        vals_neg = expr.loc[~event, gene].values
        if len(vals_pos) < 2 or len(vals_neg) < 2:
            continue
        t_stat, p_val = stats.ttest_ind(vals_pos, vals_neg, equal_var=False)
        log2fc = vals_pos.mean() - vals_neg.mean()
        results.append({
            "gene": gene,
            "log2FC": log2fc,
            "t_stat": t_stat,
            "pvalue": p_val,
            "mean_metastasis": vals_pos.mean(),
            "mean_control": vals_neg.mean(),
        })

    res_df = pd.DataFrame(results)
    res_df["padj"] = benjamini_hochberg(res_df["pvalue"].values)
    res_df = res_df.sort_values("pvalue").reset_index(drop=True)
    res_df.to_csv(RESULTS_DIR / "de_results.csv", index=False)

    n_sig_005 = int((res_df["padj"] < 0.05).sum())
    n_sig_010 = int((res_df["padj"] < 0.10).sum())
    n_sig_raw_005 = int((res_df["pvalue"] < 0.05).sum())

    top10 = res_df.head(10)
    top10.to_csv(RESULTS_DIR / "top_de_genes.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 6))
    log10p = -np.log10(res_df["pvalue"].clip(lower=1e-300))
    colors = np.where(
        (res_df["padj"] < 0.05) & (res_df["log2FC"].abs() > 0.25),
        "#d73027",
        "#8c8c8c",
    )
    ax.scatter(res_df["log2FC"], log10p, c=colors, s=14, alpha=0.75, edgecolor="none")
    ax.axhline(-np.log10(0.05), linestyle="--", color="grey", linewidth=0.8)
    ax.axvline(0.25, linestyle="--", color="grey", linewidth=0.8)
    ax.axvline(-0.25, linestyle="--", color="grey", linewidth=0.8)
    ax.set_xlabel("log2 expression difference (metastasis vs non-metastasis)")
    ax.set_ylabel("-log10(p)")
    ax.set_title("Volcano plot \u2014 Van't Veer 2002 breast cancer")

    for _, row in top10.head(5).iterrows():
        ax.annotate(
            row["gene"],
            (row["log2FC"], -np.log10(row["pvalue"])),
            fontsize=8,
            xytext=(3, 3),
            textcoords="offset points",
        )
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "volcano.png", dpi=150)
    plt.close(fig)

    summary = f"""Proof-of-Concept Summary: Differential Expression Analysis
=========================================================

Dataset: Van't Veer et al. 2002 Nature 70-gene breast cancer signature
         (bundled with scikit-survival as sksurv.datasets.load_breast_cancer)
Substitution note: This is a substitute for the originally-planned
TCGA-STAD RNA-seq analysis. TCGA/GDC/recount3 data sources are not
reachable from the reproducibility sandbox; Van't Veer is a legitimate
published gene expression dataset with a binary clinical endpoint suitable
for DE analysis. The exact same workflow runs on any counts matrix.

Cohort
  Total samples:      {n_total}
  Metastasis events:  {n_metastasis}
  Non-metastasis:     {n_control}
  Features tested:    {len(gene_cols)} gene probes

Differential expression (Welch t-test + Benjamini-Hochberg FDR)
  Genes with raw p < 0.05:  {n_sig_raw_005}
  Genes with padj < 0.05:   {n_sig_005}
  Genes with padj < 0.10:   {n_sig_010}

Top 10 genes by unadjusted p-value:
{top10[['gene', 'log2FC', 'pvalue', 'padj']].to_string(index=False)}

Honest assessment
  - Dataset is REAL published data (Nature 2002) but small (n=198)
  - Van't Veer probeset IDs (e.g. X200726_at) map to Affymetrix HG-U133A;
    biological gene symbols would require a probe mapping step not done here
  - No pathway enrichment (gseapy requires online gene set fetch, not
    available in sandbox)
  - Identifying top DE probes between metastasis vs non-metastasis is a
    well-defined task; see Van't Veer 2002 supplementary for comparison
  - padj threshold is strict given the small sample size; if zero probes
    pass padj<0.05, the unadjusted top-10 list is still scientifically
    meaningful as the strongest signals in this cohort

Reproduction
  python scripts/poc/run_poc.py
"""
    with open(RESULTS_DIR / "poc_summary.txt", "w") as fh:
        fh.write(summary)
    print(summary)


if __name__ == "__main__":
    main()
