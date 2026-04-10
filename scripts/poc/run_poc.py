#!/usr/bin/env python3
"""
Proof-of-Concept v2: Differential Expression on Van't Veer 2002 Breast Cancer

Two well-powered comparisons that exercise the Welch t-test + BH-FDR + pathway
scoring pipeline on real published gene expression data:

  1. ER+ vs ER-  (134 vs 64) -- the strongest biological contrast in breast
     cancer, expected to show tens of significant probes.
  2. Poorly vs well differentiated grade (83 vs 30) -- secondary contrast,
     correlated with ER but biologically distinct.

Also runs a single-sample pathway signature score derived from the top DEGs,
and reports the Mann-Whitney AUC of that score for recovering ER status
(a training-fold sanity check).

Substitution note: Original pipeline targets TCGA-STAD RNA-seq counts
(recount3/GDC), which is not reachable from the reproducibility sandbox.
Van't Veer 2002 is a landmark published breast-cancer gene-expression
cohort bundled with scikit-survival, and is scientifically defensible as
a substitute for validating the statistical pipeline.

Outputs:
    results/poc/de_er.csv          DEGs ER+ vs ER-
    results/poc/de_grade.csv       DEGs grade
    results/poc/top_de_er.csv      Top 10 by padj
    results/poc/pathway_scores.csv Signature score per sample
    results/poc/volcano_er.png
    results/poc/volcano_grade.png
    results/poc/er_score_separation.png
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
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adj, 1.0)
    return out


def run_de(expr, labels, name):
    results = []
    for g in expr.columns:
        vals_pos = expr.loc[labels, g].values
        vals_neg = expr.loc[~labels, g].values
        if len(vals_pos) < 2 or len(vals_neg) < 2:
            continue
        t, p = stats.ttest_ind(vals_pos, vals_neg, equal_var=False)
        lfc = vals_pos.mean() - vals_neg.mean()
        results.append({
            "gene": g, "log2FC": lfc, "t_stat": t, "pvalue": p,
            "mean_pos": vals_pos.mean(), "mean_neg": vals_neg.mean(),
        })
    df = pd.DataFrame(results)
    df["padj"] = benjamini_hochberg(df["pvalue"].values)
    df = df.sort_values("pvalue").reset_index(drop=True)
    df["contrast"] = name
    return df


def pathway_signature_scores(expr, de_df, top_n=10):
    up = de_df[de_df["log2FC"] > 0].head(top_n)["gene"].tolist()
    dn = de_df[de_df["log2FC"] < 0].head(top_n)["gene"].tolist()
    z = (expr - expr.mean()) / expr.std()
    sig_up = z[up].mean(axis=1) if up else pd.Series(0, index=expr.index)
    sig_dn = z[dn].mean(axis=1) if dn else pd.Series(0, index=expr.index)
    return sig_up - sig_dn, up, dn


def volcano_plot(de_df, title, out_path):
    fig, ax = plt.subplots(figsize=(8, 6))
    log10p = -np.log10(de_df["pvalue"].clip(lower=1e-300))
    sig_mask = (de_df["padj"] < 0.05) & (de_df["log2FC"].abs() > 0.25)
    ax.scatter(de_df.loc[~sig_mask, "log2FC"], log10p[~sig_mask],
               c="#8c8c8c", s=14, alpha=0.5, edgecolor="none", label="n.s.")
    ax.scatter(de_df.loc[sig_mask, "log2FC"], log10p[sig_mask],
               c="#d73027", s=20, alpha=0.85, edgecolor="none",
               label=f"padj<0.05 & |log2FC|>0.25 (n={int(sig_mask.sum())})")
    ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=0.8)
    ax.axvline(0.25, ls="--", c="grey", lw=0.8)
    ax.axvline(-0.25, ls="--", c="grey", lw=0.8)
    ax.set_xlabel("log2 expression difference")
    ax.set_ylabel("-log10(p)")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=9)
    for _, row in de_df.head(5).iterrows():
        ax.annotate(row["gene"], (row["log2FC"], -np.log10(row["pvalue"])),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    X, y = load_breast_cancer()
    gene_cols = [c for c in X.columns if X[c].dtype.kind == "f"]
    expr = X[gene_cols]
    n_total = len(expr)

    er_pos = (X["er"] == "positive").values
    n_er_pos, n_er_neg = int(er_pos.sum()), int((~er_pos).sum())
    de_er = run_de(expr, er_pos, "ER_pos_vs_neg")
    de_er.to_csv(RESULTS_DIR / "de_er.csv", index=False)
    de_er.head(10).to_csv(RESULTS_DIR / "top_de_er.csv", index=False)
    n_er_005 = int((de_er["padj"] < 0.05).sum())
    n_er_001 = int((de_er["padj"] < 0.01).sum())
    n_er_0001 = int((de_er["padj"] < 0.001).sum())

    grade = X["grade"].values
    keep = (grade == "poorly differentiated") | (grade == "well differentiated")
    is_poor = (grade[keep] == "poorly differentiated")
    de_grade = run_de(expr.loc[keep], is_poor, "poor_vs_well_grade")
    de_grade.to_csv(RESULTS_DIR / "de_grade.csv", index=False)
    n_g_005 = int((de_grade["padj"] < 0.05).sum())
    n_g_001 = int((de_grade["padj"] < 0.01).sum())
    n_poor = int(is_poor.sum())
    n_well = int((~is_poor).sum())

    er_score, up_er, dn_er = pathway_signature_scores(expr, de_er, top_n=10)
    t_val, p_val = stats.ttest_ind(er_score[er_pos], er_score[~er_pos], equal_var=False)
    u_stat, u_p = stats.mannwhitneyu(er_score[er_pos], er_score[~er_pos])
    auc_mw = u_stat / (n_er_pos * n_er_neg)

    scores_df = pd.DataFrame({
        "sample_idx": np.arange(n_total),
        "er_status": X["er"].values,
        "grade": X["grade"].values,
        "er_signature_score": er_score.values,
    })
    scores_df.to_csv(RESULTS_DIR / "pathway_scores.csv", index=False)

    volcano_plot(de_er, "Van't Veer 2002 - ER+ vs ER-",
                 RESULTS_DIR / "volcano_er.png")
    volcano_plot(de_grade, "Van't Veer 2002 - poor vs well-differentiated grade",
                 RESULTS_DIR / "volcano_grade.png")

    fig, ax = plt.subplots(figsize=(6, 5))
    for i, (label, mask) in enumerate([("ER-", ~er_pos), ("ER+", er_pos)]):
        y_vals = er_score[mask].values
        x_jit = np.random.default_rng(i).normal(i, 0.05, size=len(y_vals))
        ax.scatter(x_jit, y_vals, alpha=0.5, s=18, color=["#2b6cb0", "#d73027"][i])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["ER-", "ER+"])
    ax.set_ylabel("ER signature score (top 10 up - top 10 down)")
    ax.set_title(f"ER signature separates ER status\nAUC={auc_mw:.3f}  p={u_p:.2e}")
    ax.axhline(0, c="grey", ls="--", lw=0.8)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "er_score_separation.png", dpi=150)
    plt.close(fig)

    top10_er = de_er[['gene', 'log2FC', 'pvalue', 'padj']].head(10)
    top5_grade = de_grade[['gene', 'log2FC', 'pvalue', 'padj']].head(5)

    summary = f"""Proof-of-Concept v2: Differential Expression + Pathway Scoring
==============================================================

Dataset: Van't Veer et al. 2002, Nature 415:530 - 70-gene breast cancer
         signature (198 patients, 78 gene probes, clinical + survival)
         Bundled with scikit-survival as sksurv.datasets.load_breast_cancer

Substitution note: Original pipeline targets TCGA-STAD RNA-seq counts
(recount3/GDC, not reachable from this reproducibility sandbox).
Van't Veer 2002 is a landmark published breast-cancer gene-expression
cohort; statistical pipeline is agnostic to dataset and runs unchanged
on any counts / log-expression matrix.

============================================================
Contrast 1 - ER+ (n={n_er_pos}) vs ER- (n={n_er_neg})
============================================================
  Genes with raw p < 0.05:   {int((de_er['pvalue'] < 0.05).sum())}
  Genes with BH padj < 0.05: {n_er_005}
  Genes with BH padj < 0.01: {n_er_001}
  Genes with BH padj < 0.001: {n_er_0001}
  Strongest p-value:         {de_er['pvalue'].min():.2e}

Top 10 DE probes (ER+ vs ER-):
{top10_er.to_string(index=False)}

============================================================
Contrast 2 - poorly (n={n_poor}) vs well-differentiated (n={n_well}) grade
============================================================
  Genes with BH padj < 0.05: {n_g_005}
  Genes with BH padj < 0.01: {n_g_001}

Top 5 DE probes (poor vs well):
{top5_grade.to_string(index=False)}

============================================================
Single-sample pathway signature scoring (validation)
============================================================
A 20-probe \"ER signature\" (top-10 up + top-10 down from DEG ranking)
was scored on each sample. Separation of ER+ from ER-:
  Mann-Whitney U AUC:        {auc_mw:.3f}
  Welch t-statistic:         {t_val:.2f}
  p-value:                   {u_p:.2e}

Honest assessment
=================
* REAL published gene expression (Nature 2002) but small (n=198, 78 probes).
* ER+/ER- finds {n_er_005} significant probes at padj<0.05 - strong signal.
* Grade comparison finds {n_g_005} significant probes at padj<0.05.
* The ER signature AUC of {auc_mw:.3f} is training-fold (signature derived
  from the same samples) - pipeline validation, not held-out classifier test.
* No external gene sets (Hallmark/KEGG/GO) used; MSigDB/Enrichr unreachable.

Reproduction
============
    pip install scikit-survival pandas numpy scipy matplotlib
    python scripts/poc/run_poc.py
"""
    (RESULTS_DIR / "poc_summary.txt").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
