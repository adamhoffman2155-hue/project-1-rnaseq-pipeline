# Project 1: RNA-seq Differential Expression Pipeline

> **A pipeline that finds which genes are "turned up" or "turned down" when cancer patients resist chemotherapy.**

## The short version

**What this project does.** Builds a reproducible workflow that takes raw RNA-sequencing data and finds genes whose activity differs between two groups of patients (for example, responders vs. non-responders to chemo).

**The question behind it.** Some stomach/esophageal cancer patients don't respond to standard chemotherapy. Is there a genetic fingerprint we can read from their tumors that would flag who will resist? To answer that honestly, you first need a pipeline that is actually reproducible — meaning someone else can re-run it and get the same numbers.

**What the proof-of-concept shows.** Because real patient RNA-seq data from TCGA isn't accessible from the reproducibility sandbox, the POC runs the same statistical workflow on a well-known 1998 breast-cancer dataset (Van't Veer, n=198, 78 genes). It correctly finds **34 genes** whose activity differs between ER-positive and ER-negative tumors (top hit: p = 1.7 × 10⁻¹²), and a 20-gene signature separates ER status with **0.87 AUC** (training-fold; see honesty note below).

**Why it matters.** If the same approach finds robust signatures in GEA patients, clinicians could potentially identify chemo-resistant tumors before treatment and avoid putting those patients through ineffective therapy.

---

_The rest of this README is technical detail for bioinformaticians, recruiters doing a deep review, or anyone reproducing the work._

## At a Glance

| | |
|---|---|
| **Stack** | Snakemake · STAR · featureCounts · DESeq2 · fgsea · Docker (full pipeline); scipy + matplotlib (POC) |
| **Data** | TCGA-STAD (full-pipeline target); Van't Veer 2002 breast cancer n=198, 78 probes (POC substitute) |
| **POC headline** | ER+ vs ER- DE: 34 probes BH padj<0.05 (15 at padj<0.001, top p=1.7e-12); 20-probe ER signature MW-AUC 0.868 (training-fold) |
| **Status** | Full pipeline: **Full-data target** (no FASTQ committed). POC: **Runnable POC** (reproducible with one pip install) |
| **Role** | Biological framing and DE interpretation; implementation AI-assisted |
| **Portfolio** | Project 1 of 7 · [full narrative](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) |

## What It Does

End-to-end bulk RNA-seq pipeline using TCGA-STAD expression data:

1. **Read alignment** — STAR aligner to GRCh38 reference genome
2. **Quantification** — featureCounts for gene-level counts
3. **Quality control** — FastQC + MultiQC reporting
4. **Differential expression** — DESeq2 normalization and statistical testing
5. **Visualization** — Volcano plots, MA plots, heatmaps of top DE genes
6. **Pathway enrichment** — GSEA/fgsea for biological interpretation

## Methods & Tools

| Category | Tools |
|----------|-------|
| Alignment | STAR |
| Quantification | featureCounts |
| QC | FastQC, MultiQC, SAMtools |
| DE Analysis | DESeq2 (R) |
| Visualization | ggplot2, ComplexHeatmap |
| Enrichment | fgsea / clusterProfiler |
| Workflow | Snakemake |
| Environment | Docker, Conda |

## Project Structure

```
project-1-rnaseq-pipeline/
├── README.md
├── Snakefile
├── config.yaml
├── Dockerfile
├── environment.yml
├── data/
│   └── metadata.csv
├── scripts/
│   ├── deseq2_analysis.R
│   ├── download_reference.sh
│   ├── gsea_analysis.R
│   ├── qc.py
│   ├── visualization.R
│   └── poc/
│       └── run_poc.py          # Proof-of-concept runner
└── results/
    └── poc/                    # POC outputs (committed)
```

## Quick Start

```bash
git clone https://github.com/adamhoffman2155-hue/project-1-rnaseq-pipeline.git
cd project-1-rnaseq-pipeline

# Using Docker
docker build -t rnaseq-pipeline .
docker run -it -v $(pwd):/workspace rnaseq-pipeline bash

# Or Conda
conda env create -f environment.yml
conda activate rnaseq-pipeline

# Run pipeline (requires FASTQ files in data/raw/)
snakemake --cores 4
```

## Proof of Concept (v2)

A lightweight, reproducible DE + pathway-scoring run on a real, published gene-expression dataset — so reviewers can verify the statistical pipeline without running the full alignment/quantification stack.

**Dataset:** Van't Veer et al. 2002, *Nature* 415:530 — 198 breast cancer patients, 78 gene probes, ER status, grade, metastasis endpoint. Accessed via `sksurv.datasets.load_breast_cancer()` so **no network or account is required**.

**Substitution note:** The full Snakemake pipeline targets TCGA-STAD counts (recount3 / GDC), but those hosts aren't reachable from the reproducibility sandbox. The Van't Veer 2002 cohort is a landmark published gene-expression dataset, and the Welch-t + BH-FDR + signature-score code runs unchanged on any expression matrix.

### What the POC runs

1. **Contrast 1 — ER+ (n=134) vs ER- (n=64):** Welch t-test per probe, Benjamini–Hochberg FDR
2. **Contrast 2 — poorly (n=83) vs well-differentiated (n=30) grade:** same test
3. **Pathway signature scoring:** 10 top-up + 10 top-down DEG probes become a "ER signature" score, evaluated on each sample
4. **Sanity check:** Mann-Whitney U / Welch t on signature score between ER+ and ER- (training-fold AUC)

### Headline numbers (actual committed outputs in `results/poc/`)

**ER+ vs ER-:**

| Threshold | N significant probes |
|---|---|
| raw p < 0.05 | 39 |
| BH padj < 0.05 | **34** |
| BH padj < 0.01 | 23 |
| BH padj < 0.001 | 15 |
| Strongest p | 1.74e-12 |

Top 5 probes by padj: X214919_s_at, X202240_at, X211762_s_at, X204540_at, X211040_x_at.

**Grade (poor vs well-differentiated):**

| Threshold | N significant probes |
|---|---|
| BH padj < 0.05 | 23 |
| BH padj < 0.01 | 17 |

**20-probe ER signature scoring:**

| Metric | Value |
|---|---|
| Mann-Whitney U AUC | **0.868** |
| Welch t-statistic | 10.49 |
| p-value | 6.16e-17 |

### Reproduce

```bash
pip install scikit-survival pandas numpy scipy matplotlib
python scripts/poc/run_poc.py
```

Outputs land in `results/poc/`:
- `de_er.csv`, `de_grade.csv` — full DE tables
- `top_de_er.csv` — top-10 by padj
- `pathway_scores.csv` — per-sample signature score
- `volcano_er.png`, `volcano_grade.png`, `er_score_separation.png`
- `poc_summary.txt` — full run log

### Honest assessment

- Van't Veer is real published gene expression but small (n=198, 78 probes). ER+/ER- is a strong contrast in breast cancer; 34 significant probes is expected.
- **The ER signature AUC of 0.868 is training-fold only** — the signature was derived from the same samples it was scored on. This is a pipeline-validation sanity check, not a held-out classifier result.
- No external gene sets (Hallmark / KEGG / GO) are used in the POC; MSigDB/Enrichr are unreachable from the sandbox.
- The POC validates the statistical pipeline, not GEA biology. The full Snakemake target (TCGA-STAD) is where GEA-specific signatures would be surfaced.

## My Role

I defined the biological question based on my thesis work, selected the TCGA-STAD dataset, and reviewed DEG outputs for pathway plausibility against known GEA biology. Implementation was heavily AI-assisted.

## Context in the Portfolio

This is **Project 1 of 7** in a portfolio that follows a single clinical question from transcriptomics through single-cell analysis, pharmacogenomics, biomarker discovery, phenomics, survival modeling, and software engineering. See the [portfolio site](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) for the full narrative.

## License

MIT

## Author

Adam Hoffman — M.Sc. Cancer Research, McGill University
