# Project 1: RNA-seq Differential Expression Pipeline

![CI](https://github.com/adamhoffman2155-hue/project-1-rnaseq-pipeline/actions/workflows/ci.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.11-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Repro](https://img.shields.io/badge/FAIR_DOME_CURE-11%2F14_%7C_N%2FA_%7C_4%2F4-brightgreen)

**Research question:** Which transcriptomic programs drive chemotherapy resistance in gastroesophageal adenocarcinoma?

This is the first project in a [computational biology portfolio](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) built to develop practical fluency in bioinformatics workflows. It grew directly from my M.Sc. thesis at McGill, where I studied chemotherapy response in MSI-high GEA patients and needed a way to identify the transcriptomic signatures behind treatment resistance.

## What It Does

End-to-end bulk RNA-seq pipeline using TCGA-STAD expression data:

1. **Quality control** — FastQC per-sample + MultiQC aggregate report
2. **Read alignment** — STAR aligner to GRCh38 (single-end reads wired in the Snakefile)
3. **Indexing** — SAMtools index on sorted BAMs
4. **Quantification** — featureCounts for gene-level counts
5. **Differential expression** — DESeq2 normalization and statistical testing
6. **Visualization** — volcano, MA, PCA, and heatmap plots of top DE genes
7. **Pathway enrichment** — `scripts/gsea_analysis.R` runs fgsea over MSigDB Hallmark + KEGG gene sets as a standalone post-analysis step (not currently part of `rule all`)

## Methods & Tools

| Category | Tools |
|----------|-------|
| Alignment | STAR |
| Quantification | featureCounts |
| QC | FastQC, MultiQC, SAMtools |
| DE Analysis | DESeq2 (R) |
| Visualization | ggplot2 (volcano / MA / PCA / heatmap) |
| Enrichment | fgsea + MSigDB Hallmark / KEGG |
| Workflow | Snakemake |
| Environment | Docker, Conda |

## Project Structure

```
project-1-rnaseq-pipeline/
├── Snakefile                            # fastqc → align → index → count → DE → plots
├── config.yaml                          # STAR, reference, metadata paths
├── Dockerfile
├── environment.yml
├── scripts/
│   ├── download_reference.sh            # Genome & annotation download
│   ├── qc.py                            # QC summary generation
│   ├── deseq2_analysis.R                # DESeq2 DE (Snakemake script rule)
│   ├── visualization.R                  # Volcano / MA / heatmap / PCA
│   ├── gsea_analysis.R                  # fgsea (standalone CLI post-analysis)
│   ├── benchmark_deseq2_vs_limma.R      # DESeq2 vs limma-voom benchmark
│   └── poc/
│       └── run_poc.py                   # Proof-of-concept runner
├── data/
│   └── metadata.csv                     # Sample sheet (raw FASTQ + reference gitignored)
└── results/
    └── poc/                             # POC outputs (committed)
```

Pipeline outputs (QC reports, BAM files, count matrices, DESeq2 tables, figures, logs) are written under `results/` and `logs/` at runtime and are gitignored.

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

# Run pipeline (alignment → DE → plots)
snakemake --cores 4

# Optional: pathway enrichment post-analysis
Rscript scripts/gsea_analysis.R results/de_analysis/deseq2_results.csv results/enrichment/
```

## Proof of Concept

A minimal end-to-end differential expression run on a real, published gene expression dataset so reviewers can verify the statistical workflow without a TCGA download.

**Dataset:** Van't Veer et al. 2002, *Nature* 415:530 — 70-gene breast cancer prognostic signature (198 samples, 78 gene probes, distant-metastasis endpoint). Accessed via `sksurv.datasets.load_breast_cancer()` so no network or account is required.

**Substitution note:** The full Snakemake pipeline targets TCGA-STAD counts from recount3/GDC, but those hosts are not reachable from this reproducibility sandbox. Van't Veer 2002 is a legitimate substitute: a peer-reviewed published gene-expression dataset with a binary clinical endpoint suitable for a two-class DE test. The same Welch-t + BH-FDR code runs unchanged on any counts matrix.

**What the POC tests:**
- Per-gene Welch t-test between metastasis and non-metastasis patients
- Benjamini-Hochberg FDR correction across all 78 probes
- Volcano plot + top DE table

**Headline numbers** (actual run output):
- Samples: 198 (51 metastasis events, 147 non-metastasis)
- Probes with raw p < 0.05: **13**
- Probes passing BH padj < 0.05: **2** (X202240_at, X203306_s_at)
- Top gene X202240_at: log2FC = +0.55, p = 1.1e-4, padj = 7.6e-3

**Limits:**
- Van't Veer probeset IDs are Affymetrix HG-U133A; biological gene symbols require a probe-to-gene mapping step not performed here
- No pathway enrichment in the POC (gseapy requires online gene-set fetch)
- Small sample size yields only 2 probes significant after multiple-testing correction; the ranked list is the primary POC artifact

**Reproduction:**
```bash
pip install scikit-survival pandas numpy scipy matplotlib
python scripts/poc/run_poc.py
```
Outputs are written to `results/poc/` (CSV tables, summary text, volcano PNG).

## My Role

I defined the biological question based on my thesis work, selected the TCGA-STAD dataset, and reviewed DEG outputs for pathway plausibility against known GEA biology. Implementation was heavily AI-assisted.

## Context in the Portfolio

This is **Project 1 of 7** in a portfolio that follows a single clinical question from transcriptomics through single-cell analysis, pharmacogenomics, biomarker discovery, phenomics, survival modeling, and software engineering. See the [portfolio site](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) for the full narrative.

### Cross-project data flow

```
Project 1 (this one — bulk RNA-seq DE + GSEA)
        │   DEGs, pathway scores
        ▼
┌───────────────┬───────────────┬───────────────┐
│ Project 3     │ Project 4     │ Project 6     │
│ (drug ML      │ (DDR biomarkers│ (survival     │
│  features)    │  pathway ctx) │  covariates)  │
└───────────────┴───────────────┴───────────────┘
```

- **Upstream** — TCGA-STAD FASTQs / pre-aligned BAMs (external).
- **Downstream** — gene-level counts and DE tables feed candidate ML features in P3, pathway context for P4's DDR biomarker interpretation, and transcriptomic covariates for P6's Cox survival model (narrative input).

## Benchmarks

| Benchmark | Script | Summary |
| --- | --- | --- |
| DESeq2 vs limma-voom | [`scripts/benchmark_deseq2_vs_limma.R`](scripts/benchmark_deseq2_vs_limma.R) | Runs both tools on the same count matrix and writes per-gene padj/log2FC plus concordance columns. Industry-standard DE-tool sanity check (nf-core/rnaseq reports > 85% inter-tool concordance). |

Run locally:

```
Rscript scripts/benchmark_deseq2_vs_limma.R --help
```

See [`results/benchmark/README.md`](results/benchmark/README.md) for usage.

## Reproducibility

See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for the FAIR-BioRS / CURE self-scorecard (11/14 · N/A · 4/4). DOME is not applicable here (statistical pipeline, not supervised ML).

## License

MIT

## Author

Adam Hoffman — M.Sc. Cancer Research, McGill University
