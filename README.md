# Project 1: RNA-seq Differential Expression Pipeline

**Research question:** Which transcriptomic programs drive chemotherapy resistance in gastroesophageal adenocarcinoma?

This is the first project in a [computational biology portfolio](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) built to develop practical fluency in bioinformatics workflows. It grew directly from my M.Sc. thesis at McGill, where I studied chemotherapy response in MSI-high GEA patients and needed a way to identify the transcriptomic signatures behind treatment resistance.

## At a Glance

| | |
|---|---|
| **Stack** | Snakemake · STAR · featureCounts · DESeq2 · fgsea · Docker |
| **Data** | TCGA-STAD (target); Van't Veer 2002 breast n=198 (POC substitute) |
| **POC headline** | Welch-t + BH-FDR over 78 probes; 2 probes padj<0.05 (X202240_at top, p=1.1e-4) |
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

# Run pipeline
snakemake --cores 4
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

## License

MIT

## Author

Adam Hoffman — M.Sc. Cancer Research, McGill University
