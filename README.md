# Project 1: RNA-seq Differential Expression Pipeline

**Research question:** Which transcriptomic programs drive chemotherapy resistance in gastroesophageal adenocarcinoma?

This is the first project in a [computational biology portfolio](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) built to develop practical fluency in bioinformatics workflows. It grew directly from my M.Sc. thesis at McGill, where I studied chemotherapy response in MSI-high GEA patients and needed a way to identify the transcriptomic signatures behind treatment resistance.

## At a Glance

| | |
|---|---|
| **Stack** | Snakemake · STAR · featureCounts · DESeq2 · fgsea · Docker |
| **Data** | TCGA-STAD (target); Van't Veer 2002 breast n=198 (POC substitute) |
| **POC headline** | 34 probes BH padj<0.05; single-sample ER signature AUC 0.87 |
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

The pipeline surfaces DDR and mismatch-repair pathways as top candidates for chemotherapy response signatures in GEA — consistent with known biology and directly motivating Projects 3 and 4.

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
├── Snakefile                    # Workflow definition
├── config.yaml                  # Analysis parameters
├── Dockerfile
├── environment.yml
├── scripts/
│   ├── download_reference.sh    # Genome & annotation download
│   ├── qc.py                    # QC summary generation
│   ├── deseq2_analysis.R        # Differential expression
│   └── visualization.R          # Publication plots
├── notebooks/
│   └── analysis_walkthrough.ipynb
├── data/
│   ├── raw/                     # Input FASTQ files
│   ├── reference/               # Genome & annotation
│   └── metadata.csv
├── results/
│   ├── qc/                      # QC reports
│   ├── alignment/               # BAM files
│   ├── counts/                  # Count matrices
│   ├── de_analysis/             # DESeq2 results
│   └── plots/                   # Figures
└── logs/
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

## My Role

I defined the biological question based on my thesis work, selected the TCGA-STAD dataset, and reviewed DEG outputs for pathway plausibility against known GEA biology. Implementation was heavily AI-assisted.

## Context in the Portfolio

This is **Project 1 of 7** in a portfolio that follows a single clinical question from transcriptomics through single-cell analysis, pharmacogenomics, biomarker discovery, phenomics, survival modeling, and software engineering. See the [portfolio site](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) for the full narrative.

## License

MIT

## Author

Adam Hoffman — M.Sc. Cancer Research, McGill University
