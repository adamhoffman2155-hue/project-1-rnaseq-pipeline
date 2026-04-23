# Project 1: RNA-seq Differential Expression Pipeline

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

The pipeline surfaces DDR and mismatch-repair pathways as top candidates for chemotherapy response signatures in GEA — consistent with known biology and directly motivating Projects 3 and 4.

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
├── Snakefile                    # fastqc → align → index → count → DE → plots
├── config.yaml                  # STAR, reference, metadata paths
├── Dockerfile
├── environment.yml
├── scripts/
│   ├── download_reference.sh    # Genome & annotation download
│   ├── qc.py                    # QC summary generation
│   ├── deseq2_analysis.R        # DESeq2 DE (Snakemake script rule)
│   ├── visualization.R          # Volcano / MA / heatmap / PCA
│   └── gsea_analysis.R          # fgsea (standalone CLI post-analysis)
├── data/
│   └── metadata.csv             # Sample sheet (raw FASTQ + reference gitignored)
└── .gitignore
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

## My Role

I defined the biological question based on my thesis work, selected the TCGA-STAD dataset, and reviewed DEG outputs for pathway plausibility against known GEA biology. Implementation was heavily AI-assisted.

## Context in the Portfolio

This is **Project 1 of 7** in a portfolio that follows a single clinical question from transcriptomics through single-cell analysis, pharmacogenomics, biomarker discovery, phenomics, survival modeling, and software engineering. See the [portfolio site](https://github.com/adamhoffman2155-hue/bioinformatics-portfolio) for the full narrative.

## License

MIT

## Author

Adam Hoffman — M.Sc. Cancer Research, McGill University
