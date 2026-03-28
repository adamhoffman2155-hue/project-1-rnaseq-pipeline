# RNA-seq Differential Expression Analysis Pipeline

A production-grade, reproducible Snakemake workflow for RNA-seq data analysis from raw FASTQ files to publication-ready differential expression results.

## Overview

This pipeline demonstrates a complete RNA-seq analysis workflow:
1. **Quality Control** — FastQC, MultiQC
2. **Read Alignment** — STAR aligner to reference genome
3. **Quantification** — featureCounts for gene expression
4. **Normalization & DE Analysis** — DESeq2 in R
5. **Visualization** — Volcano plots, heatmaps, MA plots
6. **Pathway Enrichment** — GSEA analysis

## Skills Demonstrated

✅ **Bioinformatics Tools:** STAR, featureCounts, FastQC, SAMtools  
✅ **Workflow Automation:** Snakemake orchestration, rule dependencies  
✅ **Python:** pandas, NumPy, data preprocessing, automation scripts  
✅ **R:** tidyverse, DESeq2, ggplot2, data visualization  
✅ **Bash/Linux:** Command-line tools, shell scripting, file manipulation  
✅ **Containerization:** Docker for reproducible environments  
✅ **Reproducibility:** Conda environment files, documented workflows  
✅ **Genomic Data:** FASTQ, BAM, GTF file handling  

## Quick Start

### Prerequisites
- Docker (recommended) or Conda
- Git
- 20+ GB disk space for reference genomes and results

### Installation

```bash
# Clone repository
git clone https://github.com/adamhoffman2155-hue/project-1-rnaseq-pipeline.git
cd project-1-rnaseq-pipeline

# Option 1: Using Docker (recommended)
docker build -t rnaseq-pipeline .
docker run -it -v $(pwd):/workspace rnaseq-pipeline bash

# Option 2: Using Conda
conda env create -f environment.yml
conda activate rnaseq-pipeline
```

### Running the Pipeline

```bash
# Configure analysis parameters
vim config.yaml

# Validate workflow
snakemake -n

# Run pipeline (locally with 4 cores)
snakemake --cores 4

# Or submit to cluster (SLURM example)
snakemake --profile slurm --jobs 100
```

## Workflow Architecture

```
Raw FASTQ Files
    ↓
[FastQC] → Quality Assessment
    ↓
[STAR] → Genome Alignment
    ↓
[featureCounts] → Gene Quantification
    ↓
[DESeq2] → Normalization & DE Analysis
    ↓
[Visualization] → Plots & Reports
    ↓
Results: DE genes, plots, enrichment analysis
```

## Project Structure

```
project-1-rnaseq-pipeline/
├── README.md                 # This file
├── Snakefile                 # Workflow definition
├── config.yaml               # Configuration parameters
├── Dockerfile                # Container specification
├── environment.yml           # Conda dependencies
├── scripts/
│   ├── download_reference.sh # Download genome & annotation
│   ├── qc.py                 # Quality control summary
│   ├── deseq2_analysis.R     # DESeq2 differential expression
│   └── visualization.R       # Generate publication plots
├── notebooks/
│   └── analysis_walkthrough.ipynb  # Jupyter notebook tutorial
├── data/
│   ├── raw/                  # Input FASTQ files
│   ├── reference/            # Genome & annotation files
│   └── metadata.csv          # Sample information
├── results/
│   ├── qc/                   # FastQC/MultiQC reports
│   ├── alignment/            # BAM files & alignment stats
│   ├── counts/               # Gene count matrices
│   ├── de_analysis/          # DESeq2 results
│   └── plots/                # Publication-quality figures
└── logs/                     # Snakemake execution logs
```

## Configuration

Edit `config.yaml` to customize analysis:

```yaml
# Sample metadata
samples:
  sample_1: data/raw/sample_1_R1.fastq.gz
  sample_2: data/raw/sample_2_R1.fastq.gz

# Reference genome
reference:
  genome: "GRCh38"
  gtf: "data/reference/gencode.v38.annotation.gtf"
  fasta: "data/reference/GRCh38.primary_assembly.genome.fa"

# Analysis parameters
star:
  threads: 8
  memory: 40

deseq2:
  design: "~condition"
  contrast: ["condition", "treated", "control"]
  padj_threshold: 0.05
  lfc_threshold: 1.0
```

## Input Data

### Sample Metadata (`data/metadata.csv`)
```csv
sample_id,condition,replicate
sample_1,control,1
sample_2,control,2
sample_3,treated,1
sample_4,treated,2
```

### FASTQ Files
- Single-end or paired-end reads
- Gzipped format (.fastq.gz)
- Location: `data/raw/`

## Output Files

### Key Results
- `results/de_analysis/deseq2_results.csv` — DESeq2 results (baseMean, log2FC, padj)
- `results/plots/volcano_plot.pdf` — Volcano plot (log2FC vs -log10(padj))
- `results/plots/heatmap_top_genes.pdf` — Heatmap of top DE genes
- `results/plots/ma_plot.pdf` — MA plot (mean vs log2FC)
- `results/qc/multiqc_report.html` — Quality control summary

### Intermediate Files
- `results/alignment/*.bam` — Aligned reads (BAM format)
- `results/counts/gene_counts.txt` — Raw gene count matrix
- `results/qc/fastqc/` — Per-sample QC reports

## Usage Examples

### Example 1: Analyze provided test data
```bash
# Download test data
bash scripts/download_test_data.sh

# Run pipeline
snakemake --cores 4

# View results
ls -lh results/de_analysis/
```

### Example 2: Analyze your own data
```bash
# 1. Place FASTQ files in data/raw/
cp /path/to/your/fastq/*.fastq.gz data/raw/

# 2. Create metadata.csv with sample information
# 3. Update config.yaml with your parameters
# 4. Run pipeline
snakemake --cores 8
```

### Example 3: Run specific rule
```bash
# Run only quality control
snakemake qc --cores 4

# Run only DESeq2 analysis
snakemake deseq2 --cores 4

# Dry run to see what will execute
snakemake -n
```

## Interpretation of Results

### Volcano Plot
- **X-axis:** log2(fold change) — magnitude of expression change
- **Y-axis:** -log10(adjusted p-value) — statistical significance
- **Red points:** Significantly DE genes (padj < 0.05, |log2FC| > 1)
- **Interpretation:** Points in upper left/right corners are significant DE genes

### MA Plot
- **X-axis:** log10(mean expression) — average abundance
- **Y-axis:** log2(fold change) — direction and magnitude of change
- **Red points:** Significantly DE genes
- **Interpretation:** Shows relationship between expression level and fold change

### Heatmap
- **Rows:** Top DE genes (sorted by padj)
- **Columns:** Samples
- **Color:** log2-normalized expression (blue=low, red=high)
- **Interpretation:** Shows which genes drive differences between conditions

## Troubleshooting

### Issue: "STAR index not found"
```bash
# Download and build STAR index
bash scripts/download_reference.sh
```

### Issue: "Insufficient memory"
Edit `config.yaml` and increase `star.memory` or reduce `star.threads`

### Issue: "Snakemake rule failed"
Check logs:
```bash
cat logs/snakemake.log
# Or specific rule:
cat logs/align_{sample}.log
```

## Performance

**Benchmark** (on 4-core machine with 16 GB RAM):
- 4 samples (paired-end, ~50M reads each)
- Time: ~2-3 hours
- Disk: ~150 GB (includes reference genome)

## Reproducibility

This pipeline ensures reproducibility through:
- **Conda environment** — locked dependency versions
- **Docker container** — isolated runtime environment
- **Snakemake** — deterministic rule execution
- **Seed management** — fixed random seeds in R scripts
- **Documentation** — all parameters logged in config.yaml

To reproduce results:
```bash
# Use exact same environment
conda env create -f environment.yml
conda activate rnaseq-pipeline

# Use exact same config
cat config.yaml

# Re-run pipeline
snakemake --cores 4
```

## References

- **STAR aligner:** Dobin et al. (2013) Bioinformatics
- **featureCounts:** Liao et al. (2014) Bioinformatics
- **DESeq2:** Love et al. (2014) Genome Biology
- **Snakemake:** Köster & Rahmann (2012) Bioinformatics

## License

MIT License — See LICENSE file

## Contact

Adam Hoffman  
Email: adamhoffman21@hotmail.ca  
GitHub: [@adamhoffman2155-hue](https://github.com/adamhoffman2155-hue)

## Acknowledgments

- Test data from [GEO](https://www.ncbi.nlm.nih.gov/geo/)
- Reference genome from [GENCODE](https://www.gencodegenes.org/)
- Workflow inspired by [nf-core/rnaseq](https://github.com/nf-core/rnaseq)

