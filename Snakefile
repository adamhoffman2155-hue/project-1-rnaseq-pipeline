"""
RNA-seq Differential Expression Analysis Pipeline
Snakemake workflow for processing raw FASTQ to DE results
"""

import os
import pandas as pd
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Load sample metadata
samples_df = pd.read_csv(config["metadata_file"])
SAMPLES = samples_df["sample_id"].tolist()
CONDITIONS = dict(zip(samples_df["sample_id"], samples_df["condition"]))

# Wildcard constraints
wildcard_constraints:
    sample="[a-zA-Z0-9_]+",
    read="R[12]"

# Default target
rule all:
    input:
        # Quality control
        expand("results/qc/fastqc/{sample}_R1_fastqc.html", sample=SAMPLES),
        "results/qc/multiqc_report.html",
        
        # Alignment
        expand("results/alignment/{sample}.bam", sample=SAMPLES),
        expand("results/alignment/{sample}.bam.bai", sample=SAMPLES),
        
        # Quantification
        "results/counts/gene_counts.txt",
        
        # Differential expression
        "results/de_analysis/deseq2_results.csv",
        
        # Visualization
        "results/plots/volcano_plot.pdf",
        "results/plots/ma_plot.pdf",
        "results/plots/heatmap_top_genes.pdf",
        "results/plots/pca_plot.pdf"

# ============================================================================
# QUALITY CONTROL
# ============================================================================

rule fastqc:
    """Run FastQC on raw FASTQ files"""
    input:
        fq="data/raw/{sample}_R1.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}_R1_fastqc.html",
        zip="results/qc/fastqc/{sample}_R1_fastqc.zip"
    log:
        "logs/fastqc_{sample}.log"
    threads: 2
    shell:
        """
        fastqc -o results/qc/fastqc/ {input.fq} 2> {log}
        """

rule multiqc:
    """Aggregate FastQC results"""
    input:
        expand("results/qc/fastqc/{sample}_R1_fastqc.zip", sample=SAMPLES)
    output:
        "results/qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    shell:
        """
        multiqc results/qc/fastqc/ -o results/qc/ -n multiqc_report 2> {log}
        """

# ============================================================================
# ALIGNMENT
# ============================================================================

rule star_align:
    """Align reads to reference genome using STAR"""
    input:
        fq="data/raw/{sample}_R1.fastq.gz",
        genome_dir=config["star"]["genome_index"]
    output:
        bam="results/alignment/{sample}.bam",
        log="results/alignment/{sample}_Log.final.out"
    log:
        "logs/align_{sample}.log"
    threads: config["star"]["threads"]
    resources:
        mem_mb=config["star"]["memory"]
    shell:
        """
        STAR \
            --genomeDir {input.genome_dir} \
            --readFilesIn {input.fq} \
            --readFilesCommand zcat \
            --outFileNamePrefix results/alignment/{wildcards.sample}_ \
            --outSAMtype BAM SortedByCoordinate \
            --runThreadN {threads} \
            --outBAMsortingThreadN {threads} \
            2> {log}
        
        # Rename output BAM
        mv results/alignment/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        """

rule samtools_index:
    """Index BAM files"""
    input:
        bam="results/alignment/{sample}.bam"
    output:
        bai="results/alignment/{sample}.bam.bai"
    log:
        "logs/index_{sample}.log"
    threads: 2
    shell:
        """
        samtools index -@ {threads} {input.bam} 2> {log}
        """

# ============================================================================
# QUANTIFICATION
# ============================================================================

rule featurecounts:
    """Count reads per gene using featureCounts"""
    input:
        bams=expand("results/alignment/{sample}.bam", sample=SAMPLES),
        gtf=config["reference"]["gtf"]
    output:
        counts="results/counts/gene_counts.txt",
        summary="results/counts/gene_counts.txt.summary"
    log:
        "logs/featurecounts.log"
    threads: 4
    shell:
        """
        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} \
            2> {log}
        """

# ============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

rule deseq2_analysis:
    """Perform DESeq2 differential expression analysis"""
    input:
        counts="results/counts/gene_counts.txt",
        metadata=config["metadata_file"]
    output:
        results="results/de_analysis/deseq2_results.csv",
        normalized="results/de_analysis/normalized_counts.csv"
    log:
        "logs/deseq2.log"
    script:
        "scripts/deseq2_analysis.R"

# ============================================================================
# VISUALIZATION
# ============================================================================

rule volcano_plot:
    """Generate volcano plot"""
    input:
        results="results/de_analysis/deseq2_results.csv"
    output:
        plot="results/plots/volcano_plot.pdf"
    log:
        "logs/volcano_plot.log"
    script:
        "scripts/visualization.R"

rule ma_plot:
    """Generate MA plot"""
    input:
        results="results/de_analysis/deseq2_results.csv"
    output:
        plot="results/plots/ma_plot.pdf"
    log:
        "logs/ma_plot.log"
    script:
        "scripts/visualization.R"

rule heatmap:
    """Generate heatmap of top DE genes"""
    input:
        counts="results/de_analysis/normalized_counts.csv",
        results="results/de_analysis/deseq2_results.csv"
    output:
        plot="results/plots/heatmap_top_genes.pdf"
    log:
        "logs/heatmap.log"
    script:
        "scripts/visualization.R"

rule pca_plot:
    """Generate PCA plot"""
    input:
        counts="results/de_analysis/normalized_counts.csv",
        metadata=config["metadata_file"]
    output:
        plot="results/plots/pca_plot.pdf"
    log:
        "logs/pca_plot.log"
    script:
        "scripts/visualization.R"

# ============================================================================
# UTILITY RULES
# ============================================================================

rule clean:
    """Remove all results"""
    shell:
        """
        rm -rf results/ logs/
        """

rule clean_all:
    """Remove all results and intermediate files"""
    shell:
        """
        rm -rf results/ logs/ .snakemake/
        """
