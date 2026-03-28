#!/usr/bin/env Rscript

"""
DESeq2 Differential Expression Analysis
Input: featureCounts output
Output: DE results, normalized counts, diagnostic plots
"""

library(DESeq2)
library(tidyverse)
library(ggplot2)

# Load configuration from snakemake
counts_file <- snakemake@input[["counts"]]
metadata_file <- snakemake@input[["metadata"]]
output_results <- snakemake@output[["results"]]
output_normalized <- snakemake@output[["normalized"]]

# Read counts
cat("Reading count data...\n")
counts <- read.table(counts_file, header=TRUE, row.names=1, skip=1)
# Remove unnecessary columns (Chr, Start, End, Strand, Length)
counts <- counts[, -c(1:5)]
# Rename columns to sample IDs
colnames(counts) <- gsub(".*results.alignment.", "", colnames(counts))
colnames(counts) <- gsub(".bam", "", colnames(counts))

# Read metadata
cat("Reading sample metadata...\n")
metadata <- read.csv(metadata_file, row.names=1)
metadata$condition <- factor(metadata$condition)

# Ensure sample order matches
counts <- counts[, rownames(metadata)]

# Create DESeq2 object
cat("Creating DESeq2 object...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = as.formula(snakemake@config[["deseq2"]][["design"]])
)

# Pre-filtering: remove genes with very low counts
cat("Pre-filtering low-count genes...\n")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)

# Extract results
contrast <- unlist(snakemake@config[["deseq2"]][["contrast"]])
res <- results(dds, contrast=contrast, alpha=snakemake@config[["deseq2"]][["padj_threshold"]])

# Shrink log fold changes (LFC shrinkage)
cat("Applying LFC shrinkage...\n")
res <- lfcShrink(dds, contrast=contrast, res=res, type="ashr")

# Convert to data frame and add gene names
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj)

# Save results
cat("Saving results...\n")
write.csv(res_df, output_results, row.names=FALSE)

# Get normalized counts
cat("Extracting normalized counts...\n")
norm_counts <- counts(dds, normalized=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

write.csv(norm_counts, output_normalized, row.names=FALSE)

# Print summary
cat("\n=== DESeq2 Analysis Summary ===\n")
cat("Total genes analyzed:", nrow(dds), "\n")
cat("Genes with padj < 0.05:", sum(res_df$padj < 0.05, na.rm=TRUE), "\n")
cat("Upregulated genes (log2FC > 1):", sum(res_df$log2FoldChange > 1 & res_df$padj < 0.05, na.rm=TRUE), "\n")
cat("Downregulated genes (log2FC < -1):", sum(res_df$log2FoldChange < -1 & res_df$padj < 0.05, na.rm=TRUE), "\n")

cat("\nAnalysis complete!\n")
