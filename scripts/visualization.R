#!/usr/bin/env Rscript

"""
RNA-seq Visualization
Generate publication-quality plots from DESeq2 results
"""

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Read results
results_file <- snakemake@input[["results"]]
res_df <- read.csv(results_file)

# Get rule name to determine which plot to generate
rule_name <- snakemake@rule

if (rule_name == "volcano_plot") {
  
  cat("Generating volcano plot...\n")
  
  # Prepare data
  plot_df <- res_df %>%
    mutate(
      significant = padj < 0.05 & abs(log2FoldChange) > 1,
      color = case_when(
        log2FoldChange > 1 & padj < 0.05 ~ "Up",
        log2FoldChange < -1 & padj < 0.05 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Create volcano plot
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(
      values = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "#95a5a6"),
      name = "Regulation"
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.5) +
    labs(
      title = "Volcano Plot: Differential Expression",
      x = "log2(Fold Change)",
      y = "-log10(Adjusted p-value)",
      subtitle = paste0("Up: ", sum(plot_df$color == "Up"), 
                       " | Down: ", sum(plot_df$color == "Down"))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  ggsave(snakemake@output[[1]], p, width = 10, height = 8, dpi = 300)
  
} else if (rule_name == "ma_plot") {
  
  cat("Generating MA plot...\n")
  
  # Prepare data
  plot_df <- res_df %>%
    mutate(
      color = case_when(
        log2FoldChange > 1 & padj < 0.05 ~ "Up",
        log2FoldChange < -1 & padj < 0.05 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Create MA plot
  p <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = color)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(
      values = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "#95a5a6"),
      name = "Regulation"
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50", alpha = 0.5) +
    labs(
      title = "MA Plot: Differential Expression",
      x = "log10(Mean Expression)",
      y = "log2(Fold Change)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  ggsave(snakemake@output[[1]], p, width = 10, height = 8, dpi = 300)
  
} else if (rule_name == "heatmap") {
  
  cat("Generating heatmap...\n")
  
  # Read normalized counts
  counts_file <- snakemake@input[["counts"]]
  norm_counts <- read.csv(counts_file, row.names=1)
  
  # Get top DE genes
  top_genes <- res_df %>%
    arrange(padj) %>%
    slice_head(n = 20) %>%
    pull(gene_id)
  
  # Prepare heatmap data
  heatmap_data <- norm_counts[top_genes, ] %>%
    as.matrix() %>%
    log2(. + 1)  # Log2 transformation
  
  # Create heatmap
  pdf(snakemake@output[[1]], width = 10, height = 8)
  pheatmap(
    heatmap_data,
    main = "Top 20 Differentially Expressed Genes",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    scale = "row",
    fontsize_row = 8,
    fontsize_col = 10
  )
  dev.off()
  
} else if (rule_name == "pca_plot") {
  
  cat("Generating PCA plot...\n")
  
  # Read normalized counts
  counts_file <- snakemake@input[["counts"]]
  norm_counts <- read.csv(counts_file, row.names=1)
  
  # Read metadata
  metadata_file <- snakemake@input[["metadata"]]
  metadata <- read.csv(metadata_file, row.names=1)
  
  # Get top variable genes
  ntop <- 500
  rv <- rowVars(as.matrix(norm_counts))
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # Perform PCA
  pca_data <- prcomp(t(norm_counts[select, ]), scale=TRUE)
  
  # Prepare plot data
  plot_df <- as.data.frame(pca_data$x[, 1:2]) %>%
    rownames_to_column("sample") %>%
    left_join(metadata %>% rownames_to_column("sample"), by="sample")
  
  # Create PCA plot
  p <- ggplot(plot_df, aes(x=PC1, y=PC2, color=condition, label=sample)) +
    geom_point(size=4, alpha=0.7) +
    geom_text(hjust=1.5, size=3) +
    labs(
      title = "PCA Plot: Sample Clustering",
      x = paste0("PC1 (", round(summary(pca_data)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_data)$importance[2,2]*100, 1), "%)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.position = "right"
    )
  
  ggsave(snakemake@output[[1]], p, width = 10, height = 8, dpi = 300)
}

cat("Plot saved successfully!\n")
