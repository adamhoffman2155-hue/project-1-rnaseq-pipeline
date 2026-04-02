#!/usr/bin/env Rscript
# =============================================================================
# GSEA / Pathway Enrichment Analysis
# Uses fgsea package to perform Gene Set Enrichment Analysis on DESeq2 results.
#
# Input : DESeq2 differential expression results CSV (gene, log2FoldChange, pvalue, padj)
# Output: enrichment_results.csv, gsea_barplot.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# ---- Configuration ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript gsea_analysis.R <deseq2_results.csv> <output_dir>")
}

deseq2_file <- args[1]
output_dir  <- args[2]

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- Load DESeq2 results ---------------------------------------------------
cat("Loading DESeq2 results from:", deseq2_file, "\n")
res <- fread(deseq2_file)

# Expect columns: gene (or first column), log2FoldChange, pvalue
if (!"gene" %in% colnames(res)) {
  colnames(res)[1] <- "gene"
}

stopifnot(all(c("gene", "log2FoldChange", "pvalue") %in% colnames(res)))

# Remove rows with NA pvalue or log2FoldChange
res <- res[!is.na(pvalue) & !is.na(log2FoldChange)]

# ---- Build ranked gene list ------------------------------------------------
# Rank metric: sign-aware combination of fold change and significance
res[, rank_metric := log2FoldChange * (-log10(pvalue))]
ranked_genes <- setNames(res$rank_metric, res$gene)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

cat("Ranked gene list size:", length(ranked_genes), "\n")

# ---- Load gene sets from MSigDB --------------------------------------------
# Attempt to load GMT files; fall back to msigdbr if available
gmt_hallmark <- file.path("data", "h.all.v2023.2.Hs.symbols.gmt")
gmt_kegg     <- file.path("data", "c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")

if (file.exists(gmt_hallmark) && file.exists(gmt_kegg)) {
  cat("Loading gene sets from local GMT files\n")
  pathways_hallmark <- gmtPathways(gmt_hallmark)
  pathways_kegg     <- gmtPathways(gmt_kegg)
  pathways <- c(pathways_hallmark, pathways_kegg)
} else if (requireNamespace("msigdbr", quietly = TRUE)) {
  cat("Loading gene sets via msigdbr\n")
  library(msigdbr)
  hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
  kegg_df     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
  combined_df <- bind_rows(hallmark_df, kegg_df)
  pathways <- split(combined_df$gene_symbol, combined_df$gs_name)
} else {
  stop("No gene set source found. Provide GMT files in data/ or install msigdbr.")
}

cat("Total gene sets loaded:", length(pathways), "\n")

# ---- Run fgsea -------------------------------------------------------------
cat("Running fgsea ...\n")
set.seed(42)
fgsea_res <- fgsea(
  pathways  = pathways,
  stats     = ranked_genes,
  minSize   = 15,
  maxSize   = 500,
  nPermSimple = 10000
)

fgsea_res <- fgsea_res[order(pval)]

cat("Significant pathways (padj < 0.05):", sum(fgsea_res$padj < 0.05, na.rm = TRUE), "\n")

# ---- Highlight DDR and mismatch repair pathways ----------------------------
ddr_keywords <- c("DNA_REPAIR", "DNA_DAMAGE", "MISMATCH_REPAIR",
                  "BASE_EXCISION", "NUCLEOTIDE_EXCISION", "HOMOLOGOUS_RECOMBINATION",
                  "P53_PATHWAY", "G2M_CHECKPOINT", "E2F_TARGETS")

fgsea_res[, ddr_related := any(sapply(ddr_keywords, function(k) grepl(k, pathway, ignore.case = TRUE))), by = pathway]

ddr_results <- fgsea_res[ddr_related == TRUE]
if (nrow(ddr_results) > 0) {
  cat("\n--- DDR / Mismatch Repair Pathway Results ---\n")
  print(ddr_results[, .(pathway, pval, padj, NES, size)])
} else {
  cat("\nNo DDR / mismatch repair pathways found in results.\n")
}

# ---- Write results ---------------------------------------------------------
out_table <- fgsea_res[, .(pathway, pval, padj, ES, NES, size, ddr_related)]
out_csv <- file.path(output_dir, "enrichment_results.csv")
fwrite(out_table, out_csv)
cat("Enrichment table written to:", out_csv, "\n")

# ---- Barplot of top pathways -----------------------------------------------
top_n <- 20
plot_data <- head(fgsea_res[order(padj)], top_n)
plot_data[, pathway_short := gsub("^HALLMARK_|^KEGG_", "", pathway)]
plot_data[, pathway_short := gsub("_", " ", pathway_short)]
plot_data[, direction := ifelse(NES > 0, "Up", "Down")]

p <- ggplot(plot_data, aes(x = reorder(pathway_short, NES), y = NES, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "#d73027", "Down" = "#4575b4")) +
  labs(
    title = "Top Enriched Pathways (GSEA)",
    x = NULL,
    y = "Normalized Enrichment Score (NES)",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))

out_plot <- file.path(output_dir, "gsea_barplot.pdf")
ggsave(out_plot, p, width = 10, height = 7)
cat("Barplot saved to:", out_plot, "\n")

cat("\nGSEA analysis complete.\n")
