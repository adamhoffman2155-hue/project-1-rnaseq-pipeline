#!/usr/bin/env Rscript
# Benchmark: DESeq2 vs limma-voom on the same count matrix.
#
# Additive-only benchmark script. Reads a user-supplied gene-level raw count
# matrix and two-group sample sheet, runs both DESeq2 and limma-voom on the
# same data, and writes a side-by-side concordance table to
# results/benchmark/deseq2_vs_limma.csv. Does not modify any POC output.
#
# Usage:
#   Rscript scripts/benchmark_deseq2_vs_limma.R \
#       --counts results/counts.csv \
#       --metadata data/metadata.csv \
#       --condition-column condition \
#       --out results/benchmark/deseq2_vs_limma.csv
#
# CI runs this only on `workflow_dispatch` (see .github/workflows/ci.yml)
# because the full R toolchain is heavy to install on every push.

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(limma)
  library(edgeR)
})

option_list <- list(
  make_option("--counts", type = "character", help = "Path to raw counts CSV (gene x sample, gene IDs in first column)"),
  make_option("--metadata", type = "character", help = "Path to sample metadata CSV with at least 'sample_id' and condition column"),
  make_option("--condition-column", type = "character", default = "condition", help = "Metadata column used as the two-group factor (default: condition)"),
  make_option("--out", type = "character", default = "results/benchmark/deseq2_vs_limma.csv", help = "Output CSV path")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$counts), !is.null(opt$metadata))
counts <- read.csv(opt$counts, row.names = 1, check.names = FALSE)
meta   <- read.csv(opt$metadata)

meta <- meta[meta$sample_id %in% colnames(counts), , drop = FALSE]
counts <- counts[, meta$sample_id]
cond <- factor(meta[[opt[["condition-column"]]]])
stopifnot(length(levels(cond)) == 2)

# DESeq2 ------------------------------------------------------------------
message("Running DESeq2 ...")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds)
deseq2_df <- data.frame(
  gene = rownames(res),
  deseq2_log2fc = res$log2FoldChange,
  deseq2_padj = res$padj
)

# limma-voom --------------------------------------------------------------
message("Running limma-voom ...")
design <- model.matrix(~ cond)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
limma_tbl <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
limma_df <- data.frame(
  gene = rownames(limma_tbl),
  limma_log2fc = limma_tbl$logFC,
  limma_padj = limma_tbl$adj.P.Val
)

# Merge + concordance -----------------------------------------------------
merged <- merge(deseq2_df, limma_df, by = "gene", all = FALSE)
merged$both_sig <- merged$deseq2_padj < 0.05 & merged$limma_padj < 0.05
merged$concordant_sign <- sign(merged$deseq2_log2fc) == sign(merged$limma_log2fc)

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
write.csv(merged, opt$out, row.names = FALSE)

both <- sum(merged$both_sig, na.rm = TRUE)
d2   <- sum(merged$deseq2_padj < 0.05, na.rm = TRUE)
lv   <- sum(merged$limma_padj  < 0.05, na.rm = TRUE)
sign_agreement <- mean(merged$concordant_sign, na.rm = TRUE)

message(sprintf(
  "DESeq2 sig: %d | limma sig: %d | both: %d | sign agreement: %.2f%%",
  d2, lv, both, 100 * sign_agreement
))
message(sprintf("wrote %s", opt$out))
