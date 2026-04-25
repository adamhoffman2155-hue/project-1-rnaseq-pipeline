# Benchmarks

This directory holds head-to-head method comparisons — additive to the
canonical pipeline outputs in `results/` proper.

## DESeq2 vs limma-voom

Script: `scripts/benchmark_deseq2_vs_limma.R`

Run on any count matrix + sample sheet:

```
Rscript scripts/benchmark_deseq2_vs_limma.R \
    --counts results/counts.csv \
    --metadata data/metadata.csv \
    --condition-column condition \
    --out results/benchmark/deseq2_vs_limma.csv
```

Outputs a per-gene table with both tools' log2FC and padj plus two
agreement columns (`both_sig`, `concordant_sign`). This is the standard
industry sanity-check when a DE pipeline is first deployed — nf-core/rnaseq's
own benchmark reports show > 85% concordance between pipelines.

Why not just DESeq2? limma-voom handles low-count genes differently and
tends to have tighter FDR control at very small sample sizes, so reporting
both is the honest default.
