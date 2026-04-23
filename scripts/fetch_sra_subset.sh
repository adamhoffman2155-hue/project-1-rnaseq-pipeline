#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# Fetch a small public RNA-seq dataset (6 samples, ~3-5 GB total) so the
# Snakemake DAG can be run end-to-end on real reads, not just the synthetic
# counts fixture.
#
# Default target: GSE109169 — a 6-sample (3 control, 3 treated) TCGA-style
# study. Override via environment variables or edit the accession list.
#
# Usage:
#     scripts/fetch_sra_subset.sh                     # default: 6 samples
#     SRR_ACCS="SRR123 SRR456" scripts/fetch_sra_subset.sh
#
# Requirements:
#   - sra-tools >= 3.0 (prefetch, fasterq-dump) OR curl + direct ENA links
#   - pigz (parallel gzip) for compression
#   - ~20 GB disk free during conversion (SRA → FASTQ → gzip)
#
# Environment:
#   SRR_ACCS         space-separated SRR accessions (default: built-in list)
#   OUT_DIR          destination for FASTQ.gz (default: data/fastq)
#   THREADS          parallel threads for fasterq-dump (default: 4)
#   USE_ENA          1 = download fastq.gz straight from ENA (no sra-tools)
#                    0 = prefetch + fasterq-dump + gzip (slower, more robust)
#
# After this completes, update data/metadata.csv so sample_id values match
# the SRR accessions here, then run the full pipeline:
#
#     snakemake --cores 4
# ----------------------------------------------------------------------------

set -euo pipefail

: "${SRR_ACCS:=SRR6435936 SRR6435937 SRR6435938 SRR6435939 SRR6435940 SRR6435941}"
: "${OUT_DIR:=data/fastq}"
: "${THREADS:=4}"
: "${USE_ENA:=1}"

mkdir -p "$OUT_DIR"

fetch_one_ena() {
    local acc=$1
    local dir6=${acc:0:6}
    local ena_dir="https://ftp.sra.ebi.ac.uk/vol1/fastq/${dir6}/00${acc: -1}/${acc}"
    local r1="${OUT_DIR}/${acc}_R1.fastq.gz"

    if [[ -s "$r1" ]]; then
        echo "[skip] $r1 already present"
        return
    fi
    echo "[ena] $acc -> $r1"
    # Single-end layout assumed (matches the Snakefile's _R1.fastq.gz naming).
    # Paired-end studies would also fetch _2.fastq.gz.
    curl -fsSL "${ena_dir}/${acc}.fastq.gz" -o "$r1"
}

fetch_one_sra() {
    local acc=$1
    local r1="${OUT_DIR}/${acc}_R1.fastq.gz"
    if [[ -s "$r1" ]]; then
        echo "[skip] $r1 already present"
        return
    fi

    command -v prefetch >/dev/null 2>&1 || {
        echo "ERROR: sra-tools 'prefetch' not on PATH. Install via 'conda install -c bioconda sra-tools'." >&2
        exit 1
    }

    echo "[sra] $acc prefetch"
    prefetch --max-size 50g -O "${OUT_DIR}/.sra_tmp" "$acc"
    echo "[sra] $acc fasterq-dump (${THREADS} threads)"
    fasterq-dump --threads "$THREADS" --split-files --skip-technical \
        -O "${OUT_DIR}/.sra_tmp" "$acc"
    # Rename to match Snakefile _R1.fastq.gz convention (single-end)
    if [[ -f "${OUT_DIR}/.sra_tmp/${acc}.fastq" ]]; then
        pigz -p "$THREADS" "${OUT_DIR}/.sra_tmp/${acc}.fastq"
        mv "${OUT_DIR}/.sra_tmp/${acc}.fastq.gz" "$r1"
    elif [[ -f "${OUT_DIR}/.sra_tmp/${acc}_1.fastq" ]]; then
        pigz -p "$THREADS" "${OUT_DIR}/.sra_tmp/${acc}_1.fastq"
        mv "${OUT_DIR}/.sra_tmp/${acc}_1.fastq.gz" "$r1"
    else
        echo "ERROR: fasterq-dump did not produce a FASTQ for $acc" >&2
        exit 1
    fi
}

for acc in $SRR_ACCS; do
    if [[ "$USE_ENA" == "1" ]]; then
        fetch_one_ena "$acc"
    else
        fetch_one_sra "$acc"
    fi
done

echo
echo "Fetched:"
ls -lh "$OUT_DIR"/*.fastq.gz 2>/dev/null | awk '{print "  " $0}'
echo
echo "Next steps:"
echo "  1. Edit data/metadata.csv so sample_id matches the SRR accessions above"
echo "     and condition columns reflect the design."
echo "  2. scripts/download_reference.sh   # GRCh38 FASTA + GTF + STAR index"
echo "  3. snakemake --cores 4             # fastqc -> STAR -> counts -> DE -> GSEA"
