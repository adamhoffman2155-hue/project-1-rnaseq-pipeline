#!/usr/bin/env bash
# =============================================================================
# Download GRCh38 Reference Genome and Build STAR Index
#
# Usage:
#   bash scripts/download_reference.sh [output_dir] [threads]
#
# Arguments:
#   output_dir  Directory to store reference files (default: reference)
#   threads     Number of threads for STAR index build (default: 8)
#
# Requirements:
#   - wget, md5sum
#   - STAR (>= 2.7) on PATH
#
# Downloads GENCODE v38 primary assembly FASTA and GTF annotation,
# validates checksums, and builds a STAR genome index.
# =============================================================================
set -euo pipefail

OUTPUT_DIR="${1:-reference}"
THREADS="${2:-8}"

GENCODE_VERSION="38"
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}"

FASTA_FILE="GRCh38.primary_assembly.genome.fa"
GTF_FILE="gencode.v${GENCODE_VERSION}.primary_assembly.annotation.gtf"

FASTA_URL="${BASE_URL}/${FASTA_FILE}.gz"
GTF_URL="${BASE_URL}/${GTF_FILE}.gz"
MD5_URL="${BASE_URL}/MD5SUMS"

STAR_INDEX_DIR="${OUTPUT_DIR}/star_index"

# ---- Helper functions ------------------------------------------------------
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

check_tool() {
  command -v "$1" >/dev/null 2>&1 || die "Required tool not found: $1"
}

# ---- Pre-flight checks -----------------------------------------------------
check_tool wget
check_tool md5sum
check_tool STAR

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${STAR_INDEX_DIR}"

# ---- Download MD5SUMS ------------------------------------------------------
log "Downloading MD5SUMS from GENCODE v${GENCODE_VERSION}"
wget -q -O "${OUTPUT_DIR}/MD5SUMS" "${MD5_URL}"

# ---- Download FASTA --------------------------------------------------------
if [[ -f "${OUTPUT_DIR}/${FASTA_FILE}" ]]; then
  log "FASTA already exists, skipping download"
else
  log "Downloading primary assembly FASTA: ${FASTA_FILE}.gz"
  wget -q --show-progress -O "${OUTPUT_DIR}/${FASTA_FILE}.gz" "${FASTA_URL}"

  log "Validating FASTA checksum"
  expected_md5=$(grep "${FASTA_FILE}.gz" "${OUTPUT_DIR}/MD5SUMS" | awk '{print $1}')
  actual_md5=$(md5sum "${OUTPUT_DIR}/${FASTA_FILE}.gz" | awk '{print $1}')
  if [[ "${expected_md5}" != "${actual_md5}" ]]; then
    die "FASTA checksum mismatch: expected ${expected_md5}, got ${actual_md5}"
  fi
  log "FASTA checksum OK"

  log "Decompressing FASTA"
  gunzip "${OUTPUT_DIR}/${FASTA_FILE}.gz"
fi

# ---- Download GTF ----------------------------------------------------------
if [[ -f "${OUTPUT_DIR}/${GTF_FILE}" ]]; then
  log "GTF already exists, skipping download"
else
  log "Downloading GTF annotation: ${GTF_FILE}.gz"
  wget -q --show-progress -O "${OUTPUT_DIR}/${GTF_FILE}.gz" "${GTF_URL}"

  log "Validating GTF checksum"
  expected_md5=$(grep "${GTF_FILE}.gz" "${OUTPUT_DIR}/MD5SUMS" | awk '{print $1}')
  actual_md5=$(md5sum "${OUTPUT_DIR}/${GTF_FILE}.gz" | awk '{print $1}')
  if [[ "${expected_md5}" != "${actual_md5}" ]]; then
    die "GTF checksum mismatch: expected ${expected_md5}, got ${actual_md5}"
  fi
  log "GTF checksum OK"

  log "Decompressing GTF"
  gunzip "${OUTPUT_DIR}/${GTF_FILE}.gz"
fi

# ---- Build STAR index ------------------------------------------------------
if [[ -f "${STAR_INDEX_DIR}/SA" ]]; then
  log "STAR index already exists, skipping build"
else
  log "Building STAR genome index with ${THREADS} threads (this may take 30+ minutes)"
  STAR \
    --runMode genomeGenerate \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX_DIR}" \
    --genomeFastaFiles "${OUTPUT_DIR}/${FASTA_FILE}" \
    --sjdbGTFfile "${OUTPUT_DIR}/${GTF_FILE}" \
    --sjdbOverhang 100

  log "STAR index build complete"
fi

# ---- Summary ---------------------------------------------------------------
log "Reference setup complete."
log "  FASTA : ${OUTPUT_DIR}/${FASTA_FILE}"
log "  GTF   : ${OUTPUT_DIR}/${GTF_FILE}"
log "  Index : ${STAR_INDEX_DIR}/"
