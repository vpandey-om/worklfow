#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${1:-test_data/human_chr22_rnaseq}"

FASTQ_DIR="${BASE_DIR}/fastq"
REF_DIR="${BASE_DIR}/refs"
OUT_DIR="${BASE_DIR}/output"

mkdir -p "${FASTQ_DIR}" "${REF_DIR}" "${OUT_DIR}"

download_file() {
  local url="$1"
  local output="$2"

  if [[ -s "${output}" ]]; then
    echo "Already exists: ${output}"
    return 0
  fi

  if command -v wget >/dev/null 2>&1; then
    wget -c "${url}" -O "${output}"
  elif command -v curl >/dev/null 2>&1; then
    curl -L "${url}" -o "${output}"
  else
    echo "ERROR: wget or curl is required"
    exit 1
  fi
}

echo "========================================"
echo "Downloading human chr22 RNA-seq data"
echo "Base folder: ${BASE_DIR}"
echo "========================================"

echo
echo "1) Downloading RNA-seq FASTQ archive..."
download_file \
  "http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar" \
  "${FASTQ_DIR}/HBR_UHR_ERCC_ds_5pc.tar"

echo
echo "2) Extracting FASTQ files..."
tar -xf "${FASTQ_DIR}/HBR_UHR_ERCC_ds_5pc.tar" -C "${FASTQ_DIR}"

echo
echo "3) Downloading chr22 + ERCC reference FASTA..."
download_file \
  "http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa" \
  "${REF_DIR}/chr22_with_ERCC92.fa"

echo
echo "4) Downloading chr22 + ERCC annotation GTF..."
download_file \
  "http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf" \
  "${REF_DIR}/chr22_with_ERCC92.gtf"

echo
echo "5) Checking files..."

echo
echo "FASTQ files:"
find "${FASTQ_DIR}" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) -ls

echo
echo "Reference files:"
ls -lh "${REF_DIR}"

FASTQ_COUNT="$(find "${FASTQ_DIR}" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) | wc -l)"

echo
echo "FASTQ file count: ${FASTQ_COUNT}"

if [[ "${FASTQ_COUNT}" -eq 0 ]]; then
  echo "ERROR: No FASTQ files found after extraction"
  exit 1
fi

if [[ ! -s "${REF_DIR}/chr22_with_ERCC92.fa" ]]; then
  echo "ERROR: Missing reference FASTA"
  exit 1
fi

if [[ ! -s "${REF_DIR}/chr22_with_ERCC92.gtf" ]]; then
  echo "ERROR: Missing annotation GTF"
  exit 1
fi

echo
echo "========================================"
echo "DONE"
echo "FASTQ dir:"
echo "  ${FASTQ_DIR}"
echo
echo "Reference FASTA:"
echo "  ${REF_DIR}/chr22_with_ERCC92.fa"
echo
echo "Annotation GTF:"
echo "  ${REF_DIR}/chr22_with_ERCC92.gtf"
echo
echo "Output dir:"
echo "  ${OUT_DIR}"
echo "========================================"