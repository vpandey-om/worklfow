#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENOMICS_BIN="${SURVOM_GENOMICS_BIN:-/home/vikash/miniconda3/envs/survom-genomics/bin}"
PYTHON_BIN="${PYTHON:-/home/vikash/miniconda3/envs/survom-preprocess/bin/python}"
GENOMICS_CLI="${ROOT_DIR}/survom-pipelines/bin/python/common/genomics_atomic.py"

SAMPLESHEET=""
REFERENCE="${ROOT_DIR}/testdatasets/test_data/human_chr22_rnaseq/refs/chr22_with_ERCC92.fa"
OUTDIR="${ROOT_DIR}/.test-runs/genomics_demo_$(date +%Y%m%d_%H%M%S)"
THREADS=2

usage() {
  cat <<'USAGE'
Run the SurvOm Genomics demo pipeline.

Usage:
  scripts/run_genomics_demo_pipeline.sh [--samplesheet samples.csv] [--reference ref.fa] [--outdir output_dir] [--threads 2]

Sample sheet columns:
  sample_id,fastq_1,fastq_2,single_end

If --samplesheet is omitted, the human chr22 demo FASTQs are used.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

require_tool() {
  local tool="$1"
  if [[ -x "${GENOMICS_BIN}/${tool}" ]]; then
    printf '%s' "${GENOMICS_BIN}/${tool}"
    return 0
  fi
  if command -v "${tool}" >/dev/null 2>&1; then
    command -v "${tool}"
    return 0
  fi
  echo "Missing required tool: ${tool}. Set SURVOM_GENOMICS_BIN or activate survom-genomics." >&2
  return 1
}

mkdir -p "${OUTDIR}"
LOG_DIR="${OUTDIR}/logs"
mkdir -p "${LOG_DIR}"

if [[ -z "${SAMPLESHEET}" ]]; then
  FASTQ_DIR="${ROOT_DIR}/testdatasets/test_data/human_chr22_genomics/fastq"
  SAMPLESHEET="${OUTDIR}/genomics_demo_samplesheet.csv"
  cat >"${SAMPLESHEET}" <<EOF
sample_id,fastq_1,fastq_2,single_end
genomics_test,${FASTQ_DIR}/genomics_test_R1.fastq.gz,${FASTQ_DIR}/genomics_test_R2.fastq.gz,false
EOF
fi

if [[ ! -s "${SAMPLESHEET}" ]]; then
  echo "Sample sheet not found: ${SAMPLESHEET}" >&2
  exit 2
fi
if [[ ! -s "${REFERENCE}" ]]; then
  echo "Reference FASTA not found: ${REFERENCE}" >&2
  exit 2
fi

BWA="$(require_tool bwa-mem2 || require_tool bwa)"
SAMTOOLS="$(require_tool samtools)"
BCFTOOLS="$(require_tool bcftools)"
MOSDEPTH="$(require_tool mosdepth)"
MULTIQC="$(require_tool multiqc)"

echo "Output: ${OUTDIR}"
echo "Sample sheet: ${SAMPLESHEET}"
echo "Reference: ${REFERENCE}"
echo "Genomics tools: ${GENOMICS_BIN}"

"${PYTHON_BIN}" "${GENOMICS_CLI}" validate-inputs \
  --samplesheet "${SAMPLESHEET}" \
  --outdir "${OUTDIR}/01_input_validation" >"${LOG_DIR}/01_input_validation.log" 2>&1

"${PYTHON_BIN}" "${GENOMICS_CLI}" reference-prepare \
  --genome-fasta "${REFERENCE}" \
  --reference-name human_chr22_demo \
  --outdir "${OUTDIR}/02_reference_prepare" >"${LOG_DIR}/02_reference_prepare.log" 2>&1

REF_WORK="${OUTDIR}/02_reference_prepare/reference.fa"
cp "${REFERENCE}" "${REF_WORK}"
"${SAMTOOLS}" faidx "${REF_WORK}" >"${LOG_DIR}/02_samtools_faidx.log" 2>&1
"${BWA}" index "${REF_WORK}" >"${LOG_DIR}/02_bwa_index.log" 2>&1

"${PYTHON_BIN}" "${GENOMICS_CLI}" raw-qc \
  --manifest "${OUTDIR}/01_input_validation/validated_genomics_manifest.csv" \
  --outdir "${OUTDIR}/03_raw_qc" >"${LOG_DIR}/03_raw_qc.log" 2>&1

"${PYTHON_BIN}" "${GENOMICS_CLI}" trim-fastq \
  --manifest "${OUTDIR}/01_input_validation/validated_genomics_manifest.csv" \
  --quality-cutoff 20 \
  --min-length 20 \
  --outdir "${OUTDIR}/04_trim" >"${LOG_DIR}/04_trim.log" 2>&1

mkdir -p "${OUTDIR}/05_bwa_align" "${OUTDIR}/06_bam_processing" "${OUTDIR}/07_mark_duplicates" "${OUTDIR}/08_coverage_qc"
mkdir -p "${OUTDIR}/09_bqsr" "${OUTDIR}/10_germline_variants" "${OUTDIR}/11_somatic_variants" "${OUTDIR}/12_variant_filtering"
mkdir -p "${OUTDIR}/13_structural_variants" "${OUTDIR}/14_cnv_calling" "${OUTDIR}/15_variant_annotation" "${OUTDIR}/16_multiqc"

tail -n +2 "${OUTDIR}/04_trim/trimmed_fastq_manifest.csv" | while IFS=, read -r sample_id fastq_1 fastq_2 single_end strandedness; do
  if [[ "${fastq_1}" = /* ]]; then
    r1="${fastq_1}"
  else
    r1="${OUTDIR}/04_trim/${fastq_1}"
  fi
  if [[ "${fastq_2}" = /* ]]; then
    r2="${fastq_2}"
  else
    r2="${OUTDIR}/04_trim/${fastq_2}"
  fi
  sam="${OUTDIR}/05_bwa_align/${sample_id}.sam"
  bam="${OUTDIR}/06_bam_processing/${sample_id}.sorted.bam"
  marked="${OUTDIR}/07_mark_duplicates/${sample_id}.marked.bam"
  raw_vcf="${OUTDIR}/10_germline_variants/${sample_id}.raw.vcf.gz"
  filtered_vcf="${OUTDIR}/12_variant_filtering/${sample_id}.filtered.vcf.gz"

  if [[ "${single_end}" == "true" || -z "${fastq_2}" ]]; then
    "${BWA}" mem -t "${THREADS}" -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA" "${REF_WORK}" "${r1}" >"${sam}" 2>"${LOG_DIR}/05_${sample_id}_bwa.log"
  else
    "${BWA}" mem -t "${THREADS}" -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA" "${REF_WORK}" "${r1}" "${r2}" >"${sam}" 2>"${LOG_DIR}/05_${sample_id}_bwa.log"
  fi
  "${SAMTOOLS}" sort -@ "${THREADS}" -o "${bam}" "${sam}" >"${LOG_DIR}/06_${sample_id}_sort.log" 2>&1
  "${SAMTOOLS}" index "${bam}" >"${LOG_DIR}/06_${sample_id}_index.log" 2>&1
  "${SAMTOOLS}" flagstat "${bam}" >"${OUTDIR}/06_bam_processing/${sample_id}.flagstat.txt" 2>"${LOG_DIR}/06_${sample_id}_flagstat.log"
  "${SAMTOOLS}" stats "${bam}" >"${OUTDIR}/06_bam_processing/${sample_id}.stats.txt" 2>"${LOG_DIR}/06_${sample_id}_stats.log"

  cp "${bam}" "${marked}"
  cp "${bam}.bai" "${marked}.bai"
  printf '{"sample_id":"%s","status":"copied_no_duplicate_removal","reason":"Demo runner keeps duplicates; use samtools markdup/Picard in production."}\n' "${sample_id}" >"${OUTDIR}/07_mark_duplicates/${sample_id}.markdup.json"

  "${MOSDEPTH}" -t "${THREADS}" "${OUTDIR}/08_coverage_qc/${sample_id}" "${marked}" >"${LOG_DIR}/08_${sample_id}_mosdepth.log" 2>&1 || \
    printf '{"sample_id":"%s","status":"mosdepth_failed","see_log":"%s"}\n' "${sample_id}" "${LOG_DIR}/08_${sample_id}_mosdepth.log" >"${OUTDIR}/08_coverage_qc/${sample_id}.coverage_status.json"

  printf '{"sample_id":"%s","status":"skipped_no_known_sites","reason":"BQSR requires known-sites VCFs."}\n' "${sample_id}" >"${OUTDIR}/09_bqsr/${sample_id}.bqsr_status.json"

  "${BCFTOOLS}" mpileup -Ou -f "${REF_WORK}" "${marked}" 2>"${LOG_DIR}/10_${sample_id}_mpileup.log" | \
    "${BCFTOOLS}" call -mv -Oz -o "${raw_vcf}" 2>"${LOG_DIR}/10_${sample_id}_bcftools_call.log" || true
  if [[ -s "${raw_vcf}" ]]; then
    "${BCFTOOLS}" index -f "${raw_vcf}" >"${LOG_DIR}/10_${sample_id}_bcftools_index.log" 2>&1 || true
    "${BCFTOOLS}" filter -Oz -o "${filtered_vcf}" -i 'QUAL>=20' "${raw_vcf}" >"${LOG_DIR}/12_${sample_id}_bcftools_filter.log" 2>&1 || true
    [[ -s "${filtered_vcf}" ]] && "${BCFTOOLS}" index -f "${filtered_vcf}" >"${LOG_DIR}/12_${sample_id}_bcftools_index.log" 2>&1 || true
  fi

  printf '{"sample_id":"%s","status":"skipped_no_tumor_normal_design","reason":"Somatic calling requires tumor/normal metadata and resources."}\n' "${sample_id}" >"${OUTDIR}/11_somatic_variants/${sample_id}.somatic_status.json"
  printf '{"sample_id":"%s","status":"skipped_in_demo","reason":"SV calling is tool-ready, but demo FASTQs are too small for meaningful SV calls."}\n' "${sample_id}" >"${OUTDIR}/13_structural_variants/${sample_id}.sv_status.json"
  printf '{"sample_id":"%s","status":"skipped_no_targets_or_access_reference","reason":"CNVkit needs target/access/reference configuration."}\n' "${sample_id}" >"${OUTDIR}/14_cnv_calling/${sample_id}.cnv_status.json"
  printf '{"sample_id":"%s","status":"skipped_no_snpeff_database","reason":"Annotation requires an installed SnpEff/VEP genome database matching the reference."}\n' "${sample_id}" >"${OUTDIR}/15_variant_annotation/${sample_id}.annotation_status.json"
done

"${MULTIQC}" "${OUTDIR}" -o "${OUTDIR}/16_multiqc" >"${LOG_DIR}/16_multiqc.log" 2>&1 || true

"${PYTHON_BIN}" "${GENOMICS_CLI}" provenance \
  --workflow-name genomics_demo_pipeline \
  --run-id "$(basename "${OUTDIR}")" \
  --parameters-json "{\"threads\": ${THREADS}}" \
  --inputs-json "[\"${SAMPLESHEET}\", \"${REFERENCE}\"]" \
  --outputs-json "[\"${OUTDIR}\"]" \
  --outdir "${OUTDIR}/99_provenance" >"${LOG_DIR}/99_provenance.log" 2>&1

cat >"${OUTDIR}/README_outputs.txt" <<EOF
Key outputs:
- 01_input_validation/validated_genomics_manifest.csv
- 02_reference_prepare/reference_bundle.json
- 03_raw_qc/raw_qc_metrics.csv
- 04_trim/trimmed_fastq_manifest.csv
- 06_bam_processing/*.sorted.bam and *.bai
- 06_bam_processing/*.flagstat.txt
- 08_coverage_qc/*.mosdepth.summary.txt if mosdepth completed
- 10_germline_variants/*.raw.vcf.gz if variants were found
- 12_variant_filtering/*.filtered.vcf.gz if variants passed filters
- 16_multiqc/multiqc_report.html if MultiQC completed
- 99_provenance/run_manifest.json
EOF

echo "Overall: PASS"
echo "Results: ${OUTDIR}"
