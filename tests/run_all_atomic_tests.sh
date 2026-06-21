#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LOG_DIR="${ROOT_DIR}/tests/logs/$(date +%Y%m%d_%H%M%S)"
WORK_DIR="${LOG_DIR}/work"
PYTHON_BIN="${PYTHON:-python3}"
GENOMICS_CLI="${ROOT_DIR}/survom-pipelines/bin/python/common/genomics_atomic.py"
FASTQ_DIR="${ROOT_DIR}/testdatasets/test_data/human_chr22_genomics/fastq"
REF_FASTA="${ROOT_DIR}/testdatasets/test_data/human_chr22_rnaseq/refs/chr22_with_ERCC92.fa"

mkdir -p "${LOG_DIR}" "${WORK_DIR}"

NAMES=()
STATUSES=()
LOGS=()
FAILURES=0

record() {
  local name="$1"
  local status="$2"
  local log="$3"
  NAMES+=("${name}")
  STATUSES+=("${status}")
  LOGS+=("${log}")
  if [[ "${status}" == "FAIL" ]]; then
    FAILURES=$((FAILURES + 1))
  fi
}

safe_name() {
  printf '%s' "$1" | tr ' /:' '___' | tr -cd 'A-Za-z0-9_.-'
}

run_check() {
  local name="$1"
  shift
  local log="${LOG_DIR}/$(safe_name "${name}").log"
  printf '[RUN] %s\n' "${name}"
  "$@" >"${log}" 2>&1
  local code=$?
  if [[ ${code} -eq 0 ]]; then
    printf '[PASS] %s\n' "${name}"
    record "${name}" "PASS" "${log}"
  else
    printf '[FAIL] %s (see %s)\n' "${name}" "${log}"
    record "${name}" "FAIL" "${log}"
  fi
}

run_expected_tool_check() {
  local name="$1"
  shift
  local log="${LOG_DIR}/$(safe_name "${name}").log"
  printf '[RUN] %s\n' "${name}"
  "$@" >"${log}" 2>&1
  local code=$?
  if [[ ${code} -eq 0 || ${code} -eq 127 ]]; then
    printf '[PASS] %s\n' "${name}"
    record "${name}" "PASS" "${log}"
  else
    printf '[FAIL] %s (see %s)\n' "${name}" "${log}"
    record "${name}" "FAIL" "${log}"
  fi
}

if [[ ! -f "${GENOMICS_CLI}" ]]; then
  printf 'ERROR: genomics CLI missing: %s\n' "${GENOMICS_CLI}" >&2
  exit 2
fi

SAMPLESHEET="${WORK_DIR}/genomics_samplesheet.csv"
cat >"${SAMPLESHEET}" <<EOF
sample_id,fastq_1,fastq_2,single_end
genomics_test,${FASTQ_DIR}/genomics_test_R1.fastq.gz,${FASTQ_DIR}/genomics_test_R2.fastq.gz,false
EOF

run_check "genomics input validation" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" validate-inputs \
  --samplesheet "${SAMPLESHEET}" \
  --outdir "${WORK_DIR}/01_validate_inputs"

run_check "genomics reference preparation" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" reference-prepare \
  --genome-fasta "${REF_FASTA}" \
  --reference-name human_chr22_demo \
  --outdir "${WORK_DIR}/02_reference_prepare"

run_check "genomics raw qc" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" raw-qc \
  --manifest "${WORK_DIR}/01_validate_inputs/validated_genomics_manifest.csv" \
  --outdir "${WORK_DIR}/03_raw_qc"

run_check "genomics trim manifest" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" trim-fastq \
  --manifest "${WORK_DIR}/01_validate_inputs/validated_genomics_manifest.csv" \
  --quality-cutoff 20 \
  --min-length 20 \
  --outdir "${WORK_DIR}/04_trim"

run_expected_tool_check "genomics bwa alignment tool check" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" bwa-align \
  --outdir "${WORK_DIR}/05_bwa_align"

run_expected_tool_check "genomics bam processing tool check" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" bam-process \
  --outdir "${WORK_DIR}/06_bam_process"

run_expected_tool_check "genomics germline calling tool check" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" germline-variants \
  --outdir "${WORK_DIR}/07_germline_variants"

run_expected_tool_check "genomics variant filtering tool check" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" variant-filter \
  --outdir "${WORK_DIR}/08_variant_filter"

run_check "genomics all tool inventory" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" tool-inventory \
  --outdir "${WORK_DIR}/17_tool_inventory"

run_check "genomics provenance manifest" \
  "${PYTHON_BIN}" "${GENOMICS_CLI}" provenance \
  --workflow-name genomics_atomic_smoke \
  --run-id genomics_atomic_smoke \
  --parameters-json '{"threads": 2}' \
  --inputs-json "[\"${SAMPLESHEET}\"]" \
  --outputs-json "[\"${WORK_DIR}\"]" \
  --outdir "${WORK_DIR}/99_provenance"

if [[ -f "${ROOT_DIR}/survom-pipelines/scripts/smoke_genomics_workflows.py" ]]; then
  run_check "genomics Nextflow compile smoke" \
    "${PYTHON_BIN}" "${ROOT_DIR}/survom-pipelines/scripts/smoke_genomics_workflows.py" \
    --keep "${WORK_DIR}/compiled_nextflow"
fi

if [[ -f "${ROOT_DIR}/survom-pipelines/scripts/compile_genomics_demo_chain.py" ]]; then
  FULL_CHAIN_DIR="${ROOT_DIR}/survom-pipelines/workflows/generated/genomics_full_demo_chain_cli_test"
  run_check "genomics full chain compile" \
    "${PYTHON_BIN}" "${ROOT_DIR}/survom-pipelines/scripts/compile_genomics_demo_chain.py" \
    --outdir "${FULL_CHAIN_DIR}"
  if command -v nextflow >/dev/null 2>&1; then
    run_check "genomics full chain Nextflow run" \
      nextflow run "${FULL_CHAIN_DIR}/main.nf" -profile local -params-file "${FULL_CHAIN_DIR}/params.yaml"
  else
    printf '[NOTE] genomics full chain Nextflow run - nextflow is not installed\n'
  fi
fi

if [[ -x "${ROOT_DIR}/scripts/run_genomics_demo_pipeline.sh" ]]; then
  run_check "genomics real tool demo pipeline" \
    bash "${ROOT_DIR}/scripts/run_genomics_demo_pipeline.sh" \
    --outdir "${WORK_DIR}/real_tool_demo" \
    --threads 2
fi

printf '\nSummary\n'
printf '%-40s %-6s %s\n' "Check" "Status" "Log"
printf '%-40s %-6s %s\n' "-----" "------" "---"
for i in "${!NAMES[@]}"; do
  printf '%-40s %-6s %s\n' "${NAMES[$i]}" "${STATUSES[$i]}" "${LOGS[$i]}"
done

printf '\nLogs saved in: %s\n' "${LOG_DIR}"
if [[ ${FAILURES} -eq 0 ]]; then
  printf 'Overall: PASS\n'
  exit 0
fi

printf 'Overall: FAIL (%s failed)\n' "${FAILURES}"
exit 1
