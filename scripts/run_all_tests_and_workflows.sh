#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="${ROOT_DIR}/.test-runs/${RUN_ID}"
mkdir -p "${LOG_DIR}"
RUN_REAL_WORKFLOWS=1

usage() {
  cat <<'USAGE'
Usage:
  ./scripts/run_all_tests_and_workflows.sh [--quick]

Default:
  Runs code checks, unit tests, workflow compile smoke tests, and real tiny
  Nextflow workflow execution tests using generated demo FASTQ/reference inputs.

Options:
  --quick        Run only fast code/unit/compile checks. Skip real Nextflow execution.
  -h, --help    Show this help message.
USAGE
}

for arg in "$@"; do
  case "${arg}" in
    --quick)
      RUN_REAL_WORKFLOWS=0
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown option: %s\n\n' "${arg}" >&2
      usage >&2
      exit 2
      ;;
  esac
done

NAMES=()
STATUSES=()
LOGS=()
FAILURES=0

record_result() {
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
  local workdir="$2"
  shift 2
  local log="${LOG_DIR}/$(safe_name "${name}").log"

  printf '[RUN] %s\n' "${name}"
  (
    cd "${workdir}" || exit 1
    "$@"
  ) >"${log}" 2>&1
  local code=$?

  if [[ ${code} -eq 0 ]]; then
    printf '[PASS] %s\n' "${name}"
    record_result "${name}" "PASS" "${log}"
  else
    printf '[FAIL] %s (see %s)\n' "${name}" "${log}"
    record_result "${name}" "FAIL" "${log}"
  fi
}

skip_check() {
  local name="$1"
  local reason="$2"
  local log="${LOG_DIR}/$(safe_name "${name}").log"
  printf '%s\n' "${reason}" >"${log}"
  printf '[NOTE] %s - %s\n' "${name}" "${reason}"
}

PYTHON_BIN="${PYTHON:-python3}"
if [[ -x "/home/vikash/miniconda3/envs/survom-preprocess/bin/python" ]]; then
  PYTHON_BIN="/home/vikash/miniconda3/envs/survom-preprocess/bin/python"
fi

run_check "python compile checks" "${ROOT_DIR}" "${PYTHON_BIN}" -m compileall -q \
  demo_dash_app survom-pipelines/builder survom-pipelines/scripts

if [[ -f "${ROOT_DIR}/survom-pipelines/tests/test_local_workflow_runner.py" ]]; then
  run_check "python unittest local runner" "${ROOT_DIR}" \
    "${PYTHON_BIN}" -m unittest discover -s survom-pipelines/tests -p test_local_workflow_runner.py
fi

if [[ -d "${ROOT_DIR}/survom-pipelines/tests" ]]; then
  if "${PYTHON_BIN}" -c "import pytest" >/dev/null 2>&1; then
    run_check "pytest survom-pipelines" "${ROOT_DIR}/survom-pipelines" "${PYTHON_BIN}" -m pytest tests
  else
    skip_check "pytest survom-pipelines" "pytest is not installed for ${PYTHON_BIN}"
  fi
fi

if [[ -f "${ROOT_DIR}/survom-pipelines/scripts/smoke_workflows.py" ]]; then
  run_check "workflow smoke compile" "${ROOT_DIR}/survom-pipelines" \
    "${PYTHON_BIN}" scripts/smoke_workflows.py
  if [[ "${RUN_REAL_WORKFLOWS}" -eq 1 ]]; then
    run_check "workflow real Nextflow smoke" "${ROOT_DIR}/survom-pipelines" \
      "${PYTHON_BIN}" scripts/smoke_workflows.py --run --keep "${LOG_DIR}/workflow_real_nextflow"
  else
    skip_check "workflow real Nextflow smoke" "skipped because --quick was used"
  fi
fi

if [[ -f "${ROOT_DIR}/Makefile" ]]; then
  run_check "make test" "${ROOT_DIR}" make test
fi

while IFS= read -r package_json; do
  package_dir="$(dirname "${package_json}")"
  if command -v npm >/dev/null 2>&1; then
    run_check "npm test ${package_dir#${ROOT_DIR}/}" "${package_dir}" npm test
  else
    skip_check "npm test ${package_dir#${ROOT_DIR}/}" "npm is not installed"
  fi
done < <(find "${ROOT_DIR}" -maxdepth 3 -name package.json -not -path '*/node_modules/*')

while IFS= read -r go_mod; do
  go_dir="$(dirname "${go_mod}")"
  if command -v go >/dev/null 2>&1; then
    run_check "go test ${go_dir#${ROOT_DIR}/}" "${go_dir}" go test ./...
  else
    skip_check "go test ${go_dir#${ROOT_DIR}/}" "go is not installed"
  fi
done < <(find "${ROOT_DIR}" -maxdepth 3 -name go.mod)

while IFS= read -r cargo_toml; do
  cargo_dir="$(dirname "${cargo_toml}")"
  if command -v cargo >/dev/null 2>&1; then
    run_check "cargo test ${cargo_dir#${ROOT_DIR}/}" "${cargo_dir}" cargo test
  else
    skip_check "cargo test ${cargo_dir#${ROOT_DIR}/}" "cargo is not installed"
  fi
done < <(find "${ROOT_DIR}" -maxdepth 3 -name Cargo.toml)

printf '\nSummary\n'
printf '%-34s %-6s %s\n' "Check" "Status" "Log"
printf '%-34s %-6s %s\n' "-----" "------" "---"
for i in "${!NAMES[@]}"; do
  printf '%-34s %-6s %s\n' "${NAMES[$i]}" "${STATUSES[$i]}" "${LOGS[$i]}"
done

printf '\nLogs saved in: %s\n' "${LOG_DIR}"
if [[ ${FAILURES} -eq 0 ]]; then
  printf 'Overall: PASS\n'
  exit 0
fi

printf 'Overall: FAIL (%s failed)\n' "${FAILURES}"
exit 1
