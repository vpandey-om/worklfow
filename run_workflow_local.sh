#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON:-python3}"

if [[ -x "/home/vikash/miniconda3/envs/survom-preprocess/bin/python" ]]; then
  PYTHON_BIN="/home/vikash/miniconda3/envs/survom-preprocess/bin/python"
fi

exec "${PYTHON_BIN}" "${ROOT_DIR}/survom-pipelines/scripts/local_workflow_runner.py" "$@"
