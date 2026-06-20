#!/usr/bin/env bash
set -u

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON:-python3}"
if [[ -x "/home/vikash/miniconda3/envs/survom-preprocess/bin/python" ]]; then
  PYTHON_BIN="/home/vikash/miniconda3/envs/survom-preprocess/bin/python"
fi

AIRWAY_DIR="${ROOT_DIR}/testdatasets/countdata/test_data/airway"
RUN_ID="airway_downstream_$(date +%Y%m%d_%H%M%S)"
OUT_ROOT="${ROOT_DIR}/.test-runs/${RUN_ID}"
COUNTS=""
METADATA=""
CONTRASTS=""

usage() {
  cat <<'USAGE'
Usage:
  ./scripts/run_airway_downstream_atomic.sh [--airway-dir PATH] [--outdir PATH]
  ./scripts/run_airway_downstream_atomic.sh --counts counts.csv --metadata metadata.csv --contrasts contrasts.csv --outdir PATH

Runs all currently implemented downstream RNA-seq atomic CLI steps on the
airway count-matrix dataset:
  merge, validate, filter, normalize, PCA, UMAP, differential expression,
  PLS-DA, GO enrichment, pathway enrichment, ranked GSEA, status/state.

Outputs and logs are written under .test-runs/airway_downstream_<timestamp>/ by default.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --airway-dir)
      AIRWAY_DIR="$2"
      shift 2
      ;;
    --outdir)
      OUT_ROOT="$2"
      shift 2
      ;;
    --counts)
      COUNTS="$2"
      shift 2
      ;;
    --metadata)
      METADATA="$2"
      shift 2
      ;;
    --contrasts)
      CONTRASTS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

CLI="${ROOT_DIR}/survom-pipelines/downstream/bin/downstream_cli.py"
COUNTS="${COUNTS:-${AIRWAY_DIR}/airway_raw_counts.csv}"
METADATA="${METADATA:-${AIRWAY_DIR}/airway_sample_metadata.csv}"
CONTRASTS="${CONTRASTS:-${AIRWAY_DIR}/airway_contrasts.csv}"
INPUTS_DIR="${OUT_ROOT}/inputs"
LOG_DIR="${OUT_ROOT}/logs"
mkdir -p "${INPUTS_DIR}" "${LOG_DIR}"

NAMES=()
STATUSES=()
LOGS=()
FAILURES=0

record() {
  NAMES+=("$1")
  STATUSES+=("$2")
  LOGS+=("$3")
  if [[ "$2" == "FAIL" ]]; then
    FAILURES=$((FAILURES + 1))
  fi
}

run_step() {
  local name="$1"
  shift
  local log="${LOG_DIR}/${name}.log"
  echo "[RUN] ${name}"
  "$@" >"${log}" 2>&1
  local code=$?
  if [[ ${code} -eq 0 ]]; then
    echo "[PASS] ${name}"
    record "${name}" "PASS" "${log}"
  else
    echo "[FAIL] ${name} (see ${log})"
    record "${name}" "FAIL" "${log}"
  fi
}

for required in "${CLI}" "${COUNTS}" "${METADATA}" "${CONTRASTS}"; do
  if [[ ! -s "${required}" ]]; then
    echo "ERROR: Missing required file: ${required}" >&2
    exit 1
  fi
done

GENE_MAPPING="${INPUTS_DIR}/airway_gene_mapping.csv"
{
  echo "gene_id,gene_symbol,ensembl_id,entrez_id"
  tail -n +2 "${COUNTS}" | cut -d, -f1 | awk 'BEGIN{n=1} {print $1 ",gene_" n "," $1 "," n; n++}'
} > "${GENE_MAPPING}"

GMT="${INPUTS_DIR}/airway_mini_sets.gmt"
mapfile -t GENES < <(awk -F, 'NR > 1 && NR <= 31 { print $1 }' "${COUNTS}")
{
  printf "AIRWAY_TOP10\tfirst ten airway genes"
  for gene in "${GENES[@]:0:10}"; do printf "\t%s" "${gene}"; done
  printf "\nAIRWAY_MID10\tmiddle ten airway genes"
  for gene in "${GENES[@]:10:10}"; do printf "\t%s" "${gene}"; done
  printf "\nAIRWAY_NEXT10\tnext ten airway genes"
  for gene in "${GENES[@]:20:10}"; do printf "\t%s" "${gene}"; done
  printf "\n"
} > "${GMT}"

MERGE="${OUT_ROOT}/results/01_merge_featurecounts"
VALIDATE="${OUT_ROOT}/results/02_validate_inputs"
FILTER="${OUT_ROOT}/results/03_filter_low_expression"
NORMALIZE="${OUT_ROOT}/results/04_normalize_transform"
PCA="${OUT_ROOT}/results/05_pca"
UMAP="${OUT_ROOT}/results/06_umap"
DE="${OUT_ROOT}/results/07_differential_expression"
PLSDA="${OUT_ROOT}/results/08_plsda"
GO="${OUT_ROOT}/results/10_go_enrichment"
PATHWAY="${OUT_ROOT}/results/11_pathway_enrichment"
GSEA="${OUT_ROOT}/results/12_gsea"
STATE_RUNS="${OUT_ROOT}/runs"
CONTRAST_ID="$(tail -n +2 "${CONTRASTS}" | head -n 1 | cut -d, -f1)"

run_step merge_featurecounts "${PYTHON_BIN}" "${CLI}" merge-featurecounts \
  --counts "${COUNTS}" \
  --outdir "${MERGE}"

run_step validate_downstream_inputs "${PYTHON_BIN}" "${CLI}" validate \
  --counts "${MERGE}/merged_raw_counts.csv" \
  --metadata "${METADATA}" \
  --contrasts "${CONTRASTS}" \
  --gene-mapping "${GENE_MAPPING}" \
  --outdir "${VALIDATE}"

run_step filter_low_expression "${PYTHON_BIN}" "${CLI}" filter \
  --counts "${VALIDATE}/validated_counts.csv" \
  --metadata "${VALIDATE}/validated_metadata.csv" \
  --outdir "${FILTER}"

run_step normalize_transform "${PYTHON_BIN}" "${CLI}" normalize \
  --counts "${FILTER}/filtered_counts.csv" \
  --metadata "${VALIDATE}/validated_metadata.csv" \
  --outdir "${NORMALIZE}"

run_step pca "${PYTHON_BIN}" "${CLI}" pca \
  --expression "${NORMALIZE}/vst_expression.csv" \
  --metadata "${VALIDATE}/validated_metadata.csv" \
  --outdir "${PCA}"

run_step umap "${PYTHON_BIN}" "${CLI}" umap \
  --expression "${NORMALIZE}/vst_expression.csv" \
  --metadata "${VALIDATE}/validated_metadata.csv" \
  --outdir "${UMAP}"

run_step differential_expression "${PYTHON_BIN}" "${CLI}" differential-expression \
  --counts "${FILTER}/filtered_counts.csv" \
  --metadata "${VALIDATE}/validated_metadata.csv" \
  --contrasts "${VALIDATE}/validated_contrasts.csv" \
  --gene-mapping "${GENE_MAPPING}" \
  --outdir "${DE}"

run_step plsda "${PYTHON_BIN}" "${CLI}" plsda \
  --expression "${NORMALIZE}/vst_expression.csv" \
  --metadata "${VALIDATE}/validated_metadata.csv" \
  --outdir "${PLSDA}"

# Build a tiny local human/Ensembl GMT from this run's significant genes so
# enrichment tests produce visible overlaps without downloading databases.
{
  printf "AIRWAY_DE_UP\tAirway significant up-regulated demo set"
  awk -F, 'NR > 1 && $4 >= 0 { printf "\t%s", $1 }' "${DE}/results/${CONTRAST_ID}/significant_genes.csv"
  printf "\nAIRWAY_DE_DOWN\tAirway significant down-regulated demo set"
  awk -F, 'NR > 1 && $4 < 0 { printf "\t%s", $1 }' "${DE}/results/${CONTRAST_ID}/significant_genes.csv"
  printf "\nAIRWAY_TOP_RANKED\tAirway top ranked demo set"
  awk -F, 'NR > 1 && NR <= 21 { printf "\t%s", $1 }' "${DE}/results/${CONTRAST_ID}/ranked_genes.csv"
  printf "\n"
} > "${GMT}"

run_step go_enrichment "${PYTHON_BIN}" "${CLI}" go-enrichment \
  --significant-genes "${DE}/results/${CONTRAST_ID}/significant_genes.csv" \
  --ranked-genes "${DE}/results/${CONTRAST_ID}/ranked_genes.csv" \
  --universe "${FILTER}/filtered_counts.csv" \
  --gmt "${GMT}" \
  --gene-mapping "${GENE_MAPPING}" \
  --organism human \
  --input-id-type ensembl_id \
  --gene-set-id-type ensembl_id \
  --outdir "${GO}"

run_step pathway_enrichment "${PYTHON_BIN}" "${CLI}" pathway-enrichment \
  --significant-genes "${DE}/results/${CONTRAST_ID}/significant_genes.csv" \
  --ranked-genes "${DE}/results/${CONTRAST_ID}/ranked_genes.csv" \
  --universe "${FILTER}/filtered_counts.csv" \
  --gmt "${GMT}" \
  --gene-mapping "${GENE_MAPPING}" \
  --organism human \
  --input-id-type ensembl_id \
  --gene-set-id-type ensembl_id \
  --outdir "${PATHWAY}"

run_step ranked_gsea "${PYTHON_BIN}" "${CLI}" gsea \
  --ranked-genes "${DE}/results/${CONTRAST_ID}/ranked_genes.csv" \
  --gmt "${GMT}" \
  --gene-mapping "${GENE_MAPPING}" \
  --organism human \
  --input-id-type ensembl_id \
  --gene-set-id-type ensembl_id \
  --outdir "${GSEA}"

run_step write_state "${PYTHON_BIN}" "${CLI}" write-state \
  --run-id "${RUN_ID}" \
  --runs-dir "${STATE_RUNS}" \
  --status completed \
  --current-step complete

run_step status "${PYTHON_BIN}" "${CLI}" status \
  --run-id "${RUN_ID}" \
  --runs-dir "${STATE_RUNS}"

echo
echo "Summary"
printf '%-32s %-6s %s\n' "Step" "Status" "Log"
printf '%-32s %-6s %s\n' "----" "------" "---"
for i in "${!NAMES[@]}"; do
  printf '%-32s %-6s %s\n' "${NAMES[$i]}" "${STATUSES[$i]}" "${LOGS[$i]}"
done
echo
echo "Output folder: ${OUT_ROOT}"

if [[ ${FAILURES} -eq 0 ]]; then
  echo "Overall: PASS"
  exit 0
fi

echo "Overall: FAIL (${FAILURES} failed)"
exit 1
