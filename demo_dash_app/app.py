from __future__ import annotations

import base64
import json
import os
import shlex
import subprocess
import shutil
import uuid
from datetime import datetime, timezone
from pathlib import Path

import dash
from dash import Input, Output, State, dcc, html
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from flask import jsonify, request, send_file

from services.job_store import JobStore
from services.upload_service import UploadService, safe_slug
from services.workflow_service import WorkflowService
from services.downstream_service import DownstreamService
from job_view import (
    log_block,
    process_status_badge,
    run_parameter_summary,
    selected_file_summary,
    strandedness_result_summary,
    workflow_plan_component,
)


APP_ROOT = Path(__file__).resolve().parent
PIPELINE_ROOT = APP_ROOT.parent / "survom-pipelines"
DEMO_REFERENCE_ROOT = PIPELINE_ROOT / "assets" / "demo_reference"
MINI_GALLUS_REFERENCE_ROOT = APP_ROOT.parent / "refs" / "gallus_gallus_ensembl116" / "mini_ref"
MAX_UPLOAD_BYTES = int(50 * 1024 * 1024 * 1024)
DEFAULT_CUTADAPT_ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
DEFAULT_CUTADAPT_ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
DEFAULT_CUTADAPT_ERROR_RATE = 0.1
DEFAULT_CUTADAPT_MINIMUM_OVERLAP = 3
ILLUMINA_TRUSEQ_HELP = (
    "Using Illumina TruSeq demo adapter defaults. For real datasets, "
    "confirm adapters from your library prep kit."
)
CUTADAPT_NOT_AVAILABLE_MESSAGE = (
    "Cutadapt is not installed in the active local environment. "
    "Choose fastp recommended, install cutadapt, or run with a container profile."
)
CUTADAPT_DOCKER_PROFILE_MESSAGE = (
    "Cutadapt is not installed in the active local environment, so this run will use "
    "the local Docker profile with the Cutadapt container."
)
SALMON_DOCKER_PROFILE_MESSAGE = (
    "Salmon is not installed in the active local environment, so this run will use "
    "the local Docker profile with the Salmon container."
)

URL_PREFIX = "/upload/"
API_PREFIX = "/upload/api"

upload_service = UploadService(APP_ROOT / "uploads", MAX_UPLOAD_BYTES)
job_store = JobStore(APP_ROOT / "runs" / "jobs.sqlite")
workflow_service = WorkflowService(APP_ROOT, PIPELINE_ROOT, upload_service, job_store)
downstream_service = DownstreamService(APP_ROOT.parent)

external_stylesheets = [dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP]

app = dash.Dash(
    __name__,
    external_stylesheets=external_stylesheets,
    suppress_callback_exceptions=True,
    requests_pathname_prefix=URL_PREFIX,
    routes_pathname_prefix=URL_PREFIX,
)

server = app.server
app.title = "SurvOm Demo Runner"


@server.before_request
def guard_dash_update_request():
    if not request.path.endswith("/_dash-update-component"):
        return

    body = request.get_json(silent=True) or {}
    output = body.get("output")
    callback_meta = app.callback_map.get(output or "")
    if output and not callback_meta:
        return jsonify({
            "error": "Dash page is using an older callback definition. Refresh the browser page and try again.",
            "output": output,
        }), 409
    if callback_meta:
        expected_count = len(callback_meta.get("inputs", [])) + len(callback_meta.get("state", []))
        actual_count = len(body.get("inputs") or []) + len(body.get("state") or [])
        if actual_count != expected_count:
            return jsonify({
                "error": "Dash page is using an older callback definition. Refresh the browser page and try again.",
                "output": output,
                "expected_inputs_and_state": expected_count,
                "received_inputs_and_state": actual_count,
            }), 409

    if os.environ.get("SURVOM_DEBUG_DASH_POST") == "1":
        print("\n--- DASH UPDATE DEBUG ---")
        print("path:", request.path)
        print("output:", body.get("output"))
        print("outputs:", body.get("outputs"))
        print("changedPropIds:", body.get("changedPropIds"))
        print("inputs count:", len(body.get("inputs") or []))
        print("state count:", len(body.get("state") or []))
        for item in body.get("inputs") or []:
            print("  Input:", item.get("id"), item.get("property"))
        for item in body.get("state") or []:
            print("  State:", item.get("id"), item.get("property"))
        print("--- END DASH UPDATE DEBUG ---\n")


@server.after_request
def prevent_dash_cache(response):
    if request.path.startswith(URL_PREFIX):
        response.headers["Cache-Control"] = "no-store, max-age=0"
    return response

steps = workflow_service.load_steps(APP_ROOT / "config" / "workflow_steps.yaml")
testers = workflow_service.load_yaml_config(APP_ROOT / "config" / "testers.yaml").get("testers", [])
omics_types = workflow_service.load_yaml_config(APP_ROOT / "config" / "omics_types.yaml").get("omics_types", {})
steps_by_id = {step["id"]: step for step in steps}


def tester_label(tester_id: str | None) -> str:
    tester = next((item for item in testers if item["id"] == tester_id), None)
    return tester["label"] if tester else (tester_id or "")


def omics_label(omics_type: str | None) -> str:
    return omics_types.get(omics_type or "", {}).get("label", omics_type or "")


def workflow_options_for_omics(omics_type: str):
    config = omics_types.get(omics_type, {})
    if not config.get("enabled"):
        return [{"label": "Coming soon", "value": "coming_soon", "disabled": True}]
    chained_ids = ["salmon_count_matrix", "star_count_matrix", "qc_trim_strandedness"]
    workflows = [workflow_id for workflow_id in config.get("workflows", []) if workflow_id in steps_by_id]
    atomic = [workflow_id for workflow_id in workflows if steps_by_id[workflow_id].get("label", "").startswith("Atomic ")]
    chained = [workflow_id for workflow_id in chained_ids if workflow_id in workflows]
    other = [workflow_id for workflow_id in workflows if workflow_id not in set(atomic) | set(chained)]

    def option(workflow_id: str) -> dict:
        return {"label": steps_by_id[workflow_id].get("label", steps_by_id[workflow_id]["name"]), "value": workflow_id}

    options = [option(workflow_id) for workflow_id in atomic + other]
    if chained:
        options.append({"label": "Chained workflows", "value": "__chained_workflows__", "disabled": True})
        options.extend(option(workflow_id) for workflow_id in chained)
    return options


def workflow_is_downstream(workflow_id: str | None) -> bool:
    return bool(workflow_id and steps_by_id.get(workflow_id, {}).get("downstream"))


def _first_path_matching(paths: list[str], keywords: tuple[str, ...]) -> str | None:
    for path in paths:
        name = Path(path).name.lower()
        if any(keyword in name for keyword in keywords):
            return path
    return None


def infer_downstream_selected_inputs(selected_files: dict[str, list[str]]) -> dict[str, str | None]:
    metadata_files = selected_files.get("metadata") or []
    other_files = selected_files.get("other") or []
    candidates = metadata_files + other_files
    return {
        "counts": _first_path_matching(candidates, ("count", "counts", "matrix", "raw_counts")),
        "metadata": _first_path_matching(candidates, ("metadata", "sample", "samples")),
        "contrasts": _first_path_matching(candidates, ("contrast", "contrasts")),
        "expression": _first_path_matching(candidates, ("vst_expression", "expression", "normalized")),
        "significant": _first_path_matching(candidates, ("significant_genes", "significant")),
        "ranked": _first_path_matching(candidates, ("ranked_genes", "ranked")),
        "universe": _first_path_matching(candidates, ("filtered_counts", "validated_counts", "universe")),
        "gmt": _first_path_matching(candidates, (".gmt", "geneset", "gene_set")),
        "gene_mapping": _first_path_matching(candidates, ("gene_mapping", "id_mapping", "annotation")),
    }


def default_airway_downstream_inputs() -> dict[str, str]:
    airway = APP_ROOT.parent / "testdatasets" / "countdata" / "test_data" / "airway"
    return {
        "counts": str(airway / "airway_raw_counts.csv"),
        "metadata": str(airway / "airway_sample_metadata.csv"),
        "contrasts": str(airway / "airway_contrasts.csv"),
    }


DOWNSTREAM_STEP_DIRS = {
    "downstream_merge_featurecounts": "01_merge_featurecounts",
    "downstream_validate_inputs": "02_validate_inputs",
    "downstream_filter_low_expression": "03_filter_low_expression",
    "downstream_normalize_transform": "04_normalize_transform",
    "downstream_pca": "05_pca",
    "downstream_umap": "06_umap",
    "downstream_differential_expression": "07_differential_expression",
    "downstream_plsda": "08_plsda",
    "downstream_go_enrichment": "10_go_enrichment",
    "downstream_pathway_enrichment": "11_pathway_enrichment",
    "downstream_ranked_gsea": "12_gsea",
}

DOWNSTREAM_INPUT_FIELDS = {
    "counts": {
        "label": "Counts matrix CSV/TSV path",
        "placeholder": "/path/to/counts.csv",
        "help": "Raw, merged, validated, or filtered count matrix.",
    },
    "metadata": {
        "label": "Sample metadata CSV/TSV path",
        "placeholder": "/path/to/metadata.csv",
        "help": "Must include sample_id and group/condition columns.",
    },
    "contrasts": {
        "label": "Contrasts CSV path",
        "placeholder": "/path/to/contrasts.csv",
        "help": "Needed only for differential-expression or full downstream runs.",
    },
    "expression": {
        "label": "Transformed expression CSV path",
        "placeholder": "/path/to/vst_expression.csv",
        "help": "Use the output from Normalize and transform for PCA, UMAP, or PLS-DA.",
    },
    "significant": {
        "label": "Significant genes CSV path",
        "placeholder": "/path/to/significant_genes.csv",
        "help": "Use differential-expression significant_genes.csv.",
    },
    "ranked": {
        "label": "Ranked genes CSV path",
        "placeholder": "/path/to/ranked_genes.csv",
        "help": "Use differential-expression ranked_genes.csv.",
    },
    "universe": {
        "label": "Tested gene universe CSV path",
        "placeholder": "/path/to/filtered_counts.csv",
        "help": "Usually filtered_counts.csv; used as the tested background.",
    },
    "gmt": {
        "label": "Local GMT gene-set file",
        "placeholder": "/path/to/gene_sets.gmt",
        "help": "Optional in Dash. If empty, a tiny demo GMT is generated for teaching tests.",
        "optional": True,
    },
    "gene_mapping": {
        "label": "Gene mapping CSV path",
        "placeholder": "/path/to/gene_mapping.csv",
        "help": "Optional. Use when count IDs and GMT IDs differ, for example Ensembl counts with symbol pathways.",
        "optional": True,
    },
}

DOWNSTREAM_REQUIRED_INPUTS = {
    "downstream_airway_all_atomics": [("counts matrix", "counts"), ("sample metadata", "metadata"), ("contrasts", "contrasts")],
    "downstream_custom_all_atomics": [("counts matrix", "counts"), ("sample metadata", "metadata"), ("contrasts", "contrasts")],
    "downstream_merge_featurecounts": [("counts matrix", "counts")],
    "downstream_validate_inputs": [("counts matrix", "counts"), ("sample metadata", "metadata"), ("contrasts", "contrasts")],
    "downstream_filter_low_expression": [("counts matrix", "counts"), ("sample metadata", "metadata")],
    "downstream_normalize_transform": [("counts matrix", "counts"), ("sample metadata", "metadata")],
    "downstream_pca": [("transformed expression", "expression"), ("sample metadata", "metadata")],
    "downstream_umap": [("transformed expression", "expression"), ("sample metadata", "metadata")],
    "downstream_differential_expression": [("counts matrix", "counts"), ("sample metadata", "metadata"), ("contrasts", "contrasts")],
    "downstream_plsda": [("transformed expression", "expression"), ("sample metadata", "metadata")],
    "downstream_go_enrichment": [("significant genes", "significant"), ("ranked genes", "ranked"), ("gene universe", "universe")],
    "downstream_pathway_enrichment": [("significant genes", "significant"), ("ranked genes", "ranked"), ("gene universe", "universe")],
    "downstream_ranked_gsea": [("ranked genes", "ranked")],
}


def downstream_visible_input_keys(workflow_id: str | None) -> set[str]:
    keys = {key for _, key in DOWNSTREAM_REQUIRED_INPUTS.get(workflow_id or "", [])}
    if workflow_id in {"downstream_go_enrichment", "downstream_pathway_enrichment", "downstream_ranked_gsea"}:
        keys.add("gmt")
        keys.add("gene_mapping")
    return keys


def latest_downstream_output(tester_id: str, omics_type: str, project_id: str, relative_names: list[str]) -> str | None:
    jobs = [
        job for job in job_store.list_jobs(limit=500)
        if job.get("tester_id") == tester_id
        and job.get("omics_type") == omics_type
        and job.get("project_id") == project_id
        and workflow_is_downstream(job.get("workflow_id"))
    ]
    for job in jobs:
        results_dir = Path(job.get("results_dir") or "")
        for relative in relative_names:
            matches = sorted(results_dir.glob(relative)) if any(char in relative for char in "*?[") else [results_dir / relative]
            for candidate in matches:
                if candidate.exists() and candidate.is_file():
                    return str(candidate)
    return None


def downstream_previous_inputs(tester_id: str, omics_type: str, project_id: str) -> dict[str, str | None]:
    return {
        "counts": latest_downstream_output(tester_id, omics_type, project_id, [
            "03_filter_low_expression/filtered_counts.csv",
            "02_validate_inputs/validated_counts.csv",
            "01_merge_featurecounts/merged_raw_counts.csv",
        ]),
        "metadata": latest_downstream_output(tester_id, omics_type, project_id, [
            "02_validate_inputs/validated_metadata.csv",
        ]),
        "contrasts": latest_downstream_output(tester_id, omics_type, project_id, [
            "02_validate_inputs/validated_contrasts.csv",
        ]),
        "expression": latest_downstream_output(tester_id, omics_type, project_id, [
            "04_normalize_transform/vst_expression.csv",
        ]),
        "significant": latest_downstream_output(tester_id, omics_type, project_id, [
            "07_differential_expression/results/*/significant_genes.csv",
        ]),
        "ranked": latest_downstream_output(tester_id, omics_type, project_id, [
            "07_differential_expression/results/*/ranked_genes.csv",
        ]),
        "universe": latest_downstream_output(tester_id, omics_type, project_id, [
            "03_filter_low_expression/filtered_counts.csv",
            "02_validate_inputs/validated_counts.csv",
        ]),
        "gmt": None,
        "gene_mapping": None,
    }


def downstream_demo_support_files(run_dir: Path, counts: str) -> dict[str, str]:
    inputs_dir = run_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    genes = []
    with Path(counts).open(encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            if len(genes) >= 30:
                break
            genes.append(line.split(",", 1)[0].strip())
    gene_mapping = inputs_dir / "dash_gene_mapping.csv"
    with gene_mapping.open("w", encoding="utf-8") as handle:
        handle.write("gene_id,gene_symbol,ensembl_id,entrez_id\n")
        for index, gene in enumerate(genes, start=1):
            handle.write(f"{gene},gene_{index},{gene},{index}\n")
    gmt = inputs_dir / "dash_mini_sets.gmt"
    with gmt.open("w", encoding="utf-8") as handle:
        handle.write("DASH_SET_1\tdash mini set 1\t" + "\t".join(genes[:10]) + "\n")
        handle.write("DASH_SET_2\tdash mini set 2\t" + "\t".join(genes[10:20]) + "\n")
        handle.write("DASH_SET_3\tdash mini set 3\t" + "\t".join(genes[20:30]) + "\n")
    return {"gene_mapping": str(gene_mapping), "gmt": str(gmt)}


def start_downstream_dash_job(
    workflow_id: str,
    tester_id: str,
    omics_type: str,
    project_id: str,
    session_id: str,
    inputs: dict[str, str],
) -> dict:
    run_id = f"{workflow_id}_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    run_dir = APP_ROOT / "runs" / safe_slug(omics_type) / safe_slug(tester_id) / safe_slug(project_id) / run_id
    results_dir = run_dir / "results"
    logs_dir = run_dir / "logs"
    run_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    request_data = {
        "run_id": run_id,
        "workflow_id": workflow_id,
        "tester_id": tester_id,
        "omics_type": omics_type,
        "project_id": project_id,
        "inputs": inputs,
        "results_dir": str(results_dir),
    }
    (run_dir / "downstream_request.json").write_text(json.dumps(request_data, indent=2), encoding="utf-8")
    support = downstream_demo_support_files(run_dir, inputs.get("counts") or default_airway_downstream_inputs()["counts"])
    cli = str(PIPELINE_ROOT / "downstream" / "bin" / "downstream_cli.py")
    step_dir = DOWNSTREAM_STEP_DIRS.get(workflow_id)
    if workflow_id in {"downstream_airway_all_atomics", "downstream_custom_all_atomics"}:
        command = [
            str(APP_ROOT.parent / "scripts" / "run_airway_downstream_atomic.sh"),
            "--counts", inputs["counts"],
            "--metadata", inputs["metadata"],
            "--contrasts", inputs["contrasts"],
            "--outdir", str(run_dir),
        ]
    elif workflow_id == "downstream_merge_featurecounts":
        command = ["python3", cli, "merge-featurecounts", "--counts", inputs["counts"], "--outdir", str(results_dir / step_dir)]
    elif workflow_id == "downstream_validate_inputs":
        command = [
            "python3", cli, "validate",
            "--counts", inputs["counts"],
            "--metadata", inputs["metadata"],
            "--contrasts", inputs["contrasts"],
            "--gene-mapping", support["gene_mapping"],
            "--outdir", str(results_dir / step_dir),
        ]
    elif workflow_id == "downstream_filter_low_expression":
        command = ["python3", cli, "filter", "--counts", inputs["counts"], "--metadata", inputs["metadata"], "--outdir", str(results_dir / step_dir)]
    elif workflow_id == "downstream_normalize_transform":
        command = ["python3", cli, "normalize", "--counts", inputs["counts"], "--metadata", inputs["metadata"], "--outdir", str(results_dir / step_dir)]
    elif workflow_id == "downstream_pca":
        command = ["python3", cli, "pca", "--expression", inputs["expression"], "--metadata", inputs["metadata"], "--outdir", str(results_dir / step_dir)]
    elif workflow_id == "downstream_umap":
        command = ["python3", cli, "umap", "--expression", inputs["expression"], "--metadata", inputs["metadata"], "--outdir", str(results_dir / step_dir)]
    elif workflow_id == "downstream_differential_expression":
        command = [
            "python3", cli, "differential-expression",
            "--counts", inputs["counts"],
            "--metadata", inputs["metadata"],
            "--contrasts", inputs["contrasts"],
            "--gene-mapping", support["gene_mapping"],
            "--outdir", str(results_dir / step_dir),
        ]
    elif workflow_id == "downstream_plsda":
        command = ["python3", cli, "plsda", "--expression", inputs["expression"], "--metadata", inputs["metadata"], "--outdir", str(results_dir / step_dir)]
    elif workflow_id == "downstream_go_enrichment":
        command = [
            "python3", cli, "go-enrichment",
            "--significant-genes", inputs["significant"],
            "--ranked-genes", inputs["ranked"],
            "--universe", inputs["universe"],
            "--gmt", inputs.get("gmt") or support["gmt"],
            *(["--gene-mapping", inputs["gene_mapping"]] if inputs.get("gene_mapping") else []),
            "--outdir", str(results_dir / step_dir),
        ]
    elif workflow_id == "downstream_pathway_enrichment":
        command = [
            "python3", cli, "pathway-enrichment",
            "--significant-genes", inputs["significant"],
            "--ranked-genes", inputs["ranked"],
            "--universe", inputs["universe"],
            "--gmt", inputs.get("gmt") or support["gmt"],
            *(["--gene-mapping", inputs["gene_mapping"]] if inputs.get("gene_mapping") else []),
            "--outdir", str(results_dir / step_dir),
        ]
    elif workflow_id == "downstream_ranked_gsea":
        command = [
            "python3", cli, "gsea",
            "--ranked-genes", inputs["ranked"],
            "--gmt", inputs.get("gmt") or support["gmt"],
            *(["--gene-mapping", inputs["gene_mapping"]] if inputs.get("gene_mapping") else []),
            "--outdir", str(results_dir / step_dir),
        ]
    else:
        raise ValueError(f"Unknown downstream workflow: {workflow_id}")
    wrapper = run_dir / "run_downstream_job.sh"
    job_record = run_dir / "job_record.json"
    quoted_command = " ".join(shlex.quote(str(part)) for part in command)
    stdout_log = shlex.quote(str(logs_dir / "stdout.log"))
    stderr_log = shlex.quote(str(logs_dir / "stderr.log"))
    wrapper.write_text(
        "#!/usr/bin/env bash\n"
        "set +e\n"
        f"{quoted_command} > {stdout_log} 2> {stderr_log}\n"
        "code=$?\n"
        f"python3 - <<'PY' \"$code\" \"{job_record}\"\n"
        "import json, sys\n"
        "from datetime import datetime, timezone\n"
        "code = int(sys.argv[1])\n"
        "path = sys.argv[2]\n"
        "json.dump({\n"
        "  'status': 'succeeded' if code == 0 else 'failed',\n"
        "  'exit_code': code,\n"
        "  'finished_at': datetime.now(timezone.utc).isoformat(),\n"
        "}, open(path, 'w'), indent=2)\n"
        "PY\n"
        "exit $code\n",
        encoding="utf-8",
    )
    wrapper.chmod(0o755)
    started_at = datetime.now(timezone.utc).isoformat()
    proc = subprocess.Popen(["bash", str(wrapper)], cwd=APP_ROOT.parent)
    job_id = f"job_{run_id}"
    metadata_json = {
        "request_path": str(run_dir / "downstream_request.json"),
        "submitted_params": {
            **inputs,
            "output_mode": "csv_json_only",
        },
        "selected_files": {
            "metadata": [value for value in inputs.values() if value],
            "fastq": [],
            "reference": [],
            "other": [],
        },
    }
    job_store.upsert_job({
        "job_id": job_id,
        "run_id": run_id,
        "session_id": session_id,
        "tester_id": tester_id,
        "tester_label": tester_label(tester_id),
        "omics_type": omics_type,
        "project_id": project_id,
        "project_label": project_label_from_id(project_id),
        "workflow_id": workflow_id,
        "status": "running",
        "pid": proc.pid,
        "run_dir": str(run_dir),
        "compiled_dir": str(run_dir),
        "results_dir": str(results_dir),
        "command": " ".join(command),
        "stdout_path": str(logs_dir / "stdout.log"),
        "stderr_path": str(logs_dir / "stderr.log"),
        "started_at": started_at,
        "metadata_json": json.dumps(metadata_json),
    })
    return {"run_id": run_id, "job_id": job_id, "run_dir": str(run_dir), "job_path": str(run_dir / "downstream_request.json"), "job_type": "downstream"}


def refresh_downstream_job(job: dict) -> dict:
    run_dir = Path(job.get("run_dir") or "")
    record_path = run_dir / "job_record.json"
    if record_path.exists():
        record = json.loads(record_path.read_text(errors="replace"))
        job_store.update_job(
            job["job_id"],
            status=record.get("status", job.get("status")),
            finished_at=record.get("finished_at"),
            exit_code=record.get("exit_code"),
        )
        return job_store.get_job(job["job_id"]) or job
    if str(job.get("status", "")).lower() == "running" and job.get("pid"):
        try:
            os.kill(int(job["pid"]), 0)
        except OSError:
            job_store.update_job(job["job_id"], status="unknown")
            return job_store.get_job(job["job_id"]) or job
    return job


def downstream_execution_record(job: dict) -> dict:
    results = Path(job.get("results_dir") or "")
    files = sorted(str(path) for path in results.rglob("*") if path.is_file()) if results.exists() else []
    steps = []
    for step_dir in sorted(path for path in results.glob("*") if path.is_dir()) if results.exists() else []:
        steps.append({
            "step_id": step_dir.name,
            "step_name": step_dir.name,
            "process_name": step_dir.name.upper(),
            "status": "completed",
            "outputs": [{"name": path.name, "files": [str(path)], "count": 1} for path in step_dir.iterdir() if path.is_file()],
            "missing_outputs": [],
            "tasks": [],
            "dependencies": [],
        })
    return {
        "run_id": job.get("run_id"),
        "status": job.get("status"),
        "step_outputs": steps,
        "all_result_files": files,
        "task_artifacts": [],
    }


def downstream_parameter_panel():
    def downstream_input_row(key: str):
        spec = DOWNSTREAM_INPUT_FIELDS[key]
        optional = " Optional." if spec.get("optional") else ""
        return html.Div(
            [
                dbc.Label(spec["label"], className="small fw-semibold"),
                dbc.Input(
                    id=f"downstream-{key}-input",
                    type="text",
                    placeholder=spec["placeholder"],
                    size="sm",
                    className="mb-1",
                ),
                html.Div(f"{spec['help']}{optional}", className="small text-muted mb-2"),
            ],
            id=f"downstream-{key}-row",
        )

    return workflow_section(
        "Downstream count matrix inputs",
        [
            html.P(
                "Runs the selected downstream atomic step and produces CSV/JSON outputs only. No plots are generated.",
                className="small text-muted",
            ),
            html.Div(id="downstream-required-help", className="small text-primary mb-2"),
            *[downstream_input_row(key) for key in DOWNSTREAM_INPUT_FIELDS],
            html.Div(
                "Tip: upload CSV/TSV files with the metadata uploader, select them in Uploaded files, and this panel will auto-fill matching paths. Previous outputs from this student/project are reused when compatible.",
                className="small text-muted",
            ),
        ],
        "downstream-params",
    )


def ensure_demo_project(tester_id: str | None, omics_type: str | None):
    if tester_id and omics_type:
        job_store.upsert_project(tester_id, omics_type, "demo_project", "demo_project")


def project_label_from_id(project_id: str | None) -> str:
    return project_id or "demo_project"


def param_input(step_id, name, spec):
    input_id = {"type": "param", "step": step_id, "name": name}
    label = spec.get("ui_label", name)
    help_text = spec.get("help")

    if spec.get("type") == "bool":
        control = dbc.Checkbox(
            id=input_id,
            value=bool(spec.get("default", False)),
            label=label,
        )
    elif spec.get("type") == "enum":
        labels = spec.get("labels", {})
        control = dbc.Select(
            id=input_id,
            options=[
                {"label": str(labels.get(choice, labels.get(str(choice), str(choice).title()))), "value": str(choice).lower()}
                for choice in spec.get("allowed", [])
            ],
            value=str(spec.get("default", "")),
            size="sm",
        )
    else:
        control = dbc.Input(
            id=input_id,
            type="number" if spec.get("type") in {"int", "number"} else "text",
            value=spec.get("default"),
            min=spec.get("min"),
            max=spec.get("max"),
            step="any" if spec.get("type") == "number" else None,
            size="sm",
        )

    if help_text:
        return html.Div([control, html.Div(help_text, className="form-text small")])
    return control


def param_group(step_id, params):
    children = []
    for name, spec in params.items():
        if spec.get("type") == "bool":
            children.append(html.Div(param_input(step_id, name, spec), className="mb-2"))
            continue
        children.append(
            html.Div(
                [
                    dbc.Label(spec.get("ui_label", name), className="small fw-semibold"),
                    param_input(step_id, name, spec),
                ],
                className="mb-2",
            )
        )
    return children


def step_card(step):
    children = [html.P(step["description"], className="text-muted small")]

    if step.get("params"):
        if "rnaseq_04_adapter_quality_trimming" in step.get("selected_steps", []):
            children.append(html.H6("Adapter and quality trimming settings", className="mt-2"))
            children.append(
                html.Div(
                    "Recommended defaults are fastp, Q20, minimum length 20, and polyG auto.",
                    className="small text-muted mb-2",
                )
            )
        children.extend(param_group(step["id"], step.get("params", {})))

    if step.get("advanced_params"):
        children.append(
            html.Details(
                [
                    html.Summary("Advanced trimming options", className="fw-semibold small mb-2"),
                    html.Div(
                        "Use Cutadapt only when you know the exact adapter sequences or need protocol-specific trimming.",
                        className="small text-muted mb-2",
                    ),
                    *param_group(step["id"], step.get("advanced_params", {})),
                ],
                className="advanced-options mt-2",
            )
        )

    children.extend(
        [
            dbc.Button(
                "Run this step",
                id={"type": "run-step", "step": step["id"]},
                color="primary",
                size="sm",
                className="mt-2",
            ),
            html.Span(
                " pending",
                id={"type": "step-status", "step": step["id"]},
                className="ms-2 badge text-bg-secondary",
            ),
        ]
    )

    return dbc.AccordionItem(children, title=step["name"], item_id=step["id"])


def upload_box(upload_id: str, title: str, detail: str, multiple: bool):
    return dcc.Upload(
        id=upload_id,
        children=html.Div(
            [
                html.Div(title, className="fw-semibold"),
                html.Div(detail, className="small text-muted"),
            ]
        ),
        className="dash-upload-zone mb-2",
        multiple=multiple,
    )


def workflow_section(title: str, children, section_id: str | None = None, wrapper_id: str | None = None):
    accordion = dbc.Accordion(
        [
            dbc.AccordionItem(
                children,
                title=title,
                item_id=section_id or title.lower().replace(" ", "-"),
            )
        ],
        start_collapsed=True,
        always_open=True,
        className="workflow-settings-accordion mb-2",
    )
    if wrapper_id:
        return html.Div(accordion, id=wrapper_id)
    return accordion


def raw_qc_parameter_panel():
    return workflow_section(
        "Raw QC settings",
        [
            html.Div(
                "FastQC runs with the default robust settings. Outputs are the FastQC HTML and ZIP reports.",
                className="small text-muted",
            ),
        ],
        "raw-qc-settings",
        "raw-qc-parameter-panel",
    )


def trimming_parameter_panel():
    return workflow_section(
        "Trimming settings",
        [
            html.Div(
                "Recommended defaults are fastp, Q20, minimum length 20, and polyG auto.",
                className="small text-muted mb-2",
            ),
            dbc.Label("Minimum base quality", className="small fw-semibold"),
            dbc.Select(
                id="quality-threshold-dropdown",
                options=[
                    {"label": "Q15", "value": "15"},
                    {"label": "Q20 recommended", "value": "20"},
                    {"label": "Q25", "value": "25"},
                    {"label": "Q30", "value": "30"},
                ],
                value="20",
                size="sm",
                className="mb-2",
            ),
            dbc.Label("Minimum read length after trimming", className="small fw-semibold"),
            dbc.Select(
                id="minimum-read-length-dropdown",
                options=[
                    {"label": "20 bp recommended", "value": "20"},
                    {"label": "25 bp", "value": "25"},
                    {"label": "30 bp", "value": "30"},
                    {"label": "50 bp", "value": "50"},
                ],
                value="20",
                size="sm",
                className="mb-2",
            ),
            dbc.Label("Trim polyG tails", className="small fw-semibold"),
            dbc.Select(
                id="trim-poly-g-dropdown",
                options=[
                    {"label": "Auto recommended", "value": "auto"},
                    {"label": "On", "value": "true"},
                    {"label": "Off", "value": "false"},
                ],
                value="auto",
                size="sm",
                className="mb-2",
            ),
            html.Details(
                [
                    html.Summary("Advanced trimming options", className="fw-semibold small mb-2"),
                    html.Div(
                        "Use Cutadapt only when you know exact adapter sequences or need protocol-specific trimming.",
                        className="small text-muted mb-2",
                    ),
                    dbc.Label("Trimming tool", className="small fw-semibold"),
                    dbc.Select(
                        id="trimming-tool-dropdown",
                        options=[
                            {"label": "fastp recommended", "value": "fastp"},
                            {"label": "Cutadapt advanced", "value": "cutadapt"},
                        ],
                        value="fastp",
                        size="sm",
                        className="mb-2",
                    ),
                    dbc.Label("Trim polyX tails", className="small fw-semibold"),
                    dbc.Select(
                        id="trim-poly-x-dropdown",
                        options=[
                            {"label": "Off recommended", "value": "false"},
                            {"label": "On", "value": "true"},
                        ],
                        value="false",
                        size="sm",
                        className="mb-2",
                    ),
                    dbc.Label("Adapter preset", className="small fw-semibold"),
                    dbc.Select(
                        id="adapter-preset-dropdown",
                        options=[
                            {"label": "Illumina TruSeq demo defaults", "value": "illumina_truseq"},
                            {"label": "Custom adapters", "value": "custom"},
                        ],
                        value="illumina_truseq",
                        size="sm",
                        className="mb-2",
                    ),
                    html.Div(id="adapter-default-message", className="small text-muted mb-2"),
                    dbc.Label("R1 adapter sequence", className="small fw-semibold"),
                    dbc.Input(id="adapter-sequence-r1-input", type="text", size="sm", className="mb-2"),
                    html.Div("Default: Illumina TruSeq R1 adapter", className="form-text small"),
                    dbc.Label("R2 adapter sequence", className="small fw-semibold"),
                    dbc.Input(id="adapter-sequence-r2-input", type="text", size="sm", className="mb-2"),
                    html.Div("Default: Illumina TruSeq R2 adapter", className="form-text small"),
                    dbc.Label("Adapter FASTA file/path", className="small fw-semibold"),
                    dbc.Input(id="adapter-fasta-input", type="text", size="sm", className="mb-2"),
                    dbc.Label("Fixed 5 prime trim R1", className="small fw-semibold"),
                    dbc.Input(id="trim-front-r1-input", type="number", min=0, max=50, value=0, size="sm", className="mb-2"),
                    dbc.Label("Fixed 5 prime trim R2", className="small fw-semibold"),
                    dbc.Input(id="trim-front-r2-input", type="number", min=0, max=50, value=0, size="sm", className="mb-2"),
                    dbc.Label("Cutadapt error rate", className="small fw-semibold"),
                    dbc.Input(id="cutadapt-error-rate-input", type="number", min=0, max=0.3, step="any", value=0.1, size="sm", className="mb-2"),
                    dbc.Label("Cutadapt minimum overlap", className="small fw-semibold"),
                    dbc.Input(id="cutadapt-minimum-overlap-input", type="number", min=1, max=20, value=3, size="sm", className="mb-2"),
                ],
                className="advanced-options mt-2",
            ),
        ],
        "trimming-settings",
        "trimming-parameter-panel",
    )


def strandedness_parameter_panel():
    return workflow_section(
        "Strandedness settings",
        [
            html.Div(
                "Runs after trimming/post-trim QC. Default is Salmon -l A on trimmed FASTQ; STAR/RSeQC is fallback validation only.",
                className="small text-muted mb-2",
            ),
            dbc.Label("Inference method", className="small fw-semibold"),
            dbc.Select(
                id="strandedness-method-dropdown",
                options=[
                    {"label": "Salmon auto-detection default", "value": "salmon_auto"},
                    {"label": "RSeQC fallback/validation", "value": "rseqc_validation"},
                ],
                value="salmon_auto",
                size="sm",
                className="mb-2",
            ),
            html.Div(
                [
                    dbc.Label("Salmon transcriptome index", className="small fw-semibold"),
                    dbc.Input(
                        id="salmon-index-input",
                        type="text",
                        placeholder="/path/to/salmon_index",
                        size="sm",
                        className="mb-2",
                    ),
                ],
                id="salmon-strandedness-fields",
            ),
            dbc.Label("Reads sampled for inference", className="small fw-semibold"),
            dbc.Input(
                id="strandedness-inference-reads-input",
                type="number",
                min=10000,
                max=5000000,
                value=1000000,
                size="sm",
                className="mb-2",
            ),
            html.Div(
                [
                    dbc.Label("Existing sorted BAM for RSeQC", className="small fw-semibold"),
                    dbc.Input(
                        id="existing-bam-input",
                        type="text",
                        placeholder="/path/to/sample.sorted.bam",
                        size="sm",
                        className="mb-2",
                    ),
                    dbc.Label("STAR genome index for RSeQC fallback", className="small fw-semibold"),
                    dbc.Input(
                        id="star-index-input",
                        type="text",
                        placeholder="/path/to/star_index",
                        size="sm",
                        className="mb-2",
                    ),
                    dbc.Label("RSeQC BED annotation", className="small fw-semibold"),
                    dbc.Input(
                        id="rseqc-ref-bed-input",
                        type="text",
                        placeholder="/path/to/genes.bed",
                        size="sm",
                        className="mb-2",
                    ),
                    dbc.Label("RSeQC strandedness threshold", className="small fw-semibold"),
                    dbc.Input(
                        id="rseqc-stranded-threshold-input",
                        type="number",
                        min=0.5,
                        max=0.95,
                        step="any",
                        value=0.6,
                        size="sm",
                        className="mb-2",
                    ),
                ],
                id="rseqc-strandedness-fields",
            ),
            dbc.Label("Manual approval / override before downstream", className="small fw-semibold"),
            dbc.Select(
                id="manual-strandedness-dropdown",
                options=[
                    {"label": "Pending auto inference", "value": "pending"},
                    {"label": "Unstranded / U", "value": "U"},
                    {"label": "Forward stranded / SF", "value": "SF"},
                    {"label": "Reverse stranded / SR", "value": "SR"},
                    {"label": "Paired inward unstranded / IU", "value": "IU"},
                    {"label": "Paired inward forward / ISF", "value": "ISF"},
                    {"label": "Paired inward reverse / ISR", "value": "ISR"},
                ],
                value="pending",
                size="sm",
                className="mb-2",
            ),
            html.Div(id="strandedness-approval-message", className="small text-muted"),
        ],
        "strandedness-settings",
        "strandedness-parameter-panel",
    )


def reference_parameter_panel():
    return dbc.Card(
        dbc.CardBody(
            [
                html.H6("Build or Validate Reference", className="mb-1"),
                html.Div(
                    "Use a lightweight demo reference for testing, select an existing bundle, or provide custom reference paths.",
                    className="small text-muted mb-2",
                ),
                html.Div(
                    "This appears here because strandedness inference needs a Salmon index by default; STAR/RSeQC validation needs a genome reference or BED annotation.",
                    className="small text-muted mb-2",
                ),
                dbc.RadioItems(
                    id="reference-mode-dropdown",
                    options=[
                        {"label": "Demo reference", "value": "demo_reference"},
                        {"label": "Existing reference bundle", "value": "prebuilt_reference"},
                        {"label": "Custom reference files", "value": "custom_reference"},
                    ],
                    value="demo_reference",
                    className="reference-choice mb-2",
                    inputClassName="me-1",
                    labelClassName="d-block mb-1",
                ),
                html.Div(
                    [
                        dbc.Label("Reference bundle path", className="small fw-semibold"),
                        dbc.Input(
                            id="reference-bundle-input",
                            type="text",
                            placeholder="/path/to/validated_reference_bundle.json or bundle directory",
                            size="sm",
                            className="mb-2",
                        ),
                    ],
                    id="reference-bundle-fields",
                ),
                html.Div(
                    [
                        dbc.Label("Route", className="small fw-semibold"),
                        dbc.Select(
                            id="selected-route-dropdown",
                            options=[
                                {"label": "Salmon route: fast transcript quantification", "value": "salmon"},
                                {"label": "STAR route: genome alignment and gene counting", "value": "star"},
                                {"label": "HISAT2 route: genome alignment", "value": "hisat2"},
                                {"label": "Custom workflow", "value": "custom"},
                            ],
                            value="salmon",
                            size="sm",
                            className="mb-2",
                        ),
                        dbc.Label("Organism", className="small fw-semibold"),
                        dbc.Input(id="organism-input", type="text", value="demo", size="sm", className="mb-2"),
                        dbc.Label("Genome build", className="small fw-semibold"),
                        dbc.Input(id="genome-build-input", type="text", value="demo_build", size="sm", className="mb-2"),
                        dbc.Label("Genome FASTA", className="small fw-semibold"),
                        dbc.Input(id="genome-fasta-input", type="text", size="sm", className="mb-2"),
                        dbc.Label("GTF annotation", className="small fw-semibold"),
                        dbc.Input(id="gtf-input", type="text", size="sm", className="mb-2"),
                        dbc.Label("Transcriptome FASTA", className="small fw-semibold"),
                        dbc.Input(id="transcriptome-fasta-input", type="text", size="sm", className="mb-2"),
                        dbc.Label("tx2gene mapping", className="small fw-semibold"),
                        dbc.Input(id="tx2gene-input", type="text", size="sm", className="mb-2"),
                        dbc.Label("HISAT2 index", className="small fw-semibold"),
                        dbc.Input(
                            id="hisat2-index-input",
                            type="text",
                            placeholder="/path/to/hisat2_index",
                            size="sm",
                            className="mb-2",
                        ),
                    ],
                    id="custom-reference-fields",
                ),
                html.Div(id="reference-status-message", className="small text-muted"),
            ],
        ),
        id="reference-parameter-panel",
        className="mb-2",
    )


def strandedness_input_panel():
    return dbc.Card(
        dbc.CardBody(
            [
                html.H6("Trimmed FASTQ input", className="mb-1"),
                html.Div(
                    "Use trimmed FASTQs from an earlier run, or manually selected trimmed FASTQ files, for atomic downstream steps.",
                    className="small text-muted mb-2",
                ),
                dbc.RadioItems(
                    id="trim-input-mode",
                    options=[
                        {"label": "Use trimmed FASTQs from previous run", "value": "previous_manifest"},
                        {"label": "Upload/select trimmed FASTQs manually", "value": "manual_trimmed_fastq"},
                    ],
                    value="previous_manifest",
                    inputClassName="me-1",
                    labelClassName="d-block mb-1",
                    className="mb-2",
                ),
                html.Div(
                    [
                        dbc.Label("Previous trim manifest", className="small fw-semibold"),
                        dbc.Select(id="trim-manifest-select", size="sm", className="mb-2"),
                        html.Div(id="trim-manifest-message", className="small text-muted"),
                    ],
                    id="previous-trim-manifest-fields",
                ),
                html.Div(
                    "For manual mode, select uploaded trimmed FASTQ files in the Uploaded files checklist above.",
                    id="manual-trimmed-fastq-message",
                    className="small text-muted",
                ),
            ]
        ),
        id="strandedness-input-panel",
        className="mb-3",
    )


def execution_parameter_panel():
    return workflow_section(
        "Execution settings",
        [
            html.Div(
                "Local is simplest. Docker is useful when Salmon, Cutadapt, STAR, or RSeQC are not installed locally.",
                className="small text-muted mb-2",
            ),
            dbc.Label("Nextflow execution profile", className="small fw-semibold"),
            dbc.Select(
                id="execution-profile-dropdown",
                options=[
                    {"label": "Local environment", "value": "local"},
                    {"label": "Local Docker containers", "value": "local_docker"},
                    {"label": "Docker profile", "value": "docker"},
                    {"label": "Singularity profile", "value": "singularity"},
                ],
                value="local",
                size="sm",
                className="mb-2",
            ),
            html.Div("Thread and memory defaults are controlled by the Nextflow process config.", className="form-text small"),
        ],
        "execution-settings",
    )


sidebar = html.Div(
    [
        html.H4("SurvOm Multi-Omics Demo", className="mb-3"),
        html.H6("Analysis setup"),
        dbc.Label("Student / tester", className="small fw-semibold"),
        dbc.Select(
            id="tester-select",
            options=[{"label": item["label"], "value": item["id"]} for item in testers],
            placeholder="Select student/tester",
            persistence=True,
            persistence_type="local",
            className="mb-2",
        ),
        dbc.Label("Omics type", className="small fw-semibold"),
        dbc.Select(
            id="omics-select",
            options=[
                {"label": config["label"], "value": omics_id}
                for omics_id, config in omics_types.items()
            ],
            value="bulk_rnaseq",
            persistence=True,
            persistence_type="local",
            className="mb-2",
        ),
        dbc.Label("Project", className="small fw-semibold"),
        dbc.Select(id="project-select", className="mb-2"),
        dbc.Input(
            id="new-project-name",
            type="text",
            placeholder="New project name",
            size="sm",
            className="mb-2",
        ),
        dbc.Button(
            "Create project",
            id="create-project",
            color="secondary",
            size="sm",
            className="mb-3 w-100",
        ),
        html.Div(id="project-message", className="small mb-2"),
        dbc.Label("Workflow", className="small fw-semibold"),
        dbc.Select(id="workflow-select", className="mb-3"),
        html.Div(id="omics-message", className="mb-3"),

        html.Div(id="upload-placeholder", className="mb-2"),
        html.Div(
            [
                html.H6("Upload FASTQ files"),
                upload_box("fastq-upload", "Select or drop FASTQ files here", ".fastq, .fastq.gz, .fq, .fq.gz", True),
                html.H6("Upload sample sheet / metadata"),
                upload_box("metadata-upload", "Select/drop sample sheet, metadata, counts, or contrasts CSV/TSV here", ".csv, .tsv, .xlsx", False),
                html.H6("Upload reference files"),
                upload_box(
                    "reference-upload",
                    "Select or drop reference FASTA/GTF/tx2gene/index archive here",
                    ".fa, .fa.gz, .gtf, .gtf.gz, .tsv, .zip, .tar.gz",
                    True,
                ),
            ],
            id="upload-panel",
        ),
        html.Div(id="upload-message", className="small mb-2"),

        dbc.Alert(
            "This visible uploader is for demo/test files. Production FASTQ uploads should use the tus plan.",
            color="info",
            className="small",
        ),

        html.H6("Uploaded files"),
        dbc.Input(
            id="file-search",
            type="text",
            placeholder="Search filename",
            size="sm",
            className="mb-2",
        ),
        dbc.Select(
            id="file-category-filter",
            options=[
                {"label": "All categories", "value": "all"},
                {"label": "FASTQ only", "value": "fastq"},
                {"label": "Metadata only", "value": "metadata"},
                {"label": "Reference only", "value": "reference"},
                {"label": "Other only", "value": "other"},
            ],
            value="all",
            size="sm",
            className="mb-2",
        ),
        html.Div(id="uploaded-files-list", className="small mb-3"),
        dcc.Checklist(
            id="selected-files-checklist",
            options=[],
            value=[],
            className="file-checklist small mb-3",
            inputClassName="me-1",
            labelClassName="d-block mb-1",
        ),

        html.H6("Workflow setup"),
        html.Div(downstream_parameter_panel(), id="downstream-panel", className="mb-3"),
        html.Div(strandedness_input_panel(), id="strandedness-input-wrapper", className="mb-3"),
        html.Div(reference_parameter_panel(), id="reference-panel", className="mb-3"),
        dbc.Switch(
            id="advanced-options-toggle",
            label="Advanced options",
            value=False,
            className="mb-2",
        ),
        html.Div(
            [
                html.Div(raw_qc_parameter_panel(), id="raw-qc-panel", className="mb-2"),
                html.Div(trimming_parameter_panel(), id="parameter-panel", className="mb-2"),
                html.Div(strandedness_parameter_panel(), id="strandedness-panel", className="mb-2"),
                html.Div(execution_parameter_panel(), id="execution-panel", className="mb-3"),
            ],
            id="advanced-workflow-options",
            className="mb-2",
        ),

        dbc.Button(
            "Run analysis",
            id="run-analysis",
            color="success",
            className="mt-2 w-100",
        ),
    ],
    id="sidebar",
    className="sidebar p-3",
)


content = html.Div(
    [
        dcc.Store(id="current-session", storage_type="local"),
        dcc.Store(id="workspace-store", storage_type="local"),
        dcc.Store(id="current-job"),

        # Faster polling gives better demo feedback while jobs are running.
        dcc.Interval(id="poller", interval=5000, n_intervals=0),

        html.H2("SurvOm Demo Execution"),
        html.P(
            "Choose a student workspace, select an omics module and workflow, upload files, then inspect logs and outputs.",
            className="text-muted",
        ),
        html.Div(id="run-message", className="mb-3"),

        dbc.Row(
            [
                dbc.Col(
                    dbc.Card(
                        [
                            dbc.CardHeader("Selected workflow details"),
                            dbc.CardBody(id="selected-step-details"),
                        ]
                    ),
                    md=6,
                ),
                dbc.Col(
                    dbc.Card(
                        [
                            dbc.CardHeader("Job status"),
                            dbc.CardBody(
                                dcc.Loading(
                                    id="job-status-loading",
                                    type="circle",
                                    children=html.Div(id="job-status-panel"),
                                )
                            ),
                        ]
                    ),
                    md=6,
                ),
            ],
            className="g-3",
        ),

        dbc.Row(
            [
                dbc.Col(
                    dbc.Card(
                        [
                            dbc.CardHeader("Logs"),
                            dbc.CardBody(
                                dcc.Loading(
                                    id="logs-loading",
                                    type="circle",
                                    children=html.Pre(
                                        id="logs-panel",
                                        className="log-panel",
                                    ),
                                )
                            ),
                        ]
                    ),
                    md=6,
                ),
                dbc.Col(
                    dbc.Card(
                        [
                            dbc.CardHeader("Output files"),
                            dbc.CardBody(
                                dcc.Loading(
                                    id="outputs-loading",
                                    type="circle",
                                    children=html.Div(id="outputs-panel"),
                                )
                            ),
                        ]
                    ),
                    md=6,
                ),
            ],
            className="g-3 mt-1",
        ),
    ],
    id="content",
    className="content p-4",
)


app.layout = html.Div([sidebar, content], className="app-shell")


@server.post(f"{API_PREFIX}/upload/chunk/<upload_type>/<session_id>")
def upload_chunk(upload_type, session_id):
    try:
        file_storage = request.files.get("file")

        if file_storage is None:
            return jsonify({"error": "Missing upload chunk"}), 400

        result = upload_service.save_chunk(
            session_id,
            upload_type,
            request.form or request.args,
            file_storage,
            request.values.get("tester_id"),
            request.values.get("omics_type"),
            request.values.get("project_id"),
        )

        return jsonify(result)

    except Exception as exc:
        return jsonify({"error": str(exc)}), 400


@server.post(f"{API_PREFIX}/upload/complete/<upload_type>/<session_id>")
def upload_complete(upload_type, session_id):
    try:
        result = upload_service.assemble_file(
            session_id,
            upload_type,
            request.form or request.args,
            request.values.get("tester_id"),
            request.values.get("omics_type"),
            request.values.get("project_id"),
        )

        return jsonify(result)

    except Exception as exc:
        return jsonify({"error": str(exc)}), 400


@server.post(f"{API_PREFIX}/upload/simple/<upload_type>/<session_id>")
def upload_simple(upload_type, session_id):
    try:
        uploaded_files = request.files.getlist("files")

        if not uploaded_files:
            return jsonify({"error": "Missing uploaded files"}), 400

        results = [
            upload_service.save_file(
                session_id,
                upload_type,
                file_storage,
                request.values.get("tester_id"),
                request.values.get("omics_type"),
                request.values.get("project_id"),
            )
            for file_storage in uploaded_files
        ]

        return jsonify({"files": results})

    except Exception as exc:
        return jsonify({"error": str(exc)}), 400


@server.get(f"{API_PREFIX}/uploads/<session_id>")
def list_uploads(session_id):
    try:
        return jsonify(
            upload_service.list_uploads(
                session_id,
                request.args.get("tester_id"),
                request.args.get("omics_type"),
                request.args.get("project_id"),
            )
        )

    except Exception as exc:
        return jsonify({"error": str(exc)}), 400


@server.post(f"{API_PREFIX}/jobs/<job_id>/start")
def start_job(job_id):
    try:
        return jsonify(workflow_service.start_job(job_id))

    except Exception as exc:
        return jsonify({"error": str(exc)}), 400


@server.get(f"{API_PREFIX}/jobs/<job_id>")
def get_job(job_id):
    job = workflow_service.refresh_job(job_id)

    if not job:
        return jsonify({"error": "Job not found"}), 404

    return jsonify(job)


@server.get(f"{API_PREFIX}/jobs")
def list_jobs():
    session_id = request.args.get("session_id")
    tester_id = request.args.get("tester_id")
    omics_type = request.args.get("omics_type")
    project_id = request.args.get("project_id")
    if session_id and tester_id and omics_type:
        return jsonify(job_store.list_workspace_jobs(session_id, tester_id, omics_type, project_id=project_id))
    return jsonify([])


@server.get(f"{API_PREFIX}/jobs/<job_id>/outputs")
def get_outputs(job_id):
    return jsonify(
        {
            "job_id": job_id,
            "files": workflow_service.output_files(job_id),
        }
    )


@server.get(f"{API_PREFIX}/file")
def get_file():
    path = Path(request.args.get("path", ""))
    session_id = request.args.get("session_id", "")
    tester_id = request.args.get("tester_id", "")
    omics_type = request.args.get("omics_type", "")
    project_id = request.args.get("project_id", "demo_project")
    resolved = path.resolve()

    allowed_roots = []
    if session_id and tester_id and omics_type:
        allowed_roots.append(upload_service.session_dir(session_id, tester_id, omics_type, project_id).resolve())
        allowed_roots.append((APP_ROOT / "uploads" / safe_slug(omics_type) / safe_slug(tester_id) / safe_slug(project_id)).resolve())
        for job in job_store.list_workspace_jobs(session_id, tester_id, omics_type, project_id=project_id, limit=500):
            if job.get("run_dir"):
                allowed_roots.append(Path(job["run_dir"]).resolve())
            if job.get("compiled_dir"):
                allowed_roots.append(Path(job["compiled_dir"]).resolve())

    if not any(is_relative_to(resolved, root) for root in allowed_roots):
        return jsonify({"error": "File path is not allowed"}), 403

    if not resolved.exists() or not resolved.is_file():
        return jsonify({"error": "File not found"}), 404

    return send_file(resolved)


def is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def selected_files_by_category(selected_paths: list[str], tester_id: str, omics_type: str, project_id: str) -> dict[str, list[str]]:
    root = (APP_ROOT / "uploads" / safe_slug(omics_type) / safe_slug(tester_id) / safe_slug(project_id)).resolve()
    selected = {"fastq": [], "metadata": [], "other": []}
    for raw_path in selected_paths or []:
        path = Path(raw_path).resolve()
        if not is_relative_to(path, root) or not path.is_file():
            raise ValueError(f"Selected file is outside this project workspace: {raw_path}")
        parts = path.relative_to(root).parts
        if len(parts) >= 2 and parts[1] in selected:
            selected[parts[1]].append(str(path))
    return selected


def collect_params(step_id, values):
    params = {}

    for item in values:
        field_id = item["id"]

        if isinstance(field_id, dict) and field_id.get("step") == step_id:
            params[field_id["name"]] = item["value"]

    return params


def workflow_has_trimming(workflow_id: str) -> bool:
    step = next((s for s in steps if s["id"] == workflow_id), None)
    return bool(step and "rnaseq_04_adapter_quality_trimming" in step.get("selected_steps", []))


def workflow_has_strandedness(workflow_id: str) -> bool:
    step = next((s for s in steps if s["id"] == workflow_id), None)
    return bool(step and "rnaseq_06_strandedness_inference" in step.get("selected_steps", []))


def workflow_has_reference(workflow_id: str) -> bool:
    step = next((s for s in steps if s["id"] == workflow_id), None)
    return bool(step and "rnaseq_06a_reference_build_validation" in step.get("selected_steps", []))


def workflow_has_raw_qc(workflow_id: str) -> bool:
    step = next((s for s in steps if s["id"] == workflow_id), None)
    return bool(step and "rnaseq_03_raw_read_qc" in step.get("selected_steps", []))


def workflow_is_strandedness_only(workflow_id: str) -> bool:
    step = next((s for s in steps if s["id"] == workflow_id), None)
    selected = set(step.get("selected_steps", [])) if step else set()
    return bool(
        "rnaseq_06_strandedness_inference" in selected
        and "rnaseq_04_adapter_quality_trimming" not in selected
        and "rnaseq_03_raw_read_qc" not in selected
    )


def workflow_uses_existing_trimmed_reads(workflow_id: str) -> bool:
    step = next((s for s in steps if s["id"] == workflow_id), None)
    selected = set(step.get("selected_steps", [])) if step else set()
    trimmed_consumers = {
        "rnaseq_05_post_trim_quality_control",
        "rnaseq_06_strandedness_inference",
        "rnaseq_07_salmon_quantification",
        "rnaseq_09_star_alignment",
        "rnaseq_09d_hisat2_alignment",
    }
    return bool(trimmed_consumers & selected) and "rnaseq_04_adapter_quality_trimming" not in selected


def workflow_needs_fastq_input(workflow_id: str) -> bool:
    step = next((s for s in steps if s["id"] == workflow_id), None)
    return workflow_service.needs_fastq_input(list(step.get("selected_steps", []))) if step else True


def normalize_empty(value):
    return None if value == "" else value


def uses_container_profile(profile: str | None) -> bool:
    return bool({token.strip() for token in (profile or "").split(",")} & {"docker", "singularity", "aws"})


def docker_available() -> bool:
    return shutil.which("docker") is not None


def resolve_execution_profile(params: dict) -> tuple[str, str | None]:
    requested = params.get("execution_profile") or "local"
    needs_cutadapt = params.get("trimming_tool") == "cutadapt"
    needs_salmon = params.get("strandedness_method") == "salmon_auto" or params.get("selected_route") == "salmon"
    needs_rseqc = params.get("strandedness_method") == "rseqc_validation"
    if not needs_cutadapt and not needs_salmon:
        if not needs_rseqc:
            return requested, None
    if uses_container_profile(requested):
        return requested, None
    missing = []
    if needs_cutadapt and not shutil.which("cutadapt"):
        missing.append(CUTADAPT_DOCKER_PROFILE_MESSAGE)
    if needs_salmon and not shutil.which("salmon"):
        missing.append(SALMON_DOCKER_PROFILE_MESSAGE)
    if needs_rseqc and (not shutil.which("infer_experiment.py") or (not params.get("existing_bam") and not shutil.which("STAR"))):
        missing.append("STAR/RSeQC are not installed in the active local environment, so this run will use the local Docker profile with STAR/RSeQC containers.")
    if not missing:
        return requested, None
    if docker_available():
        return "local_docker", " ".join(missing)
    return requested, None


def resolve_cutadapt_defaults(params: dict) -> dict:
    resolved = dict(params)
    if resolved.get("trimming_tool") != "cutadapt":
        if resolved.get("adapter_preset") == "illumina_truseq":
            resolved["adapter_sequence_r1"] = None
            resolved["adapter_sequence_r2"] = None
            resolved["adapter_source"] = None
        return resolved
    preset = resolved.get("adapter_preset") or "illumina_truseq"
    resolved["adapter_preset"] = preset
    if preset == "illumina_truseq":
        adapter_r1 = resolved.get("adapter_sequence_r1")
        adapter_r2 = resolved.get("adapter_sequence_r2")
        resolved["adapter_sequence_r1"] = adapter_r1 or DEFAULT_CUTADAPT_ADAPTER_R1
        resolved["adapter_sequence_r2"] = adapter_r2 or DEFAULT_CUTADAPT_ADAPTER_R2
        resolved["adapter_source"] = (
            "user"
            if adapter_r1 and adapter_r1 != DEFAULT_CUTADAPT_ADAPTER_R1
            or adapter_r2 and adapter_r2 != DEFAULT_CUTADAPT_ADAPTER_R2
            else "default"
        )
    else:
        resolved["adapter_source"] = "user" if resolved.get("adapter_sequence_r1") or resolved.get("adapter_sequence_r2") else "empty"
    resolved["cutadapt_error_rate"] = resolved.get("cutadapt_error_rate") or DEFAULT_CUTADAPT_ERROR_RATE
    resolved["cutadapt_minimum_overlap"] = resolved.get("cutadapt_minimum_overlap") or DEFAULT_CUTADAPT_MINIMUM_OVERLAP
    return resolved


def resolve_reference_defaults(params: dict) -> dict:
    resolved = dict(params)
    if resolved.get("reference_mode") != "demo_reference":
        return resolved

    mini_fasta = MINI_GALLUS_REFERENCE_ROOT / "Gallus_gallus.mini.fa.gz"
    mini_gtf = MINI_GALLUS_REFERENCE_ROOT / "Gallus_gallus.mini.gtf.gz"
    mini_tx2gene = MINI_GALLUS_REFERENCE_ROOT / "tx2gene.tsv"
    resolved["organism"] = resolved.get("organism") or "Gallus_gallus"
    resolved["genome_build"] = resolved.get("genome_build") or "Ensembl_116_mini_GRCg7b"
    resolved["genome_fasta"] = resolved.get("genome_fasta") or str(mini_fasta)
    resolved["gtf"] = resolved.get("gtf") or str(mini_gtf)
    resolved["transcriptome_fasta"] = resolved.get("transcriptome_fasta") or str(mini_fasta)
    resolved["tx2gene"] = resolved.get("tx2gene") or str(mini_tx2gene)
    resolved["salmon_index"] = resolved.get("salmon_index") or None
    resolved["star_index"] = resolved.get("star_index") or None
    resolved["hisat2_index"] = resolved.get("hisat2_index") or None
    return resolved


def validate_run_params(
    workflow_id: str,
    params: dict,
    session_id: str,
    tester_id: str,
    omics_type: str,
    selected_files: dict[str, list[str]] | None = None,
) -> list[str]:
    params = resolve_cutadapt_defaults(params)
    params = resolve_reference_defaults(params)

    errors = []
    has_trimming = workflow_has_trimming(workflow_id)
    if has_trimming:
        try:
            quality = int(params.get("quality_threshold", 20))
            if quality not in {15, 20, 25, 30}:
                errors.append("Minimum base quality must be Q15, Q20, Q25, or Q30.")
        except (TypeError, ValueError):
            errors.append("Minimum base quality must be an integer.")

        try:
            min_len = int(params.get("minimum_read_length", 20))
            if min_len not in {20, 25, 30, 50}:
                errors.append("Minimum read length must be 20, 25, 30, or 50 bp.")
        except (TypeError, ValueError):
            errors.append("Minimum read length must be an integer.")

        if str(params.get("trim_poly_g", "auto")).lower() not in {"auto", "true", "false"}:
            errors.append("Trim polyG tails must be Auto, True, or False.")
        if str(params.get("trimming_tool", "fastp")).lower() not in {"fastp", "cutadapt"}:
            errors.append("Trimming tool must be fastp or Cutadapt.")

        for field in ("trim_front_r1", "trim_front_r2"):
            try:
                if int(params.get(field, 0) or 0) < 0:
                    errors.append(f"{field} must be zero or greater.")
            except (TypeError, ValueError):
                errors.append(f"{field} must be an integer.")

    if has_trimming and params.get("trimming_tool") == "cutadapt":
        if not shutil.which("cutadapt") and not uses_container_profile(params.get("execution_profile")) and not docker_available():
            errors.append(CUTADAPT_NOT_AVAILABLE_MESSAGE)
        if not params.get("adapter_sequence_r1"):
            errors.append("Cutadapt requires an R1 adapter sequence.")
        selected_fastqs = [Path(path).name for path in (selected_files or {}).get("fastq", [])]
        if not selected_fastqs:
            uploads = upload_service.list_uploads(session_id, tester_id, omics_type, params.get("project_id"))
            selected_fastqs = [item["name"] for item in uploads["fastq"]]
        has_r2 = any("_R2" in name or "_2" in name for name in selected_fastqs)
        if has_r2 and not params.get("adapter_sequence_r2"):
            errors.append("Cutadapt requires an R2 adapter sequence for paired-end uploads.")
        try:
            error_rate = float(params.get("cutadapt_error_rate", 0.1) or 0.1)
            if not 0 <= error_rate <= 0.3:
                errors.append("Cutadapt error rate must be between 0 and 0.3.")
        except (TypeError, ValueError):
            errors.append("Cutadapt error rate must be numeric.")
        try:
            overlap = int(params.get("cutadapt_minimum_overlap", 3) or 3)
            if not 1 <= overlap <= 20:
                errors.append("Cutadapt minimum overlap must be between 1 and 20.")
        except (TypeError, ValueError):
            errors.append("Cutadapt minimum overlap must be an integer.")

    if workflow_has_strandedness(workflow_id):
        method = params.get("strandedness_method") or "salmon_auto"
        if method not in {"salmon_auto", "rseqc_validation"}:
            errors.append("Strandedness method must be Salmon auto-detection or RSeQC fallback/validation.")
        if method == "salmon_auto" and not params.get("salmon_index") and not workflow_has_reference(workflow_id):
            errors.append("Provide a Salmon transcriptome index path for Salmon -l A strandedness inference.")
        if method == "salmon_auto" and params.get("salmon_index"):
            salmon_index = Path(str(params["salmon_index"])).expanduser()
            if not (salmon_index / "versionInfo.json").exists():
                errors.append(
                    "Selected Salmon index is not valid: missing versionInfo.json. "
                    "Build a real Salmon index or choose a valid existing reference bundle."
                )
        if method == "rseqc_validation":
            if not params.get("rseqc_ref_bed"):
                errors.append("Provide an RSeQC BED annotation file for infer_experiment.py.")
            if not params.get("existing_bam") and not params.get("star_index"):
                errors.append("Provide an existing sorted BAM, or provide a STAR genome index so the fallback can align trimmed FASTQ first.")
            try:
                threshold = float(params.get("rseqc_stranded_threshold", 0.6) or 0.6)
                if not 0.5 <= threshold <= 0.95:
                    errors.append("RSeQC strandedness threshold must be between 0.5 and 0.95.")
            except (TypeError, ValueError):
                errors.append("RSeQC strandedness threshold must be numeric.")
        try:
            read_count = int(params.get("strandedness_inference_reads", 1000000) or 1000000)
            if not 10000 <= read_count <= 5000000:
                errors.append("Reads sampled for strandedness inference must be between 10,000 and 5,000,000.")
        except (TypeError, ValueError):
            errors.append("Reads sampled for strandedness inference must be an integer.")

    if workflow_has_reference(workflow_id):
        mode = params.get("reference_mode") or "demo_reference"
        route = params.get("selected_route") or "salmon"
        if workflow_id in {"salmon_quant_only", "salmon_count_matrix", "qc_trim_strandedness"}:
            route = "salmon"
        elif workflow_id in {"star_align_only", "star_count_matrix", "star_htseq_route"}:
            route = "star"
        elif workflow_id in {"hisat2_align_only", "hisat2_alignment_route", "hisat2_featurecounts_route", "hisat2_htseq_route"}:
            route = "hisat2"
        if route not in {"salmon", "star", "hisat2", "custom"}:
            errors.append("Reference route must be Salmon, STAR, HISAT2, or Custom workflow.")
        salmon_route = route == "salmon" or workflow_id in {"salmon_quant_only", "salmon_count_matrix", "qc_trim_strandedness"}
        star_route = route == "star" or workflow_id in {"star_align_only", "star_count_matrix", "star_htseq_route"}
        hisat2_route = route == "hisat2" or workflow_id in {"hisat2_align_only", "hisat2_alignment_route", "hisat2_featurecounts_route", "hisat2_htseq_route"}
        if salmon_route and params.get("salmon_index"):
            salmon_index = Path(str(params["salmon_index"])).expanduser()
            if not (salmon_index / "versionInfo.json").exists():
                errors.append("Salmon reference index is missing versionInfo.json; build/provide a real Salmon index.")
        if star_route and params.get("star_index"):
            star_index = Path(str(params["star_index"])).expanduser()
            if not (star_index / "Genome").exists():
                errors.append("STAR reference index is missing Genome; build/provide a real STAR index.")
        if hisat2_route and params.get("hisat2_index"):
            hisat2_index = Path(str(params["hisat2_index"])).expanduser()
            if not ((hisat2_index / "genome.1.ht2").exists() or (hisat2_index / "genome.1.ht2l").exists()):
                errors.append("HISAT2 reference index is missing genome.1.ht2/genome.1.ht2l; build/provide a real HISAT2 index.")
        if mode == "demo_reference":
            return errors
        if not params.get("organism"):
            errors.append("Reference organism metadata is required.")
        if not params.get("genome_build"):
            errors.append("Reference genome build metadata is required.")
        if salmon_route and not (params.get("transcriptome_fasta") and params.get("tx2gene")):
            errors.append("Salmon route requires transcriptome FASTA and tx2gene. Salmon index is optional if it can be built.")
        if star_route and not (params.get("genome_fasta") and params.get("gtf")):
            errors.append("STAR route requires genome FASTA and GTF. STAR index is optional if it can be built.")
        if hisat2_route and not params.get("genome_fasta"):
            errors.append("HISAT2 route requires genome FASTA. HISAT2 index is optional if it can be built.")

    return list(dict.fromkeys(errors))


def decode_dash_upload(contents: str) -> bytes:
    if not contents or "," not in contents:
        raise ValueError("Invalid upload content")

    _, encoded = contents.split(",", 1)
    return base64.b64decode(encoded)


def file_list_component(
    tester_id: str,
    omics_type: str,
    project_id: str,
    search: str | None,
    category_filter: str,
    current_selection: list[str] | None,
):
    records = job_store.list_upload_records(tester_id, omics_type, project_id)
    search = (search or "").lower()
    if category_filter and category_filter != "all":
        records = [record for record in records if record.get("category") == category_filter]
    if search:
        records = [record for record in records if search in (record.get("filename") or "").lower()]

    if not records:
        return html.Div("No uploaded files yet.", className="text-muted"), [], []

    options = []
    default_values = []
    available_values = set()
    for record in records:
        path = record["stored_path"]
        available_values.add(path)
        label = (
            f"{record['filename']} | {record['category']} | {record['size_bytes']} bytes | "
            f"{record['created_at']} | session {record['session_id']}"
        )
        options.append({"label": label, "value": path})
        if record["category"] in {"fastq", "metadata"}:
            default_values.append(path)

    info = html.Div(f"{len(records)} file(s) in this project. Select files to use for the next run.", className="text-muted")
    preserved = [path for path in (current_selection or []) if path in available_values]
    return info, options, preserved or default_values


def infer_reference_inputs_from_uploads(tester_id: str, omics_type: str, project_id: str) -> dict[str, str | None]:
    records = [
        record
        for record in job_store.list_upload_records(tester_id, omics_type, project_id)
        if record.get("category") == "reference"
    ]
    inferred = {
        "reference_bundle": None,
        "genome_fasta": None,
        "gtf": None,
        "transcriptome_fasta": None,
        "tx2gene": None,
    }
    fasta_candidates = []
    for record in records:
        name = (record.get("filename") or "").lower()
        path = record.get("stored_path")
        if not path:
            continue
        if name.endswith((".json", ".zip", ".tar", ".tar.gz", ".tgz")) and not inferred["reference_bundle"]:
            inferred["reference_bundle"] = path
        if name.endswith((".gtf", ".gtf.gz", ".gff", ".gff.gz", ".gff3", ".gff3.gz")) and not inferred["gtf"]:
            inferred["gtf"] = path
        if name.endswith((".tsv", ".csv")) and any(token in name for token in ("tx2gene", "transcript", "gene")) and not inferred["tx2gene"]:
            inferred["tx2gene"] = path
        if name.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
            fasta_candidates.append((name, path))
    for name, path in fasta_candidates:
        if any(token in name for token in ("transcript", "cdna", "rna")) and not inferred["transcriptome_fasta"]:
            inferred["transcriptome_fasta"] = path
        elif not inferred["genome_fasta"]:
            inferred["genome_fasta"] = path
    if not inferred["transcriptome_fasta"] and inferred["genome_fasta"]:
        inferred["transcriptome_fasta"] = inferred["genome_fasta"]
    if not inferred["genome_fasta"] and inferred["transcriptome_fasta"]:
        inferred["genome_fasta"] = inferred["transcriptome_fasta"]
    return inferred


app.clientside_callback(
    """
    function(_, current) {
      const key = "survom_demo_session_id";
      let sid = window.localStorage.getItem(key);

      if (!sid) {
        sid = (window.crypto && window.crypto.randomUUID)
          ? window.crypto.randomUUID()
          : "session_" + Date.now() + "_" + Math.random().toString(16).slice(2);

        window.localStorage.setItem(key, sid);
      }

      return {"session_id": sid};
    }
    """,
    Output("current-session", "data"),
    Input("poller", "n_intervals"),
    State("current-session", "data"),
)


app.clientside_callback(
    """
    function(testerId, omicsType) {
      if (testerId) {
        window.localStorage.setItem("survom_demo_tester_id", testerId);
      }
      if (omicsType) {
        window.localStorage.setItem("survom_demo_omics_type", omicsType);
      }
      return {"tester_id": testerId || null, "omics_type": omicsType || "bulk_rnaseq"};
    }
    """,
    Output("workspace-store", "data"),
    Input("tester-select", "value"),
    Input("omics-select", "value"),
)


@app.callback(
    Output("project-select", "options"),
    Output("project-select", "value"),
    Output("project-message", "children"),
    Input("tester-select", "value"),
    Input("omics-select", "value"),
    Input("create-project", "n_clicks"),
    State("new-project-name", "value"),
    State("project-select", "value"),
)
def configure_projects(tester_id, omics_type, create_clicks, new_project_name, current_project):
    if not tester_id or not omics_type:
        return [], None, "Select a student/tester and omics type first."
    ensure_demo_project(tester_id, omics_type)
    message = ""
    ctx = dash.callback_context
    if ctx.triggered and ctx.triggered[0]["prop_id"].startswith("create-project") and new_project_name:
        project_id = safe_slug(new_project_name, default="demo_project")
        job_store.upsert_project(tester_id, omics_type, project_id, new_project_name.strip())
        current_project = project_id
        message = f"Created project: {new_project_name.strip()}"
    projects = job_store.list_projects(tester_id, omics_type)
    options = [
        {"label": item["project_label"], "value": item["project_id"]}
        for item in projects
    ]
    value = current_project if current_project in {item["project_id"] for item in projects} else "demo_project"
    return options, value, message


@app.callback(
    Output("workflow-select", "options"),
    Output("workflow-select", "value"),
    Output("omics-message", "children"),
    Output("upload-placeholder", "children"),
    Output("upload-panel", "style"),
    Input("omics-select", "value"),
)
def configure_omics(omics_type):
    omics_type = omics_type or "bulk_rnaseq"
    config = omics_types.get(omics_type, {})
    enabled = bool(config.get("enabled"))
    options = workflow_options_for_omics(omics_type)
    workflow_value = config.get("default_workflow") if enabled else "coming_soon"
    message = dbc.Alert(
        config.get("message") or "This omics module is prepared in the UI but the workflow is not active yet.",
        color="success" if enabled else "warning",
        className="small",
    )
    if not enabled:
        upload_placeholder = dbc.Alert(
            f"{config.get('label', omics_type)} upload support is under development.",
            color="secondary",
            className="small",
        )
        upload_style = {"display": "none"}
    else:
        upload_placeholder = None
        upload_style = {}
    return options, workflow_value, message, upload_placeholder, upload_style


@app.callback(
    Output("downstream-panel", "style"),
    Output("downstream-required-help", "children"),
    Output("downstream-counts-row", "style"),
    Output("downstream-metadata-row", "style"),
    Output("downstream-contrasts-row", "style"),
    Output("downstream-expression-row", "style"),
    Output("downstream-significant-row", "style"),
    Output("downstream-ranked-row", "style"),
    Output("downstream-universe-row", "style"),
    Output("downstream-gmt-row", "style"),
    Output("downstream-gene_mapping-row", "style"),
    Input("workflow-select", "value"),
    Input("omics-select", "value"),
)
def show_downstream_parameters(workflow_id, omics_type):
    if omics_type == "bulk_rnaseq" and workflow_is_downstream(workflow_id):
        visible_keys = downstream_visible_input_keys(workflow_id)
        help_text = "Required for this step: " + ", ".join(
            DOWNSTREAM_INPUT_FIELDS[key]["label"].replace(" path", "")
            for key in DOWNSTREAM_INPUT_FIELDS
            if key in visible_keys and not DOWNSTREAM_INPUT_FIELDS[key].get("optional")
        )
        styles = [({} if key in visible_keys else {"display": "none"}) for key in DOWNSTREAM_INPUT_FIELDS]
        return {"display": "block"}, help_text, *styles
    hidden = {"display": "none"}
    return hidden, "", hidden, hidden, hidden, hidden, hidden, hidden, hidden, hidden, hidden


@app.callback(
    Output("downstream-counts-input", "value"),
    Output("downstream-metadata-input", "value"),
    Output("downstream-contrasts-input", "value"),
    Output("downstream-expression-input", "value"),
    Output("downstream-significant-input", "value"),
    Output("downstream-ranked-input", "value"),
    Output("downstream-universe-input", "value"),
    Output("downstream-gmt-input", "value"),
    Output("downstream-gene_mapping-input", "value"),
    Input("workflow-select", "value"),
    Input("selected-files-checklist", "value"),
    State("tester-select", "value"),
    State("omics-select", "value"),
    State("project-select", "value"),
    State("downstream-counts-input", "value"),
    State("downstream-metadata-input", "value"),
    State("downstream-contrasts-input", "value"),
    State("downstream-expression-input", "value"),
    State("downstream-significant-input", "value"),
    State("downstream-ranked-input", "value"),
    State("downstream-universe-input", "value"),
    State("downstream-gmt-input", "value"),
    State("downstream-gene_mapping-input", "value"),
)
def autofill_downstream_inputs(
    workflow_id,
    selected_file_paths,
    tester_id,
    omics_type,
    project_id,
    current_counts,
    current_metadata,
    current_contrasts,
    current_expression,
    current_significant,
    current_ranked,
    current_universe,
    current_gmt,
    current_gene_mapping,
):
    if not workflow_is_downstream(workflow_id):
        raise PreventUpdate
    if workflow_id == "downstream_airway_all_atomics":
        defaults = default_airway_downstream_inputs()
        return defaults["counts"], defaults["metadata"], defaults["contrasts"], current_expression, current_significant, current_ranked, current_universe, current_gmt, current_gene_mapping
    inferred = {}
    if tester_id and project_id:
        inferred.update({key: value for key, value in downstream_previous_inputs(tester_id, omics_type or "bulk_rnaseq", project_id).items() if value})
    if selected_file_paths and tester_id and project_id:
        try:
            selected_files = selected_files_by_category(selected_file_paths or [], tester_id, omics_type or "bulk_rnaseq", project_id)
            inferred.update({key: value for key, value in infer_downstream_selected_inputs(selected_files).items() if value})
        except ValueError:
            pass
    return (
        inferred.get("counts") or current_counts,
        inferred.get("metadata") or current_metadata,
        inferred.get("contrasts") or current_contrasts,
        inferred.get("expression") or current_expression,
        inferred.get("significant") or current_significant,
        inferred.get("ranked") or current_ranked,
        inferred.get("universe") or current_universe,
        inferred.get("gmt") or current_gmt,
        inferred.get("gene_mapping") or current_gene_mapping,
    )


@app.callback(
    Output("advanced-workflow-options", "style"),
    Input("advanced-options-toggle", "value"),
)
def toggle_advanced_options(show_advanced):
    return {"display": "block"} if show_advanced else {"display": "none"}


@app.callback(
    Output("raw-qc-panel", "style"),
    Input("workflow-select", "value"),
    Input("omics-select", "value"),
)
def show_raw_qc_parameters(workflow_id, omics_type):
    if omics_type == "bulk_rnaseq" and workflow_has_raw_qc(workflow_id):
        return {"display": "block"}
    return {"display": "none"}


@app.callback(
    Output("parameter-panel", "style"),
    Input("workflow-select", "value"),
    Input("omics-select", "value"),
)
def show_parameters(workflow_id, omics_type):
    if omics_type == "bulk_rnaseq" and workflow_has_trimming(workflow_id):
        return {"display": "block"}
    return {"display": "none"}


@app.callback(
    Output("strandedness-panel", "style"),
    Input("workflow-select", "value"),
    Input("omics-select", "value"),
)
def show_strandedness_parameters(workflow_id, omics_type):
    if omics_type == "bulk_rnaseq" and workflow_has_strandedness(workflow_id):
        return {"display": "block"}
    return {"display": "none"}


@app.callback(
    Output("reference-panel", "style"),
    Input("workflow-select", "value"),
    Input("omics-select", "value"),
)
def show_reference_parameters(workflow_id, omics_type):
    if omics_type == "bulk_rnaseq" and workflow_has_reference(workflow_id):
        return {"display": "block"}
    return {"display": "none"}


@app.callback(
    Output("strandedness-input-wrapper", "style"),
    Input("workflow-select", "value"),
    Input("omics-select", "value"),
)
def show_strandedness_input_panel(workflow_id, omics_type):
    if omics_type == "bulk_rnaseq" and workflow_uses_existing_trimmed_reads(workflow_id):
        return {"display": "block"}
    return {"display": "none"}


@app.callback(
    Output("previous-trim-manifest-fields", "style"),
    Output("manual-trimmed-fastq-message", "style"),
    Input("trim-input-mode", "value"),
)
def show_trim_input_mode_fields(trim_input_mode):
    if trim_input_mode == "manual_trimmed_fastq":
        return {"display": "none"}, {"display": "block"}
    return {"display": "block"}, {"display": "none"}


@app.callback(
    Output("trim-manifest-select", "options"),
    Output("trim-manifest-select", "value"),
    Output("trim-manifest-message", "children"),
    Input("workflow-select", "value"),
    Input("tester-select", "value"),
    Input("omics-select", "value"),
    Input("project-select", "value"),
    Input("current-job", "data"),
    State("current-session", "data"),
    State("trim-manifest-select", "value"),
)
def configure_trim_manifest_options(workflow_id, tester_id, omics_type, project_id, current_job, session, current_value):
    if not workflow_uses_existing_trimmed_reads(workflow_id) or not tester_id or not project_id:
        return [], None, ""
    session_id = (session or {}).get("session_id")
    if not session_id:
        return [], None, "Browser session is not ready yet."
    options = workflow_service.trim_manifest_options(session_id, tester_id, omics_type or "bulk_rnaseq", project_id)
    values = {item["value"] for item in options}
    value = current_value if current_value in values else (options[0]["value"] if options else None)
    message = "No previous trim manifest found for this student/project yet." if not options else f"{len(options)} previous trim manifest(s) available."
    return options, value, message


@app.callback(
    Output("execution-panel", "style"),
    Input("workflow-select", "value"),
    Input("omics-select", "value"),
)
def show_execution_parameters(workflow_id, omics_type):
    if omics_type == "bulk_rnaseq" and workflow_id in steps_by_id:
        return {"display": "block"}
    return {"display": "none"}


@app.callback(
    Output("reference-bundle-fields", "style"),
    Output("custom-reference-fields", "style"),
    Input("reference-mode-dropdown", "value"),
)
def show_reference_mode_fields(reference_mode):
    if reference_mode == "prebuilt_reference":
        return {"display": "block"}, {"display": "none"}
    if reference_mode == "custom_reference":
        return {"display": "none"}, {"display": "block"}
    return {"display": "none"}, {"display": "none"}


@app.callback(
    Output("reference-status-message", "children"),
    Input("reference-mode-dropdown", "value"),
    Input("selected-route-dropdown", "value"),
)
def describe_reference_mode(reference_mode, selected_route):
    mode_label = {
        "demo_reference": "Demo reference is for workflow testing only, not biological interpretation.",
        "prebuilt_reference": "Prebuilt reference mode validates existing FASTA/GTF/index files.",
        "custom_reference": "Custom reference mode can build missing Salmon, STAR, or HISAT2 indexes from uploaded/provided reference files.",
    }.get(reference_mode, "")
    route_label = {
        "salmon": "Salmon index is required or will be built from transcriptome FASTA.",
        "star": "STAR index is required or will be built from genome FASTA plus GTF.",
        "hisat2": "HISAT2 index is required or will be built from genome FASTA.",
        "custom": "Choose the branch steps carefully; validation checks only selected downstream requirements.",
    }.get(selected_route, "")
    return f"{mode_label} {route_label}".strip()


@app.callback(
    Output("salmon-strandedness-fields", "style"),
    Output("rseqc-strandedness-fields", "style"),
    Input("strandedness-method-dropdown", "value"),
    Input("workflow-select", "value"),
)
def show_strandedness_method_fields(strandedness_method, workflow_id):
    if strandedness_method == "rseqc_validation":
        return {"display": "none"}, {"display": "block"}
    if workflow_has_reference(workflow_id):
        return {"display": "none"}, {"display": "none"}
    return {"display": "block"}, {"display": "none"}


@app.callback(
    Output("selected-route-dropdown", "value"),
    Input("workflow-select", "value"),
    State("selected-route-dropdown", "value"),
)
def sync_reference_route_with_workflow(workflow_id, current_route):
    if workflow_id in {"salmon_quant_only", "salmon_count_matrix", "qc_trim_strandedness"}:
        return "salmon"
    if workflow_id in {"star_align_only", "star_count_matrix", "star_htseq_route"}:
        return "star"
    if workflow_id in {"hisat2_align_only", "hisat2_alignment_route", "hisat2_featurecounts_route", "hisat2_htseq_route"}:
        return "hisat2"
    return current_route or "salmon"


@app.callback(
    Output("adapter-sequence-r1-input", "value"),
    Output("adapter-sequence-r2-input", "value"),
    Output("cutadapt-error-rate-input", "value"),
    Output("cutadapt-minimum-overlap-input", "value"),
    Output("adapter-default-message", "children"),
    Input("trimming-tool-dropdown", "value"),
    Input("adapter-preset-dropdown", "value"),
    State("adapter-sequence-r1-input", "value"),
    State("adapter-sequence-r2-input", "value"),
)
def populate_cutadapt_adapter_defaults(trimming_tool, adapter_preset, current_r1, current_r2):
    if trimming_tool != "cutadapt":
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update, ""
    if adapter_preset == "illumina_truseq":
        return (
            DEFAULT_CUTADAPT_ADAPTER_R1,
            DEFAULT_CUTADAPT_ADAPTER_R2,
            DEFAULT_CUTADAPT_ERROR_RATE,
            DEFAULT_CUTADAPT_MINIMUM_OVERLAP,
            dbc.Alert(ILLUMINA_TRUSEQ_HELP, color="info", className="py-1 px-2 mb-2"),
        )
    return (
        current_r1 or "",
        current_r2 or "",
        DEFAULT_CUTADAPT_ERROR_RATE,
        DEFAULT_CUTADAPT_MINIMUM_OVERLAP,
        dbc.Alert("Custom adapter mode: please provide adapter sequences.", color="warning", className="py-1 px-2 mb-2"),
    )


@app.callback(
    Output("strandedness-approval-message", "children"),
    Input("manual-strandedness-dropdown", "value"),
)
def describe_strandedness_approval(value):
    if value and value != "pending":
        return f"Manual override selected: {value}. Downstream Salmon/STAR may reuse this approved library type."
    return "Pending: downstream Salmon quantification or STAR alignment should stay locked until inference completes or you approve an override."


@app.callback(
    Output("run-analysis", "disabled"),
    Input("tester-select", "value"),
    Input("omics-select", "value"),
    Input("workflow-select", "value"),
)
def toggle_run_button(tester_id, omics_type, workflow_id):
    return not (
        tester_id
        and omics_type == "bulk_rnaseq"
        and workflow_id in steps_by_id
        and omics_types.get(omics_type, {}).get("enabled")
    )


@app.callback(
    Output("upload-message", "children"),
    Input("fastq-upload", "contents"),
    Input("metadata-upload", "contents"),
    Input("reference-upload", "contents"),
    State("fastq-upload", "filename"),
    State("metadata-upload", "filename"),
    State("reference-upload", "filename"),
    State("current-session", "data"),
    State("tester-select", "value"),
    State("omics-select", "value"),
    State("project-select", "value"),
    prevent_initial_call=True,
)
def save_dash_upload(
    fastq_contents,
    metadata_contents,
    reference_contents,
    fastq_names,
    metadata_name,
    reference_names,
    session,
    tester_id,
    omics_type,
    project_id,
):
    session_id = (session or {}).get("session_id")

    if not session_id:
        return dbc.Alert(
            "Browser session is not ready yet. Refresh the page and try again.",
            color="warning",
        )
    if not tester_id:
        return dbc.Alert("Select a student/tester before uploading files.", color="warning")
    if not project_id:
        return dbc.Alert("Select or create a project before uploading files.", color="warning")
    if (omics_type or "bulk_rnaseq") != "bulk_rnaseq":
        return dbc.Alert("Uploads for this omics type are under development.", color="warning")
    job_store.upsert_project(tester_id, omics_type or "bulk_rnaseq", safe_slug(project_id), project_label_from_id(project_id))

    ctx = dash.callback_context

    if not ctx.triggered:
        raise PreventUpdate

    triggered = ctx.triggered[0]["prop_id"].split(".")[0]

    if triggered == "fastq-upload":
        contents = fastq_contents or []
        names = fastq_names or []

        if isinstance(contents, str):
            contents = [contents]
            names = [names]

        upload_type = "fastq"

    elif triggered == "metadata-upload":
        contents = (
            [metadata_contents]
            if isinstance(metadata_contents, str)
            else (metadata_contents or [])
        )

        names = (
            [metadata_name]
            if isinstance(metadata_name, str)
            else (metadata_name or [])
        )

        upload_type = "metadata"
    else:
        contents = reference_contents or []
        names = reference_names or []
        if isinstance(contents, str):
            contents = [contents]
            names = [names]
        upload_type = "reference"

    saved = []

    try:
        for content, filename in zip(contents, names):
            saved.append(
                upload_service.save_bytes(
                    session_id,
                    upload_type,
                    filename,
                    decode_dash_upload(content),
                    tester_id,
                    omics_type,
                    project_id,
                )
            )
            job_store.upsert_upload_record({
                "project_id": safe_slug(project_id),
                "tester_id": tester_id,
                "omics_type": omics_type,
                "session_id": session_id,
                "category": upload_type,
                "filename": saved[-1]["filename"],
                "stored_path": saved[-1]["path"],
                "size_bytes": saved[-1]["size"],
            })

    except Exception as exc:
        return dbc.Alert(str(exc), color="danger")

    return dbc.Alert(f"Uploaded {len(saved)} {upload_type} file(s).", color="success")


@app.callback(
    Output("uploaded-files-list", "children"),
    Output("selected-files-checklist", "options"),
    Output("selected-files-checklist", "value"),
    Input("poller", "n_intervals"),
    Input("upload-message", "children"),
    Input("file-search", "value"),
    Input("file-category-filter", "value"),
    State("current-session", "data"),
    State("tester-select", "value"),
    State("omics-select", "value"),
    State("project-select", "value"),
    State("selected-files-checklist", "value"),
)
def show_uploaded_files(_, __, search, category_filter, session, tester_id, omics_type, project_id, current_selection):
    session_id = (session or {}).get("session_id")

    if not session_id:
        return html.Div("Session is loading...", className="text-muted"), [], []
    if not tester_id:
        return html.Div("Select a student/tester to view workspace uploads.", className="text-muted"), [], []
    if not project_id:
        return html.Div("Select or create a project to view files.", className="text-muted"), [], []

    return file_list_component(
        tester_id,
        omics_type or "bulk_rnaseq",
        project_id,
        search,
        category_filter or "all",
        current_selection,
    )


@app.callback(
    Output("reference-mode-dropdown", "value"),
    Output("reference-bundle-input", "value"),
    Output("genome-fasta-input", "value"),
    Output("gtf-input", "value"),
    Output("transcriptome-fasta-input", "value"),
    Output("tx2gene-input", "value"),
    Input("upload-message", "children"),
    Input("project-select", "value"),
    State("tester-select", "value"),
    State("omics-select", "value"),
    State("reference-mode-dropdown", "value"),
    State("reference-bundle-input", "value"),
    State("genome-fasta-input", "value"),
    State("gtf-input", "value"),
    State("transcriptome-fasta-input", "value"),
    State("tx2gene-input", "value"),
)
def autofill_reference_inputs(
    _,
    project_id,
    tester_id,
    omics_type,
    current_mode,
    current_bundle,
    current_genome,
    current_gtf,
    current_transcriptome,
    current_tx2gene,
):
    if not tester_id or not project_id:
        raise PreventUpdate
    inferred = infer_reference_inputs_from_uploads(tester_id, omics_type or "bulk_rnaseq", project_id)
    if not any(inferred.values()):
        raise PreventUpdate
    return (
        "custom_reference" if any(inferred.get(key) for key in ("genome_fasta", "gtf", "transcriptome_fasta", "tx2gene")) else current_mode,
        current_bundle or inferred["reference_bundle"],
        current_genome or inferred["genome_fasta"],
        current_gtf or inferred["gtf"],
        current_transcriptome or inferred["transcriptome_fasta"],
        current_tx2gene or inferred["tx2gene"],
    )


@app.callback(
    Output("selected-step-details", "children"),
    Input("workflow-select", "value"),
    Input("omics-select", "value"),
    Input("tester-select", "value"),
    Input("project-select", "value"),
)
def show_step(workflow_id, omics_type, tester_id, project_id):
    if omics_type != "bulk_rnaseq":
        return [
            html.H5(omics_label(omics_type)),
            html.P(omics_types.get(omics_type, {}).get("message", "This omics module is under development.")),
        ]
    step = steps_by_id.get(workflow_id or "qc_trim", steps_by_id["qc_trim"])

    return [
        html.H5(step.get("label", step["name"])),
        html.P([html.Strong("Student: "), tester_label(tester_id) if tester_id else "Not selected"]),
        html.P([html.Strong("Omics type: "), omics_label(omics_type)]),
        html.P([html.Strong("Project: "), project_label_from_id(project_id)]),
        html.P(step["description"]),
        html.Code(", ".join(step.get("selected_steps", []))),
    ]


@app.callback(
    Output("current-job", "data"),
    Output("run-message", "children"),
    Input("run-analysis", "n_clicks"),
    State("quality-threshold-dropdown", "value"),
    State("minimum-read-length-dropdown", "value"),
    State("trim-poly-g-dropdown", "value"),
    State("trimming-tool-dropdown", "value"),
    State("trim-poly-x-dropdown", "value"),
    State("adapter-preset-dropdown", "value"),
    State("adapter-sequence-r1-input", "value"),
    State("adapter-sequence-r2-input", "value"),
    State("adapter-fasta-input", "value"),
    State("trim-front-r1-input", "value"),
    State("trim-front-r2-input", "value"),
    State("cutadapt-error-rate-input", "value"),
    State("cutadapt-minimum-overlap-input", "value"),
    State("reference-mode-dropdown", "value"),
    State("reference-bundle-input", "value"),
    State("selected-route-dropdown", "value"),
    State("organism-input", "value"),
    State("genome-build-input", "value"),
    State("genome-fasta-input", "value"),
    State("gtf-input", "value"),
    State("transcriptome-fasta-input", "value"),
    State("tx2gene-input", "value"),
    State("strandedness-method-dropdown", "value"),
    State("salmon-index-input", "value"),
    State("strandedness-inference-reads-input", "value"),
    State("existing-bam-input", "value"),
    State("star-index-input", "value"),
    State("hisat2-index-input", "value"),
    State("rseqc-ref-bed-input", "value"),
    State("rseqc-stranded-threshold-input", "value"),
    State("manual-strandedness-dropdown", "value"),
    State("execution-profile-dropdown", "value"),
    State("trim-input-mode", "value"),
    State("trim-manifest-select", "value"),
    State("downstream-counts-input", "value"),
    State("downstream-metadata-input", "value"),
    State("downstream-contrasts-input", "value"),
    State("downstream-expression-input", "value"),
    State("downstream-significant-input", "value"),
    State("downstream-ranked-input", "value"),
    State("downstream-universe-input", "value"),
    State("downstream-gmt-input", "value"),
    State("downstream-gene_mapping-input", "value"),
    State("current-session", "data"),
    State("tester-select", "value"),
    State("omics-select", "value"),
    State("project-select", "value"),
    State("workflow-select", "value"),
    State("selected-files-checklist", "value"),
    State("current-job", "data"),
    prevent_initial_call=True,
)
def run_selected_step(
    run_clicks,
    quality_threshold,
    minimum_read_length,
    trim_poly_g,
    trimming_tool,
    trim_poly_x,
    adapter_preset,
    adapter_sequence_r1,
    adapter_sequence_r2,
    adapter_fasta,
    trim_front_r1,
    trim_front_r2,
    cutadapt_error_rate,
    cutadapt_minimum_overlap,
    reference_mode,
    reference_bundle,
    selected_route,
    organism,
    genome_build,
    genome_fasta,
    gtf,
    transcriptome_fasta,
    tx2gene,
    strandedness_method,
    salmon_index,
    strandedness_inference_reads,
    existing_bam,
    star_index,
    hisat2_index,
    rseqc_ref_bed,
    rseqc_stranded_threshold,
    manual_strandedness,
    execution_profile,
    trim_input_mode,
    trim_manifest,
    downstream_counts,
    downstream_metadata,
    downstream_contrasts,
    downstream_expression,
    downstream_significant,
    downstream_ranked,
    downstream_universe,
    downstream_gmt,
    downstream_gene_mapping,
    session,
    tester_id,
    omics_type,
    project_id,
    workflow_id,
    selected_file_paths,
    current_job,
):
    if not run_clicks:
        return current_job, ""

    session_id = (session or {}).get("session_id")
    omics_type = omics_type or "bulk_rnaseq"

    setup_errors = []
    if not session_id:
        setup_errors.append("Browser session is not ready yet.")
    if not tester_id:
        setup_errors.append("Select a student/tester.")
    if not project_id:
        setup_errors.append("Select or create a project.")
    if omics_type != "bulk_rnaseq":
        setup_errors.append("Only Bulk RNA-seq workflows are active right now.")
    if workflow_id not in steps_by_id:
        setup_errors.append("Select an active workflow.")
    if setup_errors:
        return current_job, dbc.Alert(html.Ul([html.Li(error) for error in setup_errors], className="mb-0"), color="danger")

    try:
        selected_files = selected_files_by_category(selected_file_paths or [], tester_id, omics_type, project_id)
    except ValueError as exc:
        return current_job, dbc.Alert(str(exc), color="danger")

    if workflow_is_downstream(workflow_id):
        if workflow_id == "downstream_airway_all_atomics":
            defaults = default_airway_downstream_inputs()
            downstream_counts = defaults["counts"]
            downstream_metadata = defaults["metadata"]
            downstream_contrasts = defaults["contrasts"]
        downstream_inputs = {
            "counts": normalize_empty(downstream_counts),
            "metadata": normalize_empty(downstream_metadata),
            "contrasts": normalize_empty(downstream_contrasts),
            "expression": normalize_empty(downstream_expression),
            "significant": normalize_empty(downstream_significant),
            "ranked": normalize_empty(downstream_ranked),
            "universe": normalize_empty(downstream_universe),
            "gmt": normalize_empty(downstream_gmt),
            "gene_mapping": normalize_empty(downstream_gene_mapping),
        }
        missing = []
        for label, key in DOWNSTREAM_REQUIRED_INPUTS.get(workflow_id, []):
            value = downstream_inputs.get(key)
            if not value or not Path(str(value)).exists():
                missing.append(f"Provide a valid {label} file path.")
        if missing:
            return current_job, dbc.Alert(html.Ul([html.Li(error) for error in missing], className="mb-0"), color="danger")
        try:
            run = start_downstream_dash_job(
                workflow_id,
                tester_id,
                omics_type,
                project_id,
                session_id,
                {key: str(Path(str(value)).resolve()) for key, value in downstream_inputs.items() if value},
            )
        except Exception as exc:
            return current_job, dbc.Alert(str(exc), color="danger")
        return run, dbc.Alert(
            [
                html.Div(f"Started downstream CSV-output run {run['run_id']}."),
                html.Div("No plots are generated; download CSV/JSON outputs from the Outputs panel.", className="small mt-1"),
            ],
            color="success",
        )

    is_strandedness_only = workflow_is_strandedness_only(workflow_id)
    uses_existing_trimmed_reads = workflow_uses_existing_trimmed_reads(workflow_id)
    trim_input_mode = trim_input_mode or "previous_manifest"
    if uses_existing_trimmed_reads and trim_input_mode == "previous_manifest":
        if not trim_manifest:
            return current_job, dbc.Alert("Select a previous trim manifest for this atomic workflow.", color="danger")
        manifest_options = workflow_service.trim_manifest_options(session_id, tester_id, omics_type, project_id)
        if trim_manifest not in {item["value"] for item in manifest_options}:
            return current_job, dbc.Alert("Selected trim manifest is not from this student/project workspace.", color="danger")
        manifest_errors = workflow_service.validate_trim_manifest(trim_manifest)
        if manifest_errors:
            return current_job, dbc.Alert(html.Ul([html.Li(error) for error in manifest_errors], className="mb-0"), color="danger")
    else:
        if workflow_needs_fastq_input(workflow_id) and not selected_files["fastq"]:
            return current_job, dbc.Alert("Select at least one FASTQ file for this workflow.", color="danger")
        if uses_existing_trimmed_reads:
            pair_errors = workflow_service.validate_fastq_pairs(selected_files["fastq"])
            if pair_errors:
                return current_job, dbc.Alert(html.Ul([html.Li(error) for error in pair_errors], className="mb-0"), color="danger")

    params = {
        "quality_threshold": quality_threshold,
        "minimum_read_length": minimum_read_length,
        "trim_poly_g": trim_poly_g,
        "trimming_tool": trimming_tool,
        "trim_poly_x": trim_poly_x,
        "adapter_preset": adapter_preset,
        "adapter_sequence_r1": adapter_sequence_r1,
        "adapter_sequence_r2": adapter_sequence_r2,
        "adapter_fasta": adapter_fasta,
        "trim_front_r1": trim_front_r1,
        "trim_front_r2": trim_front_r2,
        "cutadapt_error_rate": cutadapt_error_rate,
        "cutadapt_minimum_overlap": cutadapt_minimum_overlap,
        "reference_mode": reference_mode,
        "reference_bundle": reference_bundle,
        "selected_route": selected_route,
        "organism": organism,
        "genome_build": genome_build,
        "genome_fasta": genome_fasta,
        "gtf": gtf,
        "transcriptome_fasta": transcriptome_fasta,
        "tx2gene": tx2gene,
        "strandedness_method": strandedness_method,
        "salmon_index": salmon_index,
        "strandedness_inference_reads": strandedness_inference_reads,
        "existing_bam": existing_bam,
        "star_index": star_index,
        "hisat2_index": hisat2_index,
        "rseqc_ref_bed": rseqc_ref_bed,
        "rseqc_stranded_threshold": rseqc_stranded_threshold,
        "manual_strandedness": None if manual_strandedness == "pending" else manual_strandedness,
        "strandedness_approval_status": "overridden" if manual_strandedness not in {None, "pending"} else "pending",
        "execution_profile": execution_profile or "local",
        "trim_input_mode": trim_input_mode,
        "trim_manifest": trim_manifest if uses_existing_trimmed_reads and trim_input_mode == "previous_manifest" else None,
    }
    params = {key: normalize_empty(value) for key, value in params.items()}
    params = resolve_cutadapt_defaults(params)
    params = resolve_reference_defaults(params)
    execution_profile, profile_message = resolve_execution_profile(params)
    params["execution_profile"] = execution_profile
    params["workflow_id"] = workflow_id
    params["omics_type"] = omics_type
    params["tester_id"] = tester_id
    params["tester_label"] = tester_label(tester_id)
    params["project_id"] = project_id
    params["project_label"] = project_label_from_id(project_id)

    errors = validate_run_params(workflow_id, params, session_id, tester_id, omics_type, selected_files)
    if errors:
        return current_job, dbc.Alert(html.Ul([html.Li(error) for error in errors], className="mb-0"), color="danger")

    try:
        run = workflow_service.create_run(
            session_id,
            workflow_id,
            params,
            tester_id=tester_id,
            tester_label=tester_label(tester_id),
            omics_type=omics_type,
            project_id=project_id,
            project_label=project_label_from_id(project_id),
            selected_files=selected_files,
        )
        workflow_service.start_job(run["job_id"])
    except Exception as exc:
        return current_job, dbc.Alert(str(exc), color="danger")

    message = [html.Div(f"Started {workflow_id} run {run['run_id']}.")]
    if profile_message:
        message.append(html.Div(profile_message, className="small mt-1"))
    return run, dbc.Alert(message, color="success")


@app.callback(
    Output("job-status-panel", "children"),
    Output("logs-panel", "children"),
    Output("outputs-panel", "children"),
    Input("poller", "n_intervals"),
    State("current-job", "data"),
    State("current-session", "data"),
    State("tester-select", "value"),
    State("omics-select", "value"),
    State("project-select", "value"),
)
def refresh_job(_, current_job, session, tester_id, omics_type, project_id):
    session_id = (session or {}).get("session_id")
    omics_type = omics_type or "bulk_rnaseq"
    if not current_job:
        jobs = (
            job_store.list_workspace_jobs(session_id, tester_id, omics_type, project_id=project_id, limit=10)
            if session_id and tester_id
            else []
        )

        return (
            html.Div(
                [
                    html.P("No active job."),
                    html.H6("Recent jobs in this workspace"),
                    html.Pre(json.dumps(jobs, indent=2)[:2000]),
                ]
            ),
            "",
            "No outputs yet.",
        )

    job = workflow_service.refresh_job(current_job["job_id"])

    if not job:
        return "Job not found.", "", "No outputs."
    if workflow_is_downstream(job.get("workflow_id")):
        job = refresh_downstream_job(job)
    if tester_id and (
        job.get("tester_id") != tester_id
        or job.get("omics_type") != omics_type
        or job.get("session_id") != session_id
        or (project_id and job.get("project_id") != project_id)
    ):
        return "This job belongs to a different workspace.", "", "No outputs."

    stdout = (
        Path(job["stdout_path"]).read_text(errors="replace")[-4000:]
        if job.get("stdout_path") and Path(job["stdout_path"]).exists()
        else ""
    )

    stderr = (
        Path(job["stderr_path"]).read_text(errors="replace")[-4000:]
        if job.get("stderr_path") and Path(job["stderr_path"]).exists()
        else ""
    )

    record = (
        downstream_execution_record(job)
        if workflow_is_downstream(job.get("workflow_id"))
        else workflow_service.execution_record(job["job_id"], refresh=False)
    )
    files = record.get("all_result_files") or workflow_service.output_files(job["job_id"])

    output_links = [
        html.Li(
            html.A(
                Path(path).name,
                href=(
                    f"{API_PREFIX}/file?path={path}"
                    f"&session_id={session_id}&tester_id={tester_id}&omics_type={omics_type}&project_id={project_id}"
                ),
                target="_blank",
            )
        )
        for path in files
    ] or [html.Li("No output files discovered yet.")]

    job_status = str(job.get("status", "unknown")).lower()

    terminal_statuses = {
        "succeeded",
        "success",
        "completed",
        "complete",
        "done",
        "failed",
        "error",
        "cancelled",
        "canceled",
    }

    is_running = job_status not in terminal_statuses

    status_color = {
        "created": "secondary",
        "queued": "secondary",
        "pending": "secondary",
        "submitted": "secondary",
        "starting": "info",
        "started": "info",
        "running": "primary",
        "succeeded": "success",
        "success": "success",
        "completed": "success",
        "complete": "success",
        "done": "success",
        "failed": "danger",
        "error": "danger",
        "cancelled": "warning",
        "canceled": "warning",
    }.get(job_status, "primary" if is_running else "secondary")

    status_badge = dbc.Badge(job_status, color=status_color)
    metadata = json.loads(job.get("metadata_json") or "{}")
    submitted_params = metadata.get("submitted_params", {})
    selected_files = metadata.get("selected_files", {})
    try:
        job_selected_steps = workflow_service.selected_steps_for_workflow(job.get("workflow_id"))
    except Exception:
        job_selected_steps = []
    failed_steps = [
        step for step in record.get("step_outputs", [])
        if str(step.get("status", "")).lower() == "failed"
    ]
    error_summary = None
    if failed_steps:
        error_items = []
        for step in failed_steps:
            task = (step.get("tasks") or [{}])[0]
            tail = task.get("stderr_tail") or task.get("log_tail") or "No task stderr captured."
            error_items.append(
                html.Li(
                    [
                        html.Strong(step.get("process_name")),
                        f": {tail[-500:]}",
                    ]
                )
            )
        error_summary = dbc.Alert(
            [html.Strong("Failed process summary"), html.Ul(error_items, className="mb-0")],
            color="danger",
        )

    running_loader = (
        dbc.Alert(
            [
                dbc.Spinner(
                    size="sm",
                    color="primary",
                    spinner_class_name="me-2",
                ),
                html.Span("Workflow is running. Logs and outputs refresh automatically."),
            ],
            color="info",
            className="d-flex align-items-center",
        )
        if is_running
        else None
    )

    return (
        html.Div(
            [
                running_loader,
                html.P([html.Strong("Job: "), job["job_id"], " ", status_badge]),
                html.P([html.Strong("Run: "), job["run_id"]]),
                html.P([html.Strong("Student: "), job.get("tester_label") or tester_label(job.get("tester_id"))]),
                html.P([html.Strong("Omics type: "), omics_label(job.get("omics_type"))]),
                html.P([html.Strong("Project: "), job.get("project_label") or job.get("project_id")]),
                html.P([html.Strong("Workflow: "), steps_by_id.get(job.get("workflow_id"), {}).get("label", job.get("workflow_id"))]),
                error_summary,
                workflow_plan_component(record, session_id, tester_id, omics_type, project_id, API_PREFIX),
                run_parameter_summary(submitted_params, job_selected_steps),
                strandedness_result_summary(files)
                if submitted_params.get("strandedness_method") or any(path.endswith(".strandedness_call.json") for path in files)
                else None,
                selected_file_summary(selected_files),
            ]
        ),
        html.Div(
            [
                log_block("STDOUT", stdout),
                log_block("STDERR", stderr),
            ]
        ),
        html.Div([html.H6("Output files"), html.Ul(output_links)]),
    )   


def audit_dash_callbacks(dash_app):
    print("\n=== REGISTERED DASH CALLBACKS ===")
    for output, meta in dash_app.callback_map.items():
        inputs = meta.get("inputs", [])
        state = meta.get("state", [])
        print(f"\nOUTPUT: {output}")
        print(f"  inputs={len(inputs)} state={len(state)}")
        for item in inputs:
            print(f"    Input: {item.get('id')}.{item.get('property')}")
        for item in state:
            print(f"    State: {item.get('id')}.{item.get('property')}")
        print(f"  callback: {getattr(meta.get('callback'), '__name__', meta.get('callback'))}")
    print("=== END REGISTERED DASH CALLBACKS ===\n")


if os.environ.get("SURVOM_DEBUG_CALLBACKS") == "1":
    audit_dash_callbacks(app)


if __name__ == "__main__":
    host = os.environ.get("SURVOM_DASH_HOST", "127.0.0.1")
    port = int(os.environ.get("SURVOM_DASH_PORT", "8057"))
    app.run(debug=False, host=host, port=port)
