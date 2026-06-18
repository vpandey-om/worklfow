from __future__ import annotations

import base64
import json
import os
import shutil
from pathlib import Path

import dash
from dash import Input, Output, State, dcc, html
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from flask import jsonify, request, send_file

from services.job_store import JobStore
from services.upload_service import UploadService, safe_slug
from services.workflow_service import WorkflowService


APP_ROOT = Path(__file__).resolve().parent
PIPELINE_ROOT = APP_ROOT.parent / "survom-pipelines"
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

URL_PREFIX = "/upload/"
API_PREFIX = "/upload/api"

upload_service = UploadService(APP_ROOT / "uploads", MAX_UPLOAD_BYTES)
job_store = JobStore(APP_ROOT / "runs" / "jobs.sqlite")
workflow_service = WorkflowService(APP_ROOT, PIPELINE_ROOT, upload_service, job_store)

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
    return [
        {"label": steps_by_id[workflow_id].get("label", steps_by_id[workflow_id]["name"]), "value": workflow_id}
        for workflow_id in config.get("workflows", [])
        if workflow_id in steps_by_id
    ]


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


def trimming_parameter_panel():
    return html.Div(
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
        id="trimming-parameter-panel",
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
                upload_box("metadata-upload", "Select or drop sample sheet here", ".csv, .tsv, .xlsx", False),
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

        html.H6("Workflow parameters"),
        html.Div(trimming_parameter_panel(), id="parameter-panel", className="mb-3"),

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


def normalize_empty(value):
    return None if value == "" else value


def uses_container_profile(profile: str | None) -> bool:
    return bool({token.strip() for token in (profile or "").split(",")} & {"docker", "singularity", "aws"})


def docker_available() -> bool:
    return shutil.which("docker") is not None


def resolve_execution_profile(params: dict) -> tuple[str, str | None]:
    requested = params.get("execution_profile") or "local"
    if params.get("trimming_tool") != "cutadapt":
        return requested, None
    if shutil.which("cutadapt") or uses_container_profile(requested):
        return requested, None
    if docker_available():
        return "local_docker", CUTADAPT_DOCKER_PROFILE_MESSAGE
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


def validate_run_params(
    workflow_id: str,
    params: dict,
    session_id: str,
    tester_id: str,
    omics_type: str,
    selected_files: dict[str, list[str]] | None = None,
) -> list[str]:
    if not workflow_has_trimming(workflow_id):
        return []
    params = resolve_cutadapt_defaults(params)

    errors = []
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

    if params.get("trimming_tool") == "cutadapt":
        if not shutil.which("cutadapt") and not uses_container_profile(params.get("execution_profile")) and not docker_available():
            errors.append(CUTADAPT_NOT_AVAILABLE_MESSAGE)
        if not params.get("adapter_sequence_r1"):
            errors.append("Cutadapt requires an R1 adapter sequence.")
        selected_fastqs = [Path(path).name for path in (selected_files or {}).get("fastq", [])]
        if not selected_fastqs:
            uploads = upload_service.list_uploads(session_id, tester_id, omics_type)
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

    return errors


def run_parameter_summary(params: dict):
    if not params:
        return None
    poly_x = "on" if params.get("trim_poly_x") else "off"
    return html.Div(
        [
            html.H6("Run parameters", className="mt-3"),
            html.Ul(
                [
                    html.Li(f"Trimming tool: {params.get('trimming_tool', 'fastp')}"),
                    html.Li(f"Minimum base quality: Q{params.get('quality_threshold', 20)}"),
                    html.Li(f"Minimum read length: {params.get('minimum_read_length', 20)} bp"),
                    html.Li(f"PolyG trimming: {params.get('trim_poly_g', 'auto')}"),
                    html.Li(f"PolyX trimming: {poly_x}"),
                    html.Li(f"Adapter preset: {params.get('adapter_preset') or 'not used'}"),
                    html.Li(f"Adapter source: {params.get('adapter_source') or 'not used'}"),
                ],
                className="small ps-3 mb-0",
            ),
        ]
    )


def selected_file_summary(selected_files: dict | None):
    if not selected_files:
        return None
    rows = []
    for category in ("fastq", "metadata", "other"):
        for path in selected_files.get(category, []):
            rows.append(html.Li(f"{category}: {Path(path).name}"))
    if not rows:
        return None
    return html.Div([html.H6("Selected files", className="mt-3"), html.Ul(rows, className="small ps-3 mb-0")])


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
    Output("parameter-panel", "style"),
    Input("workflow-select", "value"),
    Input("omics-select", "value"),
)
def show_parameters(workflow_id, omics_type):
    if omics_type == "bulk_rnaseq" and workflow_has_trimming(workflow_id):
        return {"display": "block"}
    return {"display": "none"}


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
    State("fastq-upload", "filename"),
    State("metadata-upload", "filename"),
    State("current-session", "data"),
    State("tester-select", "value"),
    State("omics-select", "value"),
    State("project-select", "value"),
    prevent_initial_call=True,
)
def save_dash_upload(
    fastq_contents,
    metadata_contents,
    fastq_names,
    metadata_name,
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

    else:
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
    if not selected_files["fastq"]:
        return current_job, dbc.Alert("Select at least one FASTQ file for bulk RNA-seq workflows.", color="danger")

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
    }
    params = {key: normalize_empty(value) for key, value in params.items()}
    params = resolve_cutadapt_defaults(params)
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

    files = workflow_service.output_files(job["job_id"])

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

    running_loader = (
        dbc.Alert(
            [
                dbc.Spinner(
                    size="sm",
                    color="primary",
                    className="me-2",
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
                run_parameter_summary(submitted_params),
                selected_file_summary(selected_files),
            ]
        ),
        f"STDOUT\n{stdout}\n\nSTDERR\n{stderr}",
        html.Ul(output_links),
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
    app.run(debug=False, host="127.0.0.1", port=8056)
