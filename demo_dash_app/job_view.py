from __future__ import annotations

import json
from pathlib import Path

from dash import dcc, html
import dash_bootstrap_components as dbc


def process_status_badge(status: str | None):
    value = (status or "pending").lower()
    color = {
        "pending": "secondary",
        "running": "primary",
        "completed": "success",
        "succeeded": "success",
        "failed": "danger",
        "error": "danger",
    }.get(value, "secondary")
    return dbc.Badge(value, color=color)


def _artifact_link(
    path: str,
    label: str,
    api_prefix: str,
    session_id: str,
    tester_id: str,
    omics_type: str,
    project_id: str,
):
    if not path:
        return None
    return html.A(
        label,
        href=(
            f"{api_prefix}/file?path={path}"
            f"&session_id={session_id}&tester_id={tester_id}&omics_type={omics_type}&project_id={project_id}"
        ),
        target="_blank",
        className="me-2",
    )


def workflow_plan_component(
    record: dict,
    session_id: str,
    tester_id: str,
    omics_type: str,
    project_id: str,
    api_prefix: str,
):
    steps = record.get("step_outputs") or []
    if not steps:
        return html.Div("Workflow steps will appear after the run starts.", className="text-muted")

    visible_steps = [
        step
        for step in steps
        if str(step.get("status", "")).lower() != "pending" or step.get("tasks") or step.get("missing_outputs")
    ]
    if not visible_steps:
        visible_steps = steps[:1]

    rows = []
    for step in visible_steps:
        tasks = step.get("tasks") or []
        first_task = tasks[0] if tasks else {}
        missing = step.get("missing_outputs") or []
        links = [
            _artifact_link(first_task.get("stdout_path"), "stdout", api_prefix, session_id, tester_id, omics_type, project_id),
            _artifact_link(first_task.get("stderr_path"), "stderr", api_prefix, session_id, tester_id, omics_type, project_id),
            _artifact_link(first_task.get("command_log_path"), "log", api_prefix, session_id, tester_id, omics_type, project_id),
        ]
        links = [link for link in links if link]
        outputs = sum(output.get("count", 0) for output in step.get("outputs", []))
        missing_text = ""
        if missing:
            missing_text = "Missing: " + ", ".join(item.get("name") or item.get("pattern", "") for item in missing)
        rows.append(
            html.Tr(
                [
                    html.Td(step.get("step_name") or step.get("process_name")),
                    html.Td(process_status_badge(step.get("status"))),
                    html.Td(str(outputs)),
                    html.Td(missing_text, className="small text-danger"),
                    html.Td(links),
                ]
            )
        )

    hidden_count = max(0, len(steps) - len(visible_steps))
    hidden_note = (
        html.Div(f"{hidden_count} pending step(s) hidden until they start.", className="small text-muted mt-1")
        if hidden_count
        else None
    )
    return html.Div(
        [
            html.H6("Workflow steps"),
            dbc.Table(
                [
                    html.Thead(
                        html.Tr(
                            [
                                html.Th("Step"),
                                html.Th("Status"),
                                html.Th("Outputs"),
                                html.Th("Issues"),
                                html.Th("Logs"),
                            ]
                        )
                    ),
                    html.Tbody(rows),
                ],
                bordered=True,
                hover=True,
                responsive=True,
                size="sm",
            ),
            hidden_note,
        ]
    )


def log_block(title: str, content: str):
    text = content or ""
    return html.Div(
        [
            html.Div(
                [
                    html.Strong(title),
                    dcc.Clipboard(content=text, title=f"Copy {title}", className="ms-2"),
                ],
                className="d-flex align-items-center mb-1",
            ),
            html.Pre(text or "(empty)", className="small bg-light border rounded p-2"),
        ],
        className="mb-3",
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


def run_parameter_summary(params: dict, workflow_steps: list[str] | None = None):
    if not params:
        return None
    workflow_steps = set(workflow_steps or [])
    items = []
    if "rnaseq_04_adapter_quality_trimming" in workflow_steps:
        poly_x = "on" if params.get("trim_poly_x") else "off"
        items.extend(
            [
                html.Li(f"Trimming tool: {params.get('trimming_tool', 'fastp')}"),
                html.Li(f"Minimum base quality: Q{params.get('quality_threshold', 20)}"),
                html.Li(f"Minimum read length: {params.get('minimum_read_length', 20)} bp"),
                html.Li(f"PolyG trimming: {params.get('trim_poly_g', 'auto')}"),
                html.Li(f"PolyX trimming: {poly_x}"),
            ]
        )
    if "rnaseq_06a_reference_build_validation" in workflow_steps:
        items.extend(
            [
                html.Li(f"Reference mode: {params.get('reference_mode')}"),
                html.Li(f"Selected route: {params.get('selected_route') or 'salmon'}"),
                html.Li(f"Organism/build: {params.get('organism') or 'not set'} / {params.get('genome_build') or 'not set'}"),
            ]
        )
    if "rnaseq_06_strandedness_inference" in workflow_steps:
        items.append(html.Li(f"Strandedness method: {params.get('strandedness_method')}"))
        if params.get("manual_strandedness"):
            items.append(html.Li(f"Manual strandedness: {params.get('manual_strandedness')}"))
    if params.get("trim_manifest"):
        items.append(html.Li(f"Trim manifest: {Path(str(params.get('trim_manifest'))).name}"))
    if not items:
        return None
    return html.Div(
        [
            html.H6("Run parameters", className="mt-3"),
            html.Ul(items, className="small ps-3 mb-0"),
        ]
    )


def strandedness_result_summary(files: list[str]):
    calls = [Path(path) for path in files if path.endswith(".strandedness_call.json")]
    if not calls:
        return None

    items = []
    for call_path in calls:
        try:
            call = json.loads(call_path.read_text())
        except Exception as exc:
            items.append(html.Li(f"{call_path.name}: failed to read call ({exc})"))
            continue
        raw_status = call.get("status", "pending")
        display_status = "ambiguous" if raw_status == "needs_review" else raw_status
        color = {
            "pending": "secondary",
            "running": "primary",
            "completed": "success",
            "ambiguous": "warning",
            "failed": "danger",
        }.get(display_status, "secondary")
        inferred = call.get("inferred_library_type") or "not inferred"
        items.append(
            html.Li(
                [
                    html.Strong(call.get("sample_id", call_path.stem)),
                    ": ",
                    dbc.Badge(display_status, color=color, className="me-2"),
                    f"inferred library type = {inferred}. ",
                    html.Span(call.get("reason", ""), className="text-muted"),
                ]
            )
        )
    return html.Div(
        [
            html.H6("Strandedness Inference", className="mt-3"),
            html.Ul(items, className="small ps-3 mb-1"),
        ]
    )
