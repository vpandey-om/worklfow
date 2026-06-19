#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PIPELINE_ROOT = ROOT / "survom-pipelines"
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(PIPELINE_ROOT / "scripts"))

from builder.compiler import NextflowCompiler
from builder.models import RunRequest
from builder.outputs import collect_run_outputs
from builder.registry import load_steps
from smoke_workflows import CASES, Case, reference_params


RUNS_ROOT = ROOT / "runs"
CASE_BY_NAME = {case.name: case for case in CASES}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def slug(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._-")
    return cleaned or "workflow"


def parse_param_value(value: str):
    lowered = value.strip().lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    if lowered in {"none", "null"}:
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def parse_param_overrides(items: list[str] | None) -> dict:
    overrides = {}
    for item in items or []:
        if "=" not in item:
            raise ValueError(f"Invalid --param '{item}'. Use --param key=value.")
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise ValueError(f"Invalid --param '{item}'. Parameter name is empty.")
        overrides[key] = parse_param_value(value)
    return overrides


def new_run_id(workflow: str) -> str:
    return f"{slug(workflow)}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"


def run_dir(run_id: str) -> Path:
    return RUNS_ROOT / run_id


def state_path(run_id: str) -> Path:
    return run_dir(run_id) / "state.json"


def load_state(run_id: str) -> dict:
    path = state_path(run_id)
    if not path.exists():
        raise FileNotFoundError(f"Cannot resume '{run_id}': missing {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def save_state(state: dict) -> None:
    path = state_path(state["run_id"])
    path.parent.mkdir(parents=True, exist_ok=True)
    state["updated_at"] = utc_now()
    path.write_text(json.dumps(state, indent=2), encoding="utf-8")


def initial_state(run_id: str, workflow: str, output_dir: Path) -> dict:
    now = utc_now()
    return {
        "run_id": run_id,
        "workflow": workflow,
        "status": "created",
        "created_at": now,
        "updated_at": now,
        "current_step": None,
        "failed_step": None,
        "completed_steps": [],
        "steps": [],
        "input_dir": str(run_dir(run_id) / "inputs"),
        "output_dir": str(output_dir),
        "run_folder": str(run_dir(run_id)),
        "compiled_dir": None,
        "logs_dir": str(run_dir(run_id) / "logs"),
        "params": {},
        "nextflow_profile": "local,docker",
    }


def copy_inputs(input_dir: Path, target: Path) -> None:
    if not input_dir.exists() or not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    target.mkdir(parents=True, exist_ok=True)
    for item in input_dir.iterdir():
        dest = target / item.name
        if item.is_dir():
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(item, dest)
        elif item.is_file():
            shutil.copy2(item, dest)


def fastq_files(input_dir: Path) -> list[Path]:
    patterns = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    files: list[Path] = []
    for pattern in patterns:
        files.extend(input_dir.rglob(pattern))
    return sorted(path for path in files if path.is_file())


def sample_id_from_r1(name: str) -> str:
    base = re.sub(r"\.(fastq|fq)(\.gz)?$", "", name)
    base = re.sub(r"([._-])R?1$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"([._-])read1$", "", base, flags=re.IGNORECASE)
    return base or "sample"


def pair_fastqs(files: list[Path]) -> list[dict[str, str]]:
    by_name = {path.name: path for path in files}
    rows = []
    used = set()
    for path in files:
        if path.name in used:
            continue
        mate_name = None
        for old, new in (("_R1", "_R2"), ("_1", "_2"), (".R1", ".R2"), ("-R1", "-R2"), ("read1", "read2")):
            if old in path.name:
                candidate = path.name.replace(old, new, 1)
                if candidate in by_name:
                    mate_name = candidate
                    break
        sample_id = sample_id_from_r1(path.name)
        if mate_name:
            used.add(path.name)
            used.add(mate_name)
            rows.append({
                "sample_id": sample_id,
                "fastq_1": str(path.resolve()),
                "fastq_2": str(by_name[mate_name].resolve()),
                "single_end": "false",
                "strandedness": "unknown",
            })
        else:
            used.add(path.name)
            rows.append({
                "sample_id": sample_id,
                "fastq_1": str(path.resolve()),
                "fastq_2": "",
                "single_end": "true",
                "strandedness": "unknown",
            })
    return rows


def find_manifest(input_dir: Path, kind: str) -> Path | None:
    preferred = (
        ["trim_manifest.tsv", "trimmed_fastq_manifest.tsv"] if kind == "trimmed"
        else ["raw_manifest.csv", "samplesheet.csv", "samplesheet.tsv"]
    )
    for name in preferred:
        candidate = input_dir / name
        if candidate.exists() and candidate.is_file():
            return candidate
    for candidate in sorted(input_dir.glob("*.tsv")) + sorted(input_dir.glob("*.csv")):
        text = candidate.read_text(errors="replace", encoding="utf-8")
        first = text.splitlines()[0] if text.splitlines() else ""
        if "fastq_1" in first:
            return candidate
    return None


def write_manifest_from_fastqs(input_dir: Path, kind: str) -> Path:
    rows = pair_fastqs(fastq_files(input_dir))
    if not rows:
        raise FileNotFoundError(
            f"No manifest and no FASTQ files were found in {input_dir}. "
            "Provide raw_manifest.csv, trim_manifest.tsv, or FASTQ files."
        )
    suffix = "tsv" if kind == "trimmed" else "csv"
    manifest = input_dir / (f"{kind}_manifest.{suffix}" if kind != "raw" else "raw_manifest.csv")
    delimiter = "\t" if manifest.suffix == ".tsv" else ","
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "fastq_1", "fastq_2", "single_end", "strandedness"], delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)
    return manifest


def input_manifest_for_case(case: Case, input_dir: Path) -> Path:
    if case.input_kind == "nofastq":
        manifest = input_dir / "no_fastq.csv"
        manifest.write_text("sample_id,fastq_1,fastq_2,single_end,strandedness,platform\n", encoding="utf-8")
        return manifest
    manifest = find_manifest(input_dir, case.input_kind)
    if manifest:
        return manifest
    return write_manifest_from_fastqs(input_dir, case.input_kind)


def compiler() -> NextflowCompiler:
    return NextflowCompiler(load_steps(PIPELINE_ROOT / "registry" / "steps.yaml"), PIPELINE_ROOT / "workflows" / "templates")


def compile_workflow(case: Case, manifest: Path, state: dict) -> Path:
    run_id = state["run_id"]
    compiled_dir = PIPELINE_ROOT / "workflows" / "generated" / f"local_{run_id}"
    output_dir = Path(state["output_dir"])
    params = reference_params(case.route, PIPELINE_ROOT / "assets" / "demo_reference")
    if case.input_kind == "trimmed":
        params.update({"trim_manifest": str(manifest), "trim_input_mode": "previous_manifest"})
    if case.extra_params:
        params.update(case.extra_params)
    params.update(state.get("params") or {})
    request = RunRequest(
        run_id=run_id,
        omics="rnaseq",
        execution_profile=state.get("nextflow_profile") or "local,docker",
        input=str(manifest),
        outdir=str(output_dir),
        selected_steps=case.selected_steps,
        params=params,
    )
    compiled = compiler().compile(request, compiled_dir)
    plan = json.loads((compiled / "workflow_plan.json").read_text(encoding="utf-8"))
    state["compiled_dir"] = str(compiled)
    state["steps"] = [
        {
            "step_id": step["step_id"],
            "process_name": step["process_name"],
            "status": "pending",
            "outputs": [],
        }
        for step in plan.get("steps", [])
    ]
    state["status"] = "compiled"
    state["current_step"] = state["steps"][0]["step_id"] if state["steps"] else None
    save_state(state)
    return compiled


def update_state_from_execution(state: dict, execution_record: Path, status: str) -> None:
    record = json.loads(execution_record.read_text(encoding="utf-8"))
    completed = []
    failed = None
    current = None
    steps = []
    for step in record.get("step_outputs", []):
        step_status = step.get("status", "pending")
        if step_status == "completed":
            completed.append(step.get("step_id"))
        elif step_status == "failed" and failed is None:
            failed = step.get("step_id")
        elif current is None:
            current = step.get("step_id")
        steps.append({
            "step_id": step.get("step_id"),
            "process_name": step.get("process_name"),
            "status": step_status,
            "outputs": step.get("outputs", []),
            "missing_outputs": step.get("missing_outputs", []),
        })
    state["steps"] = steps
    state["completed_steps"] = [step for step in completed if step]
    state["failed_step"] = failed
    state["current_step"] = failed or current
    state["status"] = "completed" if status == "succeeded" and not failed else status
    save_state(state)


def outputs_valid(state: dict) -> bool:
    if state.get("status") != "completed":
        return False
    output_dir = Path(state.get("output_dir") or "")
    return output_dir.exists() and any(output_dir.rglob("*"))


def run_nextflow(compiled: Path, state: dict) -> int:
    logs_dir = Path(state["logs_dir"])
    logs_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = logs_dir / "nextflow.stdout.log"
    stderr_path = logs_dir / "nextflow.stderr.log"
    profile = state.get("nextflow_profile") or "local,docker"
    command = ["nextflow", "run", "main.nf", "-params-file", "params.yaml", "-profile", profile, "-work-dir", "work", "-resume"]
    state["status"] = "running"
    state["failed_step"] = None
    save_state(state)
    with stdout_path.open("a", encoding="utf-8") as stdout, stderr_path.open("a", encoding="utf-8") as stderr:
        proc = subprocess.Popen(command, cwd=compiled, stdout=stdout, stderr=stderr, text=True)
        while proc.poll() is None:
            time.sleep(10)
            try:
                execution_record = collect_run_outputs(compiled, status="running")
                update_state_from_execution(state, execution_record, "running")
            except Exception:
                pass
    status = "succeeded" if proc.returncode == 0 else "failed"
    execution_record = collect_run_outputs(compiled, status=status)
    shutil.copy2(compiled / ".nextflow.log", logs_dir / ".nextflow.log") if (compiled / ".nextflow.log").exists() else None
    update_state_from_execution(state, execution_record, status)
    return proc.returncode


def start_run(args: argparse.Namespace) -> int:
    case = CASE_BY_NAME.get(args.workflow)
    if case is None:
        allowed = ", ".join(sorted(CASE_BY_NAME))
        raise ValueError(f"Unknown workflow '{args.workflow}'. Available workflows: {allowed}")
    input_dir = Path(args.input).expanduser().resolve()
    if not input_dir.exists() or not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    run_id = args.run_id or new_run_id(case.name)
    folder = run_dir(run_id)
    inputs_dir = folder / "inputs"
    output_dir = Path(args.output).expanduser().resolve() if args.output else folder / "outputs"
    state = initial_state(run_id, case.name, output_dir)
    state["params"] = parse_param_overrides(args.param)
    state["nextflow_profile"] = args.profile
    save_state(state)
    copy_inputs(input_dir, inputs_dir)
    manifest = input_manifest_for_case(case, inputs_dir)
    state["input_manifest"] = str(manifest)
    save_state(state)
    compiled = compile_workflow(case, manifest, state)
    code = run_nextflow(compiled, state)
    print(f"Run ID: {run_id}")
    print(f"State: {state_path(run_id)}")
    print(f"Outputs: {output_dir}")
    return code


def resume_run(run_id: str) -> int:
    state = load_state(run_id)
    if outputs_valid(state):
        print(f"Run '{run_id}' is already completed and outputs are present.")
        print(f"State: {state_path(run_id)}")
        print(f"Outputs: {state['output_dir']}")
        return 0
    compiled = Path(state.get("compiled_dir") or "")
    if not compiled.exists():
        raise FileNotFoundError(f"Cannot resume '{run_id}': compiled workflow is missing: {compiled}")
    print(f"Resuming run '{run_id}' from {state.get('failed_step') or state.get('current_step') or 'first incomplete step'}")
    return run_nextflow(compiled, state)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a SurvOm workflow locally with saved state and resume support.")
    parser.add_argument("--workflow", help="Workflow name, for example trim_only or qc_trim_strandedness.")
    parser.add_argument("--input", help="Directory containing FASTQ files or an input manifest.")
    parser.add_argument("--output", help="Output directory. Defaults to runs/<run_id>/outputs.")
    parser.add_argument("--resume", help="Resume an existing run ID from runs/<run_id>.")
    parser.add_argument("--run-id", help="Optional fixed run ID for a new run.")
    parser.add_argument("--profile", default="local,docker", help="Nextflow profile for a new run. Default: local,docker.")
    parser.add_argument(
        "--param",
        action="append",
        help="Override a workflow parameter. Repeatable. Example: --param quality_threshold=25",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.resume:
            return resume_run(args.resume)
        missing = [name for name in ("workflow", "input") if not getattr(args, name)]
        if missing:
            raise ValueError(f"Missing required option(s): {', '.join('--' + name for name in missing)}")
        return start_run(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
