import json
import re
from datetime import datetime, timezone
from pathlib import Path


def _read_text_tail(path: Path, limit: int = 4000) -> str:
    if not path.exists() or not path.is_file():
        return ""
    return path.read_text(errors="replace")[-limit:]


def _parse_task_name(command_run: Path) -> tuple[str | None, str | None]:
    for line in _read_text_tail(command_run, 2000).splitlines():
        match = re.match(r"### name: '([^']+)'", line)
        if match:
            label = match.group(1)
            process = label.split(" (", 1)[0].strip()
            return process, label
    return None, None


def _task_artifacts(compiled: Path) -> list[dict]:
    tasks = []
    work_dir = compiled / "work"
    if not work_dir.exists():
        return tasks
    for command_run in sorted(work_dir.glob("*/*/.command.run")):
        task_dir = command_run.parent
        process, label = _parse_task_name(command_run)
        exit_path = task_dir / ".exitcode"
        exit_code = None
        if exit_path.exists():
            try:
                exit_code = int(exit_path.read_text().strip())
            except ValueError:
                exit_code = None
        if exit_code is None:
            status = "running" if (task_dir / ".command.begin").exists() else "pending"
        else:
            status = "completed" if exit_code == 0 else "failed"
        tasks.append({
            "process_name": process,
            "task_name": label,
            "status": status,
            "exit_code": exit_code,
            "work_dir": str(task_dir),
            "extra_log_paths": sorted(
                str(path)
                for pattern in ("*.log", "*.out", "*.err", "*.txt")
                for path in task_dir.glob(pattern)
                if path.is_file() and not path.name.startswith(".command")
            ),
            "stdout_path": str(task_dir / ".command.out"),
            "stderr_path": str(task_dir / ".command.err"),
            "command_log_path": str(task_dir / ".command.log"),
            "command_sh_path": str(task_dir / ".command.sh"),
            "stdout_tail": _read_text_tail(task_dir / ".command.out"),
            "stderr_tail": _read_text_tail(task_dir / ".command.err"),
            "log_tail": _read_text_tail(task_dir / ".command.log"),
        })
    return tasks


def collect_run_outputs(compiled_dir: str | Path, status: str = "unknown") -> Path:
    compiled = Path(compiled_dir)
    manifest_path = compiled / "workflow_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing workflow manifest: {manifest_path}")

    manifest = json.loads(manifest_path.read_text())
    outdir = Path(manifest["params"]["outdir"])
    tasks = _task_artifacts(compiled)
    tasks_by_process: dict[str, list[dict]] = {}
    for task in tasks:
        if task.get("process_name"):
            tasks_by_process.setdefault(task["process_name"], []).append(task)

    plan = manifest.get("workflow_plan") or {}
    planned_steps = plan.get("steps", [])
    step_outputs = []
    plan_by_process = {step.get("process_name"): step for step in planned_steps}
    manifest_steps = list(manifest.get("step_records", []))
    for plan_step in planned_steps:
        if not any(step.get("process_name") == plan_step.get("process_name") for step in manifest_steps):
            manifest_steps.append({
                "step_id": plan_step.get("step_id"),
                "step_name": plan_step.get("step_name"),
                "process_name": plan_step.get("process_name"),
                "expected_outputs": plan_step.get("outputs", []),
            })
    for step in manifest_steps:
        found = []
        raw_missing = []
        for output in step.get("expected_outputs", []):
            pattern = output.get("pattern")
            if not pattern:
                continue
            matches = sorted(str(path) for path in outdir.glob(pattern) if path.exists())
            if not matches:
                raw_missing.append({
                    "name": output.get("name"),
                    "pattern": pattern,
                })
            found.append({
                "name": output.get("name"),
                "type": output.get("type"),
                "pattern": pattern,
                "files": matches,
                "count": len(matches),
                "ui_visible": output.get("ui_visible", True),
            })
        process_tasks = tasks_by_process.get(step["process_name"], [])
        terminal_workflow = status in {"succeeded", "completed", "success", "failed", "error"}
        if any(task.get("status") == "failed" for task in process_tasks):
            step_status = "failed"
            missing = raw_missing
        elif process_tasks and all(task.get("status") == "completed" for task in process_tasks):
            step_status = "failed" if raw_missing else "completed"
            missing = raw_missing
        elif process_tasks:
            step_status = "running"
            missing = []
        elif terminal_workflow and raw_missing:
            step_status = "failed"
            missing = raw_missing
        else:
            step_status = "pending"
            missing = []
        step_outputs.append({
            "step_id": step["step_id"],
            "step_name": step["step_name"],
            "process_name": step["process_name"],
            "status": step_status,
            "outputs": found,
            "missing_outputs": missing,
            "tasks": process_tasks,
            "dependencies": (plan_by_process.get(step["process_name"]) or {}).get("dependencies", []),
        })

    execution_record = {
        "run_id": manifest["run_id"],
        "status": status,
        "collected_at": datetime.now(timezone.utc).isoformat(),
        "compiled_dir": str(compiled),
        "results_dir": str(outdir),
        "nextflow_command": manifest.get("nextflow_command"),
        "nextflow_log_path": str(compiled / ".nextflow.log"),
        "workflow_plan_path": str(compiled / "workflow_plan.json"),
        "task_artifacts": tasks,
        "step_outputs": step_outputs,
        "all_result_files": sorted(str(path) for path in outdir.rglob("*") if path.is_file()) if outdir.exists() else [],
        "all_result_dirs": sorted(str(path) for path in outdir.rglob("*") if path.is_dir()) if outdir.exists() else [],
    }
    output_path = compiled / "execution_record.json"
    output_path.write_text(json.dumps(execution_record, indent=2), encoding="utf-8")
    return output_path
