import json
from datetime import datetime, timezone
from pathlib import Path


def collect_run_outputs(compiled_dir: str | Path, status: str = "unknown") -> Path:
    compiled = Path(compiled_dir)
    manifest_path = compiled / "workflow_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing workflow manifest: {manifest_path}")

    manifest = json.loads(manifest_path.read_text())
    outdir = Path(manifest["params"]["outdir"])
    step_outputs = []
    for step in manifest.get("step_records", []):
        found = []
        for output in step.get("expected_outputs", []):
            pattern = output.get("pattern")
            if not pattern:
                continue
            matches = sorted(str(path) for path in outdir.glob(pattern) if path.is_file())
            found.append({
                "name": output.get("name"),
                "type": output.get("type"),
                "pattern": pattern,
                "files": matches,
                "count": len(matches),
                "ui_visible": output.get("ui_visible", True),
            })
        step_outputs.append({
            "step_id": step["step_id"],
            "step_name": step["step_name"],
            "process_name": step["process_name"],
            "outputs": found,
        })

    execution_record = {
        "run_id": manifest["run_id"],
        "status": status,
        "collected_at": datetime.now(timezone.utc).isoformat(),
        "compiled_dir": str(compiled),
        "results_dir": str(outdir),
        "nextflow_command": manifest.get("nextflow_command"),
        "nextflow_log_path": ".nextflow.log",
        "step_outputs": step_outputs,
        "all_result_files": sorted(str(path) for path in outdir.rglob("*") if path.is_file()) if outdir.exists() else [],
    }
    output_path = compiled / "execution_record.json"
    output_path.write_text(json.dumps(execution_record, indent=2), encoding="utf-8")
    return output_path
