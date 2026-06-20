from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
from typing import Any


class DownstreamService:
    def __init__(self, repo_root: Path):
        self.repo_root = Path(repo_root)
        self.pipeline_root = self.repo_root / "survom-pipelines"
        self.downstream_root = self.pipeline_root / "downstream"

    def create_downstream_run(
        self,
        run_id: str,
        featurecounts_input: str,
        metadata: str,
        runs_dir: str,
        contrasts: str | None = None,
        gene_mapping: str | None = None,
        group_column: str = "condition",
    ) -> dict[str, Any]:
        run_dir = Path(runs_dir) / run_id
        lock = run_dir / ".lock"
        run_dir.mkdir(parents=True, exist_ok=True)
        if lock.exists():
            raise RuntimeError(f"Downstream run is locked: {run_id}")
        request = {
            "run_id": run_id,
            "featurecounts_input": featurecounts_input,
            "metadata": metadata,
            "contrasts": contrasts,
            "gene_mapping": gene_mapping,
            "runs_dir": str(Path(runs_dir)),
            "group_column": group_column,
        }
        (run_dir / "downstream_request.json").write_text(json.dumps(request, indent=2) + "\n", encoding="utf-8")
        return {"run_id": run_id, "run_dir": str(run_dir), "request": request}

    def validate_downstream_inputs(self, featurecounts_input: str, metadata: str) -> dict[str, Any]:
        errors = []
        for label, value in (("featurecounts_input", featurecounts_input), ("metadata", metadata)):
            path = Path(value)
            if not path.exists():
                errors.append(f"{label} does not exist: {value}")
        return {"valid": not errors, "errors": errors}

    def launch_downstream_nextflow_run(self, run_id: str, runs_dir: str, profile: str = "local") -> dict[str, Any]:
        run_dir = Path(runs_dir) / run_id
        request_path = run_dir / "downstream_request.json"
        if not request_path.exists():
            raise FileNotFoundError(f"Missing downstream request: {request_path}")
        request = json.loads(request_path.read_text(encoding="utf-8"))
        lock = run_dir / ".lock"
        lock.write_text(str(os.getpid()), encoding="utf-8")
        command = [
            "nextflow",
            "-C",
            "nextflow_downstream.config",
            "run",
            "main_downstream.nf",
            "-profile",
            profile,
            "-resume",
            "--featurecounts_input",
            request["featurecounts_input"],
            "--metadata",
            request["metadata"],
            "--run_id",
            run_id,
            "--runs_dir",
            str(Path(runs_dir).resolve()),
            "--group_column",
            request.get("group_column") or "condition",
        ]
        if request.get("contrasts"):
            command.extend(["--contrasts", request["contrasts"]])
        if request.get("gene_mapping"):
            command.extend(["--gene_mapping", request["gene_mapping"]])
        logs = run_dir / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        with (logs / "stdout.log").open("ab") as stdout, (logs / "stderr.log").open("ab") as stderr:
            proc = subprocess.Popen(command, cwd=self.downstream_root, stdout=stdout, stderr=stderr)
        (run_dir / "nextflow_command.txt").write_text(" ".join(command) + "\n", encoding="utf-8")
        (run_dir / "pid.txt").write_text(str(proc.pid), encoding="utf-8")
        return {"run_id": run_id, "pid": proc.pid, "command": command}

    def get_downstream_run_status(self, run_id: str, runs_dir: str) -> dict[str, Any]:
        state = Path(runs_dir) / run_id / "state" / "run_state.json"
        if state.exists():
            return json.loads(state.read_text(encoding="utf-8"))
        return {"run_id": run_id, "status": "unknown", "state_path": str(state)}

    def list_downstream_outputs(self, run_id: str, runs_dir: str) -> list[str]:
        results = Path(runs_dir) / run_id / "results"
        if not results.exists():
            return []
        return sorted(str(path) for path in results.rglob("*") if path.is_file())

    def read_downstream_csv(self, run_id: str, runs_dir: str, step_name: str, filename: str) -> str:
        path = Path(runs_dir) / run_id / "results" / step_name / filename
        if not path.exists() or not path.is_file():
            raise FileNotFoundError(f"Missing downstream CSV: {path}")
        return path.read_text(encoding="utf-8")

    def resume_downstream_run(self, run_id: str, runs_dir: str, profile: str = "local") -> dict[str, Any]:
        return self.launch_downstream_nextflow_run(run_id, runs_dir, profile=profile)


def create_downstream_run(*args, **kwargs):
    return DownstreamService(Path(__file__).resolve().parents[2]).create_downstream_run(*args, **kwargs)


def validate_downstream_inputs(*args, **kwargs):
    return DownstreamService(Path(__file__).resolve().parents[2]).validate_downstream_inputs(*args, **kwargs)


def launch_downstream_nextflow_run(*args, **kwargs):
    return DownstreamService(Path(__file__).resolve().parents[2]).launch_downstream_nextflow_run(*args, **kwargs)


def get_downstream_run_status(*args, **kwargs):
    return DownstreamService(Path(__file__).resolve().parents[2]).get_downstream_run_status(*args, **kwargs)


def list_downstream_outputs(*args, **kwargs):
    return DownstreamService(Path(__file__).resolve().parents[2]).list_downstream_outputs(*args, **kwargs)


def read_downstream_csv(*args, **kwargs):
    return DownstreamService(Path(__file__).resolve().parents[2]).read_downstream_csv(*args, **kwargs)


def resume_downstream_run(*args, **kwargs):
    return DownstreamService(Path(__file__).resolve().parents[2]).resume_downstream_run(*args, **kwargs)
