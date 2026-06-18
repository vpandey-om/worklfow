from __future__ import annotations

import csv
import json
import os
import subprocess
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from .job_store import JobStore
from .upload_service import UploadService, safe_filename, safe_session_id, safe_slug


class WorkflowService:
    def __init__(self, app_root: Path, pipeline_root: Path, upload_service: UploadService, job_store: JobStore):
        self.app_root = Path(app_root)
        self.pipeline_root = Path(pipeline_root)
        self.upload_service = upload_service
        self.job_store = job_store
        self.runs_root = self.app_root / "runs"
        self.runs_root.mkdir(parents=True, exist_ok=True)

    def load_steps(self, config_path: Path) -> list[dict[str, Any]]:
        data = yaml.safe_load(Path(config_path).read_text())
        return data.get("steps", [])

    def load_yaml_config(self, config_path: Path) -> dict[str, Any]:
        return yaml.safe_load(Path(config_path).read_text()) or {}

    def create_or_find_samplesheet(
        self,
        session_id: str,
        run_dir: Path,
        tester_id: str,
        omics_type: str,
        project_id: str,
        selected_files: dict[str, list[str]] | None = None,
    ) -> Path:
        selected_files = selected_files or {}
        metadata_files = selected_files.get("metadata") or []
        if not metadata_files:
            uploads = self.upload_service.list_uploads(session_id, tester_id, omics_type, project_id)
            metadata_files = [item["path"] for item in uploads["metadata"]]
        if metadata_files:
            source = Path(metadata_files[0])
            target = run_dir / safe_filename(source.name)
            target.write_bytes(source.read_bytes())
            if target.suffix.lower() == ".csv":
                return target
            raise ValueError("For this demo runner, please upload a CSV sample sheet.")

        fastq_files = selected_files.get("fastq") or []
        if not fastq_files:
            uploads = self.upload_service.list_uploads(session_id, tester_id, omics_type, project_id)
            fastq_files = [item["path"] for item in uploads["fastq"]]
        fastqs = [Path(path) for path in fastq_files]
        if not fastqs:
            raise ValueError("Upload FASTQ files or a sample sheet before running.")
        r1 = next((p for p in fastqs if "_R1" in p.name or "_1" in p.name), fastqs[0])
        r2 = next((p for p in fastqs if p != r1 and ("_R2" in p.name or "_2" in p.name)), None)
        samplesheet = run_dir / "samplesheet.csv"
        with samplesheet.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample_id", "fastq_1", "fastq_2", "library_layout", "platform"])
            writer.writeheader()
            writer.writerow({
                "sample_id": "uploaded_sample",
                "fastq_1": str(r1),
                "fastq_2": str(r2) if r2 else "",
                "library_layout": "paired" if r2 else "single",
                "platform": "unknown",
            })
        return samplesheet

    def selected_steps_for_workflow(self, workflow_id: str) -> list[str]:
        if workflow_id == "fastq_qc":
            return ["rnaseq_03_raw_read_qc"]
        if workflow_id == "trim_only":
            return ["rnaseq_04_adapter_quality_trimming"]
        if workflow_id == "qc_trim":
            return ["rnaseq_03_raw_read_qc", "rnaseq_04_adapter_quality_trimming"]
        raise ValueError(f"Unknown workflow id: {workflow_id}")

    @staticmethod
    def _bool_param(value: Any, default: bool = False) -> bool:
        if value is None:
            return default
        if isinstance(value, bool):
            return value
        return str(value).lower() in {"1", "true", "yes", "on"}

    @staticmethod
    def _optional_param(value: Any) -> Any:
        return None if value == "" else value

    def create_run(
        self,
        session_id: str,
        workflow_id: str,
        params: dict[str, Any],
        tester_id: str | None = None,
        tester_label: str | None = None,
        omics_type: str = "bulk_rnaseq",
        project_id: str = "demo_project",
        project_label: str | None = None,
        selected_files: dict[str, list[str]] | None = None,
    ) -> dict[str, Any]:
        sid = safe_session_id(session_id)
        tester = safe_slug(tester_id or "unknown_tester")
        omics = safe_slug(omics_type or "bulk_rnaseq")
        project = safe_slug(project_id or "demo_project")
        run_id = f"{workflow_id}_{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
        run_dir = self.runs_root / omics / tester / project / run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        selected_files = selected_files or {}
        samplesheet = self.create_or_find_samplesheet(sid, run_dir, tester, omics, project, selected_files)

        selected_steps = self.selected_steps_for_workflow(workflow_id)
        outdir = run_dir / "results"
        resolved_params = {
            "qc_tool": params.get("qc_tool", "fastqc"),
            "trimming_tool": params.get("trimming_tool", "fastp"),
            "quality_threshold": int(params.get("quality_threshold", 20)),
            "minimum_read_length": int(params.get("minimum_read_length", 20)),
            "trim_poly_g": str(params.get("trim_poly_g", "auto")).lower(),
            "trim_poly_x": self._bool_param(params.get("trim_poly_x"), False),
            "adapter_preset": self._optional_param(params.get("adapter_preset")),
            "adapter_source": self._optional_param(params.get("adapter_source")),
            "adapter_sequence_r1": self._optional_param(params.get("adapter_sequence_r1")),
            "adapter_sequence_r2": self._optional_param(params.get("adapter_sequence_r2")),
            "adapter_fasta": self._optional_param(params.get("adapter_fasta")),
            "trim_front_r1": int(params.get("trim_front_r1", 0) or 0),
            "trim_front_r2": int(params.get("trim_front_r2", 0) or 0),
            "detect_adapter_for_pe": self._bool_param(params.get("detect_adapter_for_pe"), True),
            "cutadapt_error_rate": float(params.get("cutadapt_error_rate", 0.1) or 0.1),
            "cutadapt_minimum_overlap": int(params.get("cutadapt_minimum_overlap", 3) or 3),
            "threads": int(params.get("threads", 4) or 4),
        }
        if "rnaseq_04_adapter_quality_trimming" not in selected_steps:
            resolved_params = {
                "qc_tool": resolved_params["qc_tool"],
                "threads": resolved_params["threads"],
            }
        run_request = {
            "run_id": run_id,
            "omics": "rnaseq",
            "omics_type": omics,
            "tester_id": tester,
            "tester_label": tester_label or tester,
            "project_id": project,
            "project_label": project_label or project,
            "workflow_id": workflow_id,
            "selected_files": selected_files,
            "execution_profile": params.get("execution_profile", "local"),
            "input": str(samplesheet),
            "outdir": str(outdir),
            "selected_steps": selected_steps,
            "params": resolved_params,
            "dataset_id": f"{omics}:{tester}:{project}:{sid}",
            "user_id": tester,
            "run_label": params.get("run_label", run_id),
        }
        request_path = run_dir / "run_request.yaml"
        request_path.write_text(yaml.safe_dump(run_request, sort_keys=False), encoding="utf-8")

        compiled_dir = self.pipeline_root / "workflows" / "generated" / run_id
        compile_cmd = [
            "python3", "-m", "builder.api", "compile",
            "--request", str(request_path),
            "--output", str(compiled_dir),
        ]
        env = os.environ.copy()
        env["JAVA_CMD"] = "/usr/bin/java"
        env.pop("JAVA_HOME", None)
        env.pop("JAVA_LD_LIBRARY_PATH", None)
        compile_proc = subprocess.run(
            compile_cmd,
            cwd=self.pipeline_root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            check=False,
        )
        (run_dir / "compile.stdout.log").write_text(compile_proc.stdout, encoding="utf-8")
        (run_dir / "compile.stderr.log").write_text(compile_proc.stderr, encoding="utf-8")
        if compile_proc.returncode != 0:
            raise RuntimeError(f"Workflow compile failed: {compile_proc.stderr}")

        job_id = f"job_{run_id}"
        job_json = {
            "job_id": job_id,
            "run_id": run_id,
            "session_id": sid,
            "tester_id": tester,
            "tester_label": tester_label or tester,
            "omics_type": omics,
            "project_id": project,
            "project_label": project_label or project,
            "workflow_id": workflow_id,
            "selected_files": selected_files,
            "compiled_dir": str(compiled_dir),
            "command": f"bash {compiled_dir / 'run_command.sh'}",
            "results_dir": str(outdir),
            "execution_profile": run_request["execution_profile"],
            "logs_dir": str(compiled_dir / "logs"),
            "params": resolved_params,
            "submitted_params": resolved_params,
        }
        job_path = run_dir / "job.json"
        job_path.write_text(json.dumps(job_json, indent=2), encoding="utf-8")
        self.job_store.upsert_job({
            "job_id": job_id,
            "run_id": run_id,
            "session_id": sid,
            "tester_id": tester,
            "tester_label": tester_label or tester,
            "omics_type": omics,
            "project_id": project,
            "project_label": project_label or project,
            "workflow_id": workflow_id,
            "status": "compiled",
            "pid": None,
            "run_dir": str(run_dir),
            "compiled_dir": str(compiled_dir),
            "results_dir": str(outdir),
            "command": job_json["command"],
            "stdout_path": str(compiled_dir / "logs" / "stdout.log"),
            "stderr_path": str(compiled_dir / "logs" / "stderr.log"),
            "metadata_json": json.dumps({
                "request_path": str(request_path),
                "job_path": str(job_path),
                "submitted_params": resolved_params,
                "tester_id": tester,
                "tester_label": tester_label or tester,
                "omics_type": omics,
                "project_id": project,
                "project_label": project_label or project,
                "workflow_id": workflow_id,
                "selected_files": selected_files,
            }),
        })
        return {"run_id": run_id, "job_id": job_id, "run_dir": str(run_dir), "job_path": str(job_path)}

    def start_job(self, job_id: str) -> dict[str, Any]:
        job = self.job_store.get_job(job_id)
        if not job:
            raise ValueError(f"Unknown job id: {job_id}")
        metadata = json.loads(job.get("metadata_json") or "{}")
        job_path = metadata.get("job_path")
        if not job_path:
            raise ValueError("Job is missing job_path metadata")
        cmd = ["python3", "-m", "builder.api", "run-job", "--job", job_path]
        env = os.environ.copy()
        env["JAVA_CMD"] = "/usr/bin/java"
        env.pop("JAVA_HOME", None)
        env.pop("JAVA_LD_LIBRARY_PATH", None)
        proc = subprocess.Popen(
            cmd,
            cwd=self.pipeline_root,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=env,
            start_new_session=True,
        )
        self.job_store.update_job(
            job_id,
            status="running",
            pid=proc.pid,
            started_at=datetime.now(timezone.utc).isoformat(),
        )
        return {"job_id": job_id, "pid": proc.pid, "status": "running"}

    def refresh_job(self, job_id: str) -> dict[str, Any] | None:
        job = self.job_store.get_job(job_id)
        if not job:
            return None
        if job.get("compiled_dir"):
            self.job_store.refresh_from_job_record(job_id, Path(job["compiled_dir"]))
            job = self.job_store.get_job(job_id)
        return job

    def output_files(self, job_id: str) -> list[str]:
        job = self.refresh_job(job_id)
        if not job:
            return []
        record = Path(job["compiled_dir"]) / "execution_record.json"
        if record.exists():
            data = json.loads(record.read_text())
            return data.get("all_result_files", [])
        result_dir = Path(job["results_dir"])
        return sorted(str(path) for path in result_dir.rglob("*") if path.is_file()) if result_dir.exists() else []
