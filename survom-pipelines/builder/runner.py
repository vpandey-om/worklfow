import subprocess
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

from .models import JobRequest
from .outputs import collect_run_outputs


MISSING_CUTADAPT_MESSAGE = (
    "Cutadapt was selected but is not available in the active environment. "
    "Install cutadapt or run with a container/conda profile."
)
MISSING_SALMON_MESSAGE = (
    "Salmon strandedness inference was selected but salmon is not available in the active environment. "
    "Install salmon or run with a container/conda profile."
)
MISSING_RSEQC_MESSAGE = (
    "RSeQC strandedness validation was selected but STAR/RSeQC tools are not available in the active environment. "
    "Install STAR and RSeQC or run with a container/conda profile."
)


def _uses_container_profile(profile: str) -> bool:
    profile_tokens = {token.strip() for token in (profile or "").split(",")}
    return bool(profile_tokens & {"docker", "singularity", "aws", "local_docker"})


def _validate_runtime_tools(manifest: dict, profile: str):
    params = manifest.get("params", {})
    if _uses_container_profile(profile):
        return
    if params.get("trimming_tool") == "cutadapt" and not shutil.which("cutadapt"):
        raise RuntimeError(MISSING_CUTADAPT_MESSAGE)
    if (params.get("strandedness_method") == "salmon_auto" or params.get("selected_route") == "salmon") and not shutil.which("salmon"):
        raise RuntimeError(MISSING_SALMON_MESSAGE)
    if params.get("strandedness_method") == "rseqc_validation":
        if not shutil.which("infer_experiment.py") or (not params.get("existing_bam") and not shutil.which("STAR")):
            raise RuntimeError(MISSING_RSEQC_MESSAGE)


def run_nextflow(compiled_dir: str | Path, profile: str = "local,docker", resume: bool = True) -> subprocess.CompletedProcess:
    compiled = Path(compiled_dir)
    manifest_path = compiled / "workflow_manifest.json"
    if manifest_path.exists():
        _validate_runtime_tools(json.loads(manifest_path.read_text()), profile)
    cmd = [
        "nextflow",
        "run",
        str(compiled / "main.nf"),
        "-profile",
        profile,
        "-params-file",
        str(compiled / "params.yaml"),
        "-work-dir",
        str(compiled / "work"),
    ]
    if resume:
        cmd.append("-resume")
    return subprocess.run(cmd, cwd=compiled, text=True, capture_output=True, check=False)


def run_job(job_json: str | Path) -> Path:
    job_path = Path(job_json)
    job = JobRequest.model_validate(json.loads(job_path.read_text()))
    compiled_dir = Path(job.compiled_dir)
    manifest_path = compiled_dir / "workflow_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing compiled workflow manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("run_id") != job.run_id:
        raise ValueError(
            f"Job run_id '{job.run_id}' does not match compiled workflow run_id "
            f"'{manifest.get('run_id')}'. Compile a new workflow for the new run_id "
            "or point the job to the matching compiled_dir."
        )
    compiled_outdir = manifest.get("params", {}).get("outdir")
    if compiled_outdir and compiled_outdir != job.results_dir:
        raise ValueError(
            f"Job results_dir '{job.results_dir}' does not match compiled workflow outdir "
            f"'{compiled_outdir}'. Recompile or update the job JSON to match."
        )
    logs_dir = Path(job.logs_dir) if job.logs_dir else compiled_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    record = {
        "job_id": job.job_id,
        "run_id": job.run_id,
        "compiled_dir": job.compiled_dir,
        "results_dir": job.results_dir,
        "execution_profile": job.execution_profile,
        "command": job.command,
        "status": "running",
        "started_at": datetime.now(timezone.utc).isoformat(),
        "finished_at": None,
        "exit_code": None,
        "stdout_path": str(logs_dir / "stdout.log"),
        "stderr_path": str(logs_dir / "stderr.log"),
    }
    job_record_path = compiled_dir / "job_record.json"
    job_record_path.write_text(json.dumps(record, indent=2), encoding="utf-8")

    env = os.environ.copy()
    env["JAVA_CMD"] = "/usr/bin/java"
    env.pop("JAVA_HOME", None)
    env.pop("JAVA_LD_LIBRARY_PATH", None)

    try:
        _validate_runtime_tools(manifest, job.execution_profile)
        proc = subprocess.run(
            job.command,
            shell=True,
            cwd=Path.cwd(),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            check=False,
        )
        stdout = proc.stdout
        stderr = proc.stderr
        returncode = proc.returncode
    except RuntimeError as exc:
        stdout = ""
        stderr = str(exc)
        returncode = 127

    (logs_dir / "stdout.log").write_text(stdout, encoding="utf-8")
    (logs_dir / "stderr.log").write_text(stderr, encoding="utf-8")

    status = "succeeded" if returncode == 0 else "failed"
    record.update({
        "status": status,
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "exit_code": returncode,
    })
    job_record_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    collect_run_outputs(compiled_dir, status=status)
    return job_record_path
