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
            if target.suffix.lower() in {".csv", ".tsv"}:
                return target
            raise ValueError("For this demo runner, please upload a CSV or TSV sample sheet/manifest.")

        fastq_files = selected_files.get("fastq") or []
        if not fastq_files:
            uploads = self.upload_service.list_uploads(session_id, tester_id, omics_type, project_id)
            fastq_files = [item["path"] for item in uploads["fastq"]]
        fastqs = [Path(path) for path in fastq_files]
        if not fastqs:
            raise ValueError("Upload FASTQ files or a sample sheet before running.")
        rows = self._fastqs_to_samplesheet_rows(fastqs)
        samplesheet = run_dir / "samplesheet.csv"
        with samplesheet.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample_id", "fastq_1", "fastq_2", "single_end", "strandedness", "platform"])
            writer.writeheader()
            writer.writerows(rows)
        return samplesheet

    @staticmethod
    def _sample_id_from_fastq_name(name: str, read_token: str | None = None) -> str:
        sample = name
        for suffix in (".trimmed.fastq.gz", ".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if sample.endswith(suffix):
                sample = sample[: -len(suffix)]
                break
        if read_token and sample.endswith(read_token):
            sample = sample[: -len(read_token)]
        return sample.rstrip("._-") or "uploaded_sample"

    def _fastqs_to_samplesheet_rows(self, fastqs: list[Path]) -> list[dict[str, str]]:
        r1_by_sample = {}
        r2_by_sample = {}
        singles = []
        for path in sorted(fastqs, key=lambda item: item.name):
            name = path.name
            if "_R1" in name:
                r1_by_sample[self._sample_id_from_fastq_name(name, "_R1")] = path
            elif "_R2" in name:
                r2_by_sample[self._sample_id_from_fastq_name(name, "_R2")] = path
            elif "_1" in name:
                r1_by_sample[self._sample_id_from_fastq_name(name, "_1")] = path
            elif "_2" in name:
                r2_by_sample[self._sample_id_from_fastq_name(name, "_2")] = path
            else:
                singles.append(path)
        rows = []
        if r1_by_sample or r2_by_sample:
            for sample_id in sorted(set(r1_by_sample) | set(r2_by_sample)):
                r1 = r1_by_sample.get(sample_id)
                r2 = r2_by_sample.get(sample_id)
                if not r1 or not r2:
                    missing = "R1" if not r1 else "R2"
                    raise ValueError(f"Missing {missing} FASTQ for paired sample {sample_id}.")
                rows.append({
                    "sample_id": sample_id,
                    "fastq_1": str(r1),
                    "fastq_2": str(r2),
                    "single_end": "false",
                    "strandedness": "unknown",
                    "platform": "unknown",
                })
        for index, path in enumerate(singles, start=1):
            rows.append({
                "sample_id": self._sample_id_from_fastq_name(path.name) or f"uploaded_sample_{index}",
                "fastq_1": str(path),
                "fastq_2": "",
                "single_end": "true",
                "strandedness": "unknown",
                "platform": "unknown",
            })
        return rows

    def selected_steps_for_workflow(self, workflow_id: str) -> list[str]:
        config_path = self.app_root / "config" / "workflow_steps.yaml"
        for workflow in self.load_steps(config_path):
            if workflow.get("id") == workflow_id:
                return list(workflow.get("selected_steps") or [])
        raise ValueError(f"Unknown workflow id: {workflow_id}")

    @staticmethod
    def needs_fastq_input(selected_steps: list[str]) -> bool:
        read_steps = {
            "rnaseq_03_raw_read_qc",
            "rnaseq_04_adapter_quality_trimming",
            "rnaseq_05_post_trim_quality_control",
            "rnaseq_06_strandedness_inference",
            "rnaseq_07_salmon_quantification",
            "rnaseq_09_star_alignment",
            "rnaseq_09d_hisat2_alignment",
            "genomics_01_input_validation",
            "genomics_03_raw_qc",
            "genomics_04_trim_fastq",
        }
        return bool(read_steps & set(selected_steps))

    def trim_manifest_options(
        self,
        session_id: str,
        tester_id: str,
        omics_type: str,
        project_id: str,
    ) -> list[dict[str, str]]:
        options = []
        seen = set()
        for job in self.job_store.list_workspace_jobs(session_id, tester_id, omics_type, project_id=project_id, limit=500):
            results_dir = job.get("results_dir")
            if not results_dir:
                continue
            for name in ("trim_manifest.tsv", "trimmed_fastq_manifest.tsv"):
                manifest = Path(results_dir) / "04_trimmed" / name
                if not manifest.exists() or str(manifest) in seen:
                    continue
                seen.add(str(manifest))
                label = f"{job.get('workflow_id', 'run')} / {job.get('run_id', manifest.parent.name)} / {name}"
                options.append({"label": label, "value": str(manifest)})
        return options

    @staticmethod
    def validate_trim_manifest(manifest_path: str) -> list[str]:
        errors = []
        path = Path(manifest_path).resolve()
        if not path.exists() or not path.is_file():
            return [f"Trim manifest does not exist: {manifest_path}"]
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fieldnames = set(reader.fieldnames or [])
            required = {"sample_id", "fastq_1", "fastq_2", "single_end", "strandedness"}
            legacy_required = {"sample_id", "library_layout", "fastq_1", "fastq_2"}
            if not required.issubset(fieldnames) and not legacy_required.issubset(fieldnames):
                return [
                    "Trim manifest must contain columns: "
                    f"{', '.join(sorted(required))}. "
                    "Legacy manifests with sample_id, library_layout, fastq_1, fastq_2 are also accepted."
                ]
            rows = list(reader)
        if not rows:
            errors.append("Trim manifest has no samples.")
        manifest_dir = path.parent
        seen = set()
        for index, row in enumerate(rows, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            if not sample_id:
                errors.append(f"Trim manifest row {index} is missing sample_id.")
                continue
            if sample_id in seen:
                errors.append(f"Duplicate sample_id in trim manifest: {sample_id}")
            seen.add(sample_id)
            single_end_raw = str(row.get("single_end") or "").strip().lower()
            layout = (row.get("library_layout") or ("single" if single_end_raw in {"true", "1", "yes"} else "paired")).strip().lower()
            fastq_1 = (row.get("fastq_1") or "").strip()
            fastq_2 = (row.get("fastq_2") or "").strip()
            if layout not in {"single", "paired"}:
                errors.append(f"{sample_id}: single_end must be true/false or library_layout must be single/paired.")
            if not fastq_1:
                errors.append(f"{sample_id}: fastq_1 is missing.")
            if layout == "paired" and not fastq_2:
                errors.append(f"{sample_id}: paired-end row is missing fastq_2.")
            for value, column in ((fastq_1, "fastq_1"), (fastq_2, "fastq_2")):
                if not value:
                    continue
                fastq_path = Path(value)
                if not fastq_path.is_absolute():
                    fastq_path = manifest_dir / fastq_path
                if not fastq_path.exists() or not fastq_path.is_file():
                    errors.append(f"{sample_id}: {column} file does not exist: {fastq_path}")
        return errors

    @staticmethod
    def validate_fastq_pairs(fastq_files: list[str]) -> list[str]:
        errors = []
        if not fastq_files:
            return ["Select at least one trimmed FASTQ file."]
        r1_samples = {}
        r2_samples = {}
        singles = []
        for raw_path in fastq_files:
            name = Path(raw_path).name
            sample = name
            for token in (".trimmed.fastq.gz", ".fastq.gz", ".fq.gz", ".fastq", ".fq"):
                sample = sample.replace(token, "")
            if "_R1" in sample:
                r1_samples[sample.replace("_R1", "")] = raw_path
            elif "_R2" in sample:
                r2_samples[sample.replace("_R2", "")] = raw_path
            elif "_1" in sample:
                r1_samples[sample.replace("_1", "")] = raw_path
            elif "_2" in sample:
                r2_samples[sample.replace("_2", "")] = raw_path
            else:
                singles.append(raw_path)
        if r1_samples or r2_samples:
            missing_r2 = sorted(set(r1_samples) - set(r2_samples))
            missing_r1 = sorted(set(r2_samples) - set(r1_samples))
            if missing_r2:
                errors.append(f"Missing R2 trimmed FASTQ for sample(s): {', '.join(missing_r2)}")
            if missing_r1:
                errors.append(f"Missing R1 trimmed FASTQ for sample(s): {', '.join(missing_r1)}")
        if singles and (r1_samples or r2_samples):
            errors.append("Do not mix paired-end R1/R2 files and single-end FASTQ files in the same strandness-only run.")
        return errors

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

    def default_genomics_samplesheet(self, run_dir: Path) -> Path:
        fastq_dir = self.pipeline_root.parent / "testdatasets" / "test_data" / "human_chr22_genomics" / "fastq"
        r1 = fastq_dir / "genomics_test_R1.fastq.gz"
        r2 = fastq_dir / "genomics_test_R2.fastq.gz"
        if not r1.exists() or not r2.exists():
            raise FileNotFoundError(f"Default genomics FASTQ demo files are missing under {fastq_dir}")
        samplesheet = run_dir / "demo_genomics_samplesheet.csv"
        samplesheet.write_text(
            "sample_id,fastq_1,fastq_2,single_end\n"
            f"genomics_test,{r1},{r2},false\n",
            encoding="utf-8",
        )
        return samplesheet

    def default_genomics_reference_fasta(self) -> str:
        fasta = self.pipeline_root.parent / "testdatasets" / "test_data" / "human_chr22_rnaseq" / "refs" / "chr22_with_ERCC92.fa"
        if not fasta.exists():
            raise FileNotFoundError(f"Default genomics reference FASTA is missing: {fasta}")
        return str(fasta)

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
        selected_steps = self.selected_steps_for_workflow(workflow_id)
        trim_manifest = self._optional_param(params.get("trim_manifest"))
        if not self.needs_fastq_input(selected_steps):
            samplesheet = run_dir / "no_fastq_input.csv"
            samplesheet.write_text("sample_id,fastq_1,fastq_2,single_end,strandedness,platform\n", encoding="utf-8")
        elif params.get("trim_input_mode") == "previous_manifest" and trim_manifest:
            manifest_errors = self.validate_trim_manifest(trim_manifest)
            if manifest_errors:
                raise ValueError("; ".join(manifest_errors))
            samplesheet = Path(trim_manifest).resolve()
        else:
            try:
                samplesheet = self.create_or_find_samplesheet(sid, run_dir, tester, omics, project, selected_files)
            except ValueError:
                if omics != "genomics":
                    raise
                samplesheet = self.default_genomics_samplesheet(run_dir)

        outdir = run_dir / "results"
        demo_reference = self.pipeline_root.parent / "refs" / "gallus_gallus_ensembl116" / "mini_ref"
        reference_mode = params.get("reference_mode", "demo_reference")
        selected_route = params.get("selected_route", "salmon")
        if workflow_id in {"salmon_quant_only", "salmon_count_matrix", "qc_trim_strandedness"}:
            selected_route = "salmon"
        elif workflow_id in {"star_align_only", "star_count_matrix", "star_htseq_route"}:
            selected_route = "star"
        elif workflow_id in {"hisat2_align_only", "hisat2_alignment_route", "hisat2_featurecounts_route", "hisat2_htseq_route"}:
            selected_route = "hisat2"
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
            "reference_mode": reference_mode,
            "reference_manifest": self._optional_param(params.get("reference_manifest")) or self._optional_param(params.get("reference_bundle")),
            "selected_route": selected_route,
            "organism": self._optional_param(params.get("organism")) or ("Gallus_gallus" if reference_mode == "demo_reference" else "demo"),
            "genome_build": self._optional_param(params.get("genome_build")) or ("Ensembl_116_mini_GRCg7b" if reference_mode == "demo_reference" else "demo_build"),
            "genome_fasta": self._optional_param(params.get("genome_fasta")) or (str(demo_reference / "Gallus_gallus.mini.fa.gz") if reference_mode == "demo_reference" else None),
            "gtf": self._optional_param(params.get("gtf")) or (str(demo_reference / "Gallus_gallus.mini.gtf.gz") if reference_mode == "demo_reference" else None),
            "transcriptome_fasta": self._optional_param(params.get("transcriptome_fasta")) or (str(demo_reference / "Gallus_gallus.mini.fa.gz") if reference_mode == "demo_reference" else None),
            "tx2gene": self._optional_param(params.get("tx2gene")) or (str(demo_reference / "tx2gene.tsv") if reference_mode == "demo_reference" else None),
            "strandedness_method": params.get("strandedness_method", "salmon_auto"),
            "salmon_index": self._optional_param(params.get("salmon_index")),
            "strandedness_inference_reads": int(params.get("strandedness_inference_reads", 1000000) or 1000000),
            "existing_bam": self._optional_param(params.get("existing_bam")),
            "star_index": self._optional_param(params.get("star_index")),
            "hisat2_index": self._optional_param(params.get("hisat2_index")),
            "rseqc_ref_bed": self._optional_param(params.get("rseqc_ref_bed")),
            "rseqc_stranded_threshold": float(params.get("rseqc_stranded_threshold", 0.6) or 0.6),
            "manual_strandedness": self._optional_param(params.get("manual_strandedness")),
            "strandedness_approval_status": params.get("strandedness_approval_status", "pending"),
            "trim_input_mode": params.get("trim_input_mode", "manual_trimmed_fastq"),
            "trim_manifest": trim_manifest,
            "threads": int(params.get("threads", 4) or 4),
        }
        if omics == "genomics":
            genomics_reference = (
                self.default_genomics_reference_fasta()
                if reference_mode == "demo_reference"
                else self._optional_param(params.get("genome_fasta")) or self.default_genomics_reference_fasta()
            )
            resolved_params = {
                "workflow_type": params.get("workflow_type", "dna_short_read"),
                "genome_fasta": genomics_reference,
                "reference_name": self._optional_param(params.get("reference_name")) or self._optional_param(params.get("genome_build")) or "human_chr22_demo",
                "workflow_name": workflow_id,
                "quality_cutoff": int(params.get("quality_threshold", params.get("quality_cutoff", 20)) or 20),
                "min_length": int(params.get("minimum_read_length", params.get("min_length", 20)) or 20),
                "threads": int(params.get("threads", 2) or 2),
            }
        elif "rnaseq_04_adapter_quality_trimming" not in selected_steps:
            if (
                "rnaseq_06_strandedness_inference" not in selected_steps
                and "rnaseq_06a_reference_build_validation" not in selected_steps
                and "rnaseq_07_salmon_quantification" not in selected_steps
                and "rnaseq_09_star_alignment" not in selected_steps
                and "rnaseq_09d_hisat2_alignment" not in selected_steps
            ):
                resolved_params = {
                    "qc_tool": resolved_params["qc_tool"],
                    "threads": resolved_params["threads"],
                }
        run_request = {
            "run_id": run_id,
            "omics": "genomics" if omics == "genomics" else "rnaseq",
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
        if job.get("compiled_dir") and str(job.get("status", "")).lower() == "running":
            try:
                subprocess.run(
                    ["python3", "-m", "builder.api", "collect-outputs", "--compiled-dir", job["compiled_dir"], "--status", "running"],
                    cwd=self.pipeline_root,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    check=False,
                    timeout=10,
                )
            except subprocess.TimeoutExpired:
                pass
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

    def execution_record(self, job_id: str, refresh: bool = True) -> dict[str, Any]:
        job = self.refresh_job(job_id) if refresh else self.job_store.get_job(job_id)
        if not job or not job.get("compiled_dir"):
            return {}
        record = Path(job["compiled_dir"]) / "execution_record.json"
        if record.exists():
            return json.loads(record.read_text())
        plan = Path(job["compiled_dir"]) / "workflow_plan.json"
        if plan.exists():
            return {
                "status": job.get("status", "unknown"),
                "workflow_plan_path": str(plan),
                "step_outputs": [
                    {
                        "step_id": step.get("step_id"),
                        "step_name": step.get("step_name"),
                        "process_name": step.get("process_name"),
                        "status": "pending",
                        "outputs": [],
                        "missing_outputs": [],
                        "tasks": [],
                        "dependencies": step.get("dependencies", []),
                    }
                    for step in json.loads(plan.read_text()).get("steps", [])
                ],
            }
        return {}
