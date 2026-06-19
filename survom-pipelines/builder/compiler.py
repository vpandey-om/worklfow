import json
import csv
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import yaml
from jinja2 import Environment, FileSystemLoader

from .models import RunRequest, StepSpec
from .planner import WorkflowPlanner


class NextflowCompiler:
    def __init__(self, registry_steps: dict[str, StepSpec], template_dir: str | Path):
        self.steps = registry_steps
        self.planner = WorkflowPlanner(registry_steps)
        self.template_dir = Path(template_dir)
        self.project_root = self.template_dir.parents[1]
        self.env = Environment(
            loader=FileSystemLoader(str(self.template_dir)),
            trim_blocks=True,
            lstrip_blocks=True,
        )

    def compile(self, request: RunRequest, output_dir: str | Path) -> Path:
        validation = self.planner.validate_selected_steps(request.selected_steps, request.omics)
        if not validation.valid:
            raise ValueError(validation.model_dump_json(indent=2))

        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        params = self.planner.resolve_effective_params(request.selected_steps, request.params)
        params = self._resolve_reference_manifest_params(params)
        param_validation = self.planner.validate_params(request.selected_steps, params, request.input)
        if not param_validation.valid:
            raise ValueError(param_validation.model_dump_json(indent=2))

        selected = [self.steps[sid] for sid in request.selected_steps]
        resolved_input, input_summary = self._prepare_samplesheet(request, out)
        params["input"] = str(resolved_input.resolve())
        params["outdir"] = request.outdir
        params["run_id"] = request.run_id
        params["execution_profile"] = request.execution_profile
        params["needs_fastq_channel"] = self._needs_fastq_channel(request.selected_steps)

        main_template = self.env.get_template("main.nf.j2")
        atomic_template = self.env.get_template("atomic_step.nf.j2")
        config_template = self.env.get_template("nextflow.config.j2")
        params_template = self.env.get_template("params.yaml.j2")

        (out / "main.nf").write_text(
            main_template.render(request=request, steps=selected, params=params),
            encoding="utf-8",
        )
        (out / "nextflow.config").write_text(
            config_template.render(request=request, steps=selected, params=params),
            encoding="utf-8",
        )
        atomic_workflows = []
        for step in selected:
            atomic_request = SimpleNamespace(**request.model_dump())
            atomic_request.selected_steps = [step.step_id]
            atomic_params = dict(params)
            atomic_params["needs_fastq_channel"] = self._needs_fastq_channel([step.step_id])
            path = out / f"atomic_{step.step_id}.nf"
            path.write_text(
                atomic_template.render(request=atomic_request, step=step, params=atomic_params),
                encoding="utf-8",
            )
            config_path = out / f"atomic_{step.step_id}.config"
            config_path.write_text(
                config_template.render(request=atomic_request, steps=[step], params=atomic_params),
                encoding="utf-8",
            )
            atomic_workflows.append({
                "step_id": step.step_id,
                "process_name": step.process_name,
                "path": str(path),
                "relative_path": path.name,
                "config_path": str(config_path),
                "config_relative_path": config_path.name,
                "command": (
                    f"cd {out} && nextflow -C {config_path.name} run {path.name} "
                    f"-params-file params.yaml -work-dir {out / ('work_' + step.step_id)} -resume"
                ),
            })
        (out / "params.yaml").write_text(params_template.render(params=params), encoding="utf-8")
        command = self._nextflow_command(out, request.execution_profile)
        run_command_path = out / "run_command.sh"
        run_command_path.write_text("#!/usr/bin/env bash\nset -euo pipefail\n\n" + command + "\n", encoding="utf-8")
        run_command_path.chmod(0o755)

        step_records = self._step_records(selected, params)
        workflow_plan = self._workflow_plan(selected, params, input_summary)
        methods_text = self._methods_text(selected, params, input_summary)
        (out / "methods.md").write_text(methods_text, encoding="utf-8")
        (out / "workflow_plan.json").write_text(json.dumps(workflow_plan, indent=2), encoding="utf-8")

        ui_summary = {
            "run_id": request.run_id,
            "run_label": request.run_label or request.run_id,
            "run_type": "single_step" if len(request.selected_steps) == 1 else "workflow_chain",
            "omics": request.omics,
            "user_id": request.user_id,
            "dataset_id": request.dataset_id,
            "input": input_summary,
            "selected_steps": step_records,
            "workflow_plan": workflow_plan,
            "atomic_workflows": atomic_workflows,
            "outdir": request.outdir,
            "status": "compiled",
            "display_tabs": self._display_tabs(step_records),
            "nextflow_command": command,
        }
        (out / "ui_run_summary.json").write_text(json.dumps(ui_summary, indent=2), encoding="utf-8")

        manifest = {
            "schema_version": "1.1",
            "run_id": request.run_id,
            "run_label": request.run_label or request.run_id,
            "run_description": request.run_description,
            "run_type": "single_step" if len(request.selected_steps) == 1 else "workflow_chain",
            "compiled_at": datetime.now(timezone.utc).isoformat(),
            "omics": request.omics,
            "user_id": request.user_id,
            "dataset_id": request.dataset_id,
            "execution_profile": request.execution_profile,
            "input": input_summary,
            "selected_steps": request.selected_steps,
            "step_records": step_records,
            "workflow_plan": workflow_plan,
            "atomic_workflows": atomic_workflows,
            "params": params,
            "validation": validation.model_dump(),
            "param_validation": param_validation.model_dump(),
            "nextflow_command": command,
            "generated_files": [
                "main.nf",
                "nextflow.config",
                "params.yaml",
                "workflow_manifest.json",
                "ui_run_summary.json",
                "workflow_plan.json",
                *[item["relative_path"] for item in atomic_workflows],
                *[item["config_relative_path"] for item in atomic_workflows],
                "methods.md",
                "run_command.sh",
            ],
        }
        (out / "workflow_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        return out

    @staticmethod
    def _needs_fastq_channel(selected_steps: list[str]) -> bool:
        read_steps = {
            "rnaseq_03_raw_read_qc",
            "rnaseq_04_adapter_quality_trimming",
            "rnaseq_05_post_trim_quality_control",
            "rnaseq_06_strandedness_inference",
            "rnaseq_07_salmon_quantification",
            "rnaseq_09_star_alignment",
            "rnaseq_09d_hisat2_alignment",
        }
        return bool(read_steps & set(selected_steps))

    @staticmethod
    def _resolve_reference_manifest_params(params: dict) -> dict:
        manifest = params.get("reference_manifest")
        if not manifest:
            return params
        path = Path(manifest).expanduser()
        if not path.exists():
            raise FileNotFoundError(f"Reference manifest does not exist: {manifest}")
        if path.is_dir():
            for candidate_name in ("reference_manifest.tsv", "validated_reference_bundle.json"):
                candidate = path / candidate_name
                if candidate.exists() and candidate.is_file():
                    path = candidate
                    break
            else:
                raise FileNotFoundError(
                    f"Reference bundle directory must contain reference_manifest.tsv or validated_reference_bundle.json: {manifest}"
                )
        values = {}
        if path.suffix.lower() == ".json":
            loaded = json.loads(path.read_text())
            if isinstance(loaded, dict):
                values = loaded
        else:
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if {"key", "value"}.issubset(reader.fieldnames or []):
                    values = {row["key"]: row.get("value") for row in reader}
                else:
                    rows = list(reader)
                    if rows:
                        values = rows[0]
        aliases = {
            "organism": "species",
            "genome_build": "assembly",
            "gtf": "annotation_gtf",
            "transcriptome_fasta": "transcript_fasta",
        }
        for canonical, alias in aliases.items():
            if not values.get(canonical) and values.get(alias):
                values[canonical] = values[alias]
        resolved = dict(params)
        for key in (
            "reference_mode",
            "selected_route",
            "organism",
            "genome_build",
            "genome_fasta",
            "gtf",
            "transcriptome_fasta",
            "tx2gene",
            "salmon_index",
            "star_index",
            "hisat2_index",
        ):
            if values.get(key) and not resolved.get(key):
                resolved[key] = values[key]
        return resolved

    @staticmethod
    def _nextflow_command(out: Path, execution_profile: str) -> str:
        profile_map = {
            "local": "local",
            "local_docker": "local,docker",
            "docker": "local,docker",
            "singularity": "local,singularity",
            "aws": "aws",
        }
        profile = profile_map.get(execution_profile, execution_profile)
        return (
            f"cd {out} && "
            f"nextflow run main.nf "
            f"-profile {profile} "
            f"-params-file params.yaml "
            f"-work-dir {out / 'work'} "
            "-resume"
        )

    def _prepare_samplesheet(self, request: RunRequest, out: Path) -> tuple[Path, dict]:
        source = Path(request.input)
        if not source.exists():
            raise FileNotFoundError(f"Input sample sheet does not exist: {source}")

        if self._should_use_existing_fastq_manifest(request):
            resolved_manifest = self._normalize_fastq_manifest(source, out / "input_manifest.resolved.tsv")
            return resolved_manifest, self._input_summary(str(resolved_manifest), original_input=str(source))

        selected_fastqs = {
            Path(path).name: Path(path).resolve()
            for path in (request.selected_files or {}).get("fastq", [])
            if path
        }
        resolved_path = out / "samplesheet.resolved.tsv"
        rows = []
        with source.open(newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            if "fastq_1" not in fieldnames:
                raise ValueError("Input sample sheet must contain a fastq_1 column.")
            for line_no, row in enumerate(reader, start=2):
                sample_id = row.get("sample_id") or f"row_{line_no}"
                fastq_1 = str(self._resolve_fastq(row.get("fastq_1"), source.parent, selected_fastqs, sample_id, "fastq_1"))
                fastq_2 = (row.get("fastq_2") or "").strip()
                resolved_fastq_2 = ""
                if fastq_2:
                    resolved_fastq_2 = str(self._resolve_fastq(fastq_2, source.parent, selected_fastqs, sample_id, "fastq_2"))
                single_end = self._row_single_end(row, has_fastq_2=bool(resolved_fastq_2))
                rows.append({
                    "sample_id": sample_id,
                    "fastq_1": fastq_1,
                    "fastq_2": resolved_fastq_2,
                    "single_end": "true" if single_end else "false",
                    "strandedness": row.get("strandedness") or "unknown",
                    "platform": row.get("platform") or "unknown",
                })

        with resolved_path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["sample_id", "fastq_1", "fastq_2", "single_end", "strandedness", "platform"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)
        return resolved_path, self._input_summary(str(resolved_path), original_input=str(source))

    def _normalize_fastq_manifest(self, source: Path, target: Path) -> Path:
        delimiter = "\t" if source.suffix.lower() == ".tsv" else ","
        rows = []
        with source.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            fieldnames = set(reader.fieldnames or [])
            if "fastq_1" not in fieldnames:
                raise ValueError("Input manifest must contain a fastq_1 column.")
            for line_no, row in enumerate(reader, start=2):
                sample_id = row.get("sample_id") or f"row_{line_no}"
                fastq_1 = (row.get("fastq_1") or "").strip()
                fastq_2 = (row.get("fastq_2") or "").strip()
                if not fastq_1:
                    raise ValueError(f"Input manifest row {line_no} is missing fastq_1.")
                single_end = self._row_single_end(row, has_fastq_2=bool(fastq_2))
                if not single_end and not fastq_2:
                    raise ValueError(f"Input manifest row {line_no} is paired-end but missing fastq_2.")
                resolved_fastq_1 = str(self._resolve_manifest_fastq(fastq_1, source.parent, sample_id, "fastq_1"))
                resolved_fastq_2 = ""
                if fastq_2:
                    resolved_fastq_2 = str(self._resolve_manifest_fastq(fastq_2, source.parent, sample_id, "fastq_2"))
                rows.append({
                    "sample_id": sample_id,
                    "fastq_1": resolved_fastq_1,
                    "fastq_2": resolved_fastq_2,
                    "single_end": "true" if single_end else "false",
                    "strandedness": row.get("strandedness") or "unknown",
                    "platform": row.get("platform") or "unknown",
                })
        with target.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["sample_id", "fastq_1", "fastq_2", "single_end", "strandedness", "platform"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)
        return target

    @staticmethod
    def _should_use_existing_fastq_manifest(request: RunRequest) -> bool:
        if "rnaseq_04_adapter_quality_trimming" in request.selected_steps:
            return False
        manifest_consumers = {
            "rnaseq_05_post_trim_quality_control",
            "rnaseq_06_strandedness_inference",
            "rnaseq_07_salmon_quantification",
            "rnaseq_09_star_alignment",
            "rnaseq_09d_hisat2_alignment",
        }
        return bool(request.params.get("trim_manifest") or manifest_consumers & set(request.selected_steps))

    @staticmethod
    def _row_single_end(row: dict, has_fastq_2: bool) -> bool:
        raw_single_end = row.get("single_end")
        if raw_single_end is not None and str(raw_single_end).strip() != "":
            return str(raw_single_end).strip().lower() in {"true", "1", "yes", "single"}
        layout = str(row.get("library_layout") or "").strip().lower()
        if layout:
            return layout == "single"
        return not has_fastq_2

    def _resolve_fastq(
        self,
        raw_value: str | None,
        samplesheet_dir: Path,
        selected_fastqs: dict[str, Path],
        sample_id: str,
        column: str,
    ) -> Path:
        raw = (raw_value or "").strip()
        if not raw:
            raise FileNotFoundError(f"Missing {column} for sample '{sample_id}' in the sample sheet.")
        path = Path(raw).expanduser()
        candidates = [path] if path.is_absolute() else [
            samplesheet_dir / path,
            self.project_root / path,
            Path.cwd() / path,
        ]
        if path.name in selected_fastqs:
            candidates.append(selected_fastqs[path.name])
        for candidate in candidates:
            resolved = candidate.resolve()
            if resolved.exists() and resolved.is_file():
                return resolved
        checked = ", ".join(str(candidate) for candidate in candidates)
        raise FileNotFoundError(
            f"FASTQ file for sample '{sample_id}' column '{column}' does not exist or cannot be read: "
            f"{raw}. Checked: {checked}"
        )

    @staticmethod
    def _resolve_manifest_fastq(raw_value: str, manifest_dir: Path, sample_id: str, column: str) -> Path:
        path = Path(raw_value).expanduser()
        candidates = [path] if path.is_absolute() else [manifest_dir / path, Path.cwd() / path]
        for candidate in candidates:
            resolved = candidate.resolve()
            if resolved.exists() and resolved.is_file():
                return resolved
        checked = ", ".join(str(candidate) for candidate in candidates)
        raise FileNotFoundError(
            f"FASTQ file for sample '{sample_id}' column '{column}' does not exist or cannot be read: "
            f"{raw_value}. Checked: {checked}"
        )

    @staticmethod
    def _input_summary(input_path: str, original_input: str | None = None) -> dict:
        path = Path(input_path)
        summary = {
            "samplesheet": input_path,
            "original_samplesheet": original_input or input_path,
            "exists_at_compile_time": path.exists(),
            "sample_count": None,
            "columns": [],
            "fastq_columns": [],
            "samples": [],
        }
        if not path.exists():
            return summary
        delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter=delimiter)
            summary["columns"] = reader.fieldnames or []
            summary["fastq_columns"] = [c for c in summary["columns"] if c.startswith("fastq")]
            for row in reader:
                sample = {
                    "sample_id": row.get("sample_id"),
                    "library_layout": row.get("library_layout"),
                    "fastq_1": row.get("fastq_1"),
                    "fastq_2": row.get("fastq_2"),
                }
                summary["samples"].append(sample)
        summary["sample_count"] = len(summary["samples"])
        return summary

    @staticmethod
    def _step_records(steps: list[StepSpec], params: dict) -> list[dict]:
        records = []
        for index, step in enumerate(steps, start=1):
            step_param_names = list(step.parameters.keys())
            records.append({
                "order": index,
                "step_id": step.step_id,
                "step_name": step.step_name,
                "category": step.category,
                "purpose": step.purpose,
                "process_name": step.process_name,
                "module_path": step.module_path,
                "default_tool": step.default_tool,
                "fallback_tool": step.fallback_tool,
                "parameters": {name: params.get(name) for name in step_param_names},
                "expected_outputs": [output.model_dump() for output in step.outputs],
                "ui": step.ui,
                "chat_guidance": step.chat_guidance,
            })
        return records

    @staticmethod
    def _workflow_plan(steps: list[StepSpec], params: dict, input_summary: dict) -> dict:
        records = []
        previous = []
        for index, step in enumerate(steps, start=1):
            dependencies = list(step.depends_on or previous[-1:])
            records.append({
                "order": len(records) + 1,
                "step_id": step.step_id,
                "step_name": step.step_name,
                "process_name": step.process_name,
                "dependencies": dependencies,
                "inputs": [item.model_dump() for item in step.required_inputs],
                "outputs": [item.model_dump() for item in step.outputs],
                "status": "pending",
            })
            previous.append(step.step_id)
            if step.step_id == "rnaseq_04_adapter_quality_trimming":
                records.extend([
                    {
                        "order": len(records) + 1,
                        "step_id": "helper_make_preprocess_manifest",
                        "step_name": "Make preprocess manifest",
                        "process_name": "MAKE_PREPROCESS_MANIFEST",
                        "dependencies": [step.step_id],
                        "inputs": [{"name": "trimmed_fastq", "type": "file"}],
                        "outputs": [{"name": "manifest_fragment", "type": "tsv", "pattern": None}],
                        "status": "pending",
                    },
                    {
                        "order": len(records) + 2,
                        "step_id": "helper_merge_preprocess_manifest",
                        "step_name": "Merge preprocess manifest",
                        "process_name": "MERGE_PREPROCESS_MANIFEST",
                        "dependencies": ["helper_make_preprocess_manifest"],
                        "inputs": [{"name": "manifest_fragment", "type": "tsv"}],
                        "outputs": [{"name": "trimmed_fastq_manifest", "type": "tsv", "pattern": "04_trimmed/trimmed_fastq_manifest.tsv"}],
                        "status": "pending",
                    },
                ])
        return {
            "schema_version": "1.0",
            "input": input_summary,
            "steps": records,
        }

    @staticmethod
    def _display_tabs(step_records: list[dict]) -> list[dict]:
        tabs = [{"id": "overview", "label": "Overview"}, {"id": "inputs", "label": "Inputs"}]
        for step in step_records:
            tabs.append({
                "id": step["step_id"],
                "label": step["step_name"],
                "process_name": step["process_name"],
                "expected_outputs": step["expected_outputs"],
            })
        tabs.extend([
            {"id": "methods", "label": "Methods"},
            {"id": "logs", "label": "Logs"},
            {"id": "files", "label": "Files"},
        ])
        return tabs

    @staticmethod
    def _methods_text(steps: list[StepSpec], params: dict, input_summary: dict) -> str:
        lines = [
            "# Methods",
            "",
            f"Input sample sheet: `{input_summary['samplesheet']}`.",
        ]
        if input_summary.get("sample_count") is not None:
            lines.append(f"Number of samples: {input_summary['sample_count']}.")
        lines.append("")
        for step in steps:
            if step.step_id == "rnaseq_03_raw_read_qc":
                lines.extend([
                    "## Raw Read Quality Control",
                    "",
                    f"Raw FASTQ quality was assessed using {params.get('qc_tool', step.default_tool)}.",
                    "Per-sample HTML and ZIP FastQC reports were generated.",
                    "",
                ])
            elif step.step_id == "rnaseq_04_adapter_quality_trimming":
                lines.extend([
                    "## Adapter and Quality Trimming",
                    "",
                    (
                        f"Adapter and quality trimming was performed using {params.get('trimming_tool', step.default_tool)} "
                        f"with minimum base quality {params.get('quality_threshold')} and minimum read length "
                        f"{params.get('minimum_read_length')}."
                    ),
                    f"PolyG trimming mode was `{params.get('trim_poly_g')}` and polyX trimming was `{params.get('trim_poly_x')}`.",
                    "Trimmed FASTQ files and tool HTML/JSON reports were generated.",
                    "",
                ])
            else:
                lines.extend([
                    f"## {step.step_name}",
                    "",
                    step.purpose,
                    "",
                ])
        return "\n".join(lines)
