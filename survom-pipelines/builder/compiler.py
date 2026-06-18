import json
import csv
from datetime import datetime, timezone
from pathlib import Path

import yaml
from jinja2 import Environment, FileSystemLoader

from .models import RunRequest, StepSpec
from .planner import WorkflowPlanner


class NextflowCompiler:
    def __init__(self, registry_steps: dict[str, StepSpec], template_dir: str | Path):
        self.steps = registry_steps
        self.planner = WorkflowPlanner(registry_steps)
        self.template_dir = Path(template_dir)
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
        param_validation = self.planner.validate_params(request.selected_steps, params, request.input)
        if not param_validation.valid:
            raise ValueError(param_validation.model_dump_json(indent=2))

        selected = [self.steps[sid] for sid in request.selected_steps]
        params["input"] = request.input
        params["outdir"] = request.outdir
        params["run_id"] = request.run_id
        params["execution_profile"] = request.execution_profile

        main_template = self.env.get_template("main.nf.j2")
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
        (out / "params.yaml").write_text(params_template.render(params=params), encoding="utf-8")
        command = self._nextflow_command(out, request.execution_profile)
        run_command_path = out / "run_command.sh"
        run_command_path.write_text("#!/usr/bin/env bash\nset -euo pipefail\n\n" + command + "\n", encoding="utf-8")
        run_command_path.chmod(0o755)

        input_summary = self._input_summary(request.input)
        step_records = self._step_records(selected, params)
        methods_text = self._methods_text(selected, params, input_summary)
        (out / "methods.md").write_text(methods_text, encoding="utf-8")

        ui_summary = {
            "run_id": request.run_id,
            "run_label": request.run_label or request.run_id,
            "run_type": "single_step" if len(request.selected_steps) == 1 else "workflow_chain",
            "omics": request.omics,
            "user_id": request.user_id,
            "dataset_id": request.dataset_id,
            "input": input_summary,
            "selected_steps": step_records,
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
                "methods.md",
                "run_command.sh",
            ],
        }
        (out / "workflow_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        return out

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
            f"nextflow run {out / 'main.nf'} "
            f"-profile {profile} "
            f"-params-file {out / 'params.yaml'} "
            "-resume"
        )

    @staticmethod
    def _input_summary(input_path: str) -> dict:
        path = Path(input_path)
        summary = {
            "samplesheet": input_path,
            "exists_at_compile_time": path.exists(),
            "sample_count": None,
            "columns": [],
            "fastq_columns": [],
            "samples": [],
        }
        if not path.exists():
            return summary
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
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
