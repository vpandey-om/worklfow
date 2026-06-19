import csv
from pathlib import Path
from typing import Any

from .models import ParameterSpec, StepSpec, StepSuggestion, ValidationResult


class WorkflowPlanner:
    def __init__(self, steps: dict[str, StepSpec]):
        self.steps = steps

    def get_steps_for_omics(self, omics: str) -> list[StepSpec]:
        return [s for s in self.steps.values() if omics in s.omics]

    def get_available_first_steps(self, omics: str) -> list[StepSpec]:
        candidates = self.get_steps_for_omics(omics)
        return [s for s in candidates if s.can_run_individually]

    def suggest_next_steps(self, selected_steps: list[str], omics: str) -> dict[str, list[StepSuggestion]]:
        if not selected_steps:
            first = self.get_available_first_steps(omics)
            return {
                "recommended": [
                    self._suggestion(step, True, "This step can start a workflow for the selected omics type.")
                    for step in first[:5]
                ],
                "allowed": [self._suggestion(step, False, "Allowed first step.") for step in first],
            }

        last_id = selected_steps[-1]
        if last_id not in self.steps:
            return {"recommended": [], "allowed": []}
        last = self.steps[last_id]
        recommended_ids = last.next_steps.get("recommended", [])
        allowed_ids = last.next_steps.get("allowed", [])
        recommended = [
            self._suggestion(self.steps[sid], True, self.steps[sid].purpose)
            for sid in recommended_ids
            if sid in self.steps and omics in self.steps[sid].omics
        ]
        allowed = [
            self._suggestion(self.steps[sid], sid in recommended_ids, self.steps[sid].purpose)
            for sid in allowed_ids
            if sid in self.steps and omics in self.steps[sid].omics
        ]
        return {"recommended": recommended, "allowed": allowed}

    def validate_selected_steps(self, selected_steps: list[str], omics: str) -> ValidationResult:
        errors = []
        warnings = []
        suggested_fixes = []
        seen = set()

        for idx, step_id in enumerate(selected_steps):
            if step_id not in self.steps:
                errors.append({"step_id": step_id, "message": "Unknown step_id"})
                suggested_fixes.append("Choose a step_id present in registry/steps.yaml.")
                continue

            step = self.steps[step_id]
            if omics not in step.omics:
                errors.append({"step_id": step_id, "message": f"Step is not compatible with omics={omics}"})

            if step_id in seen:
                errors.append({"step_id": step_id, "message": "Duplicate step in workflow"})
            seen.add(step_id)

            if idx > 0:
                prev_id = selected_steps[idx - 1]
                prev = self.steps.get(prev_id)
                if prev and step_id not in prev.next_steps.get("allowed", []) and step_id not in prev.next_steps.get("recommended", []):
                    warnings.append({
                        "step_id": step_id,
                        "message": f"{step_id} is not listed as a normal next step after {prev_id}; verify input/output compatibility.",
                    })

            missing_deps = [d for d in step.depends_on if d not in selected_steps[:idx]]
            if missing_deps and not step.can_run_individually:
                errors.append({"step_id": step_id, "message": f"Missing dependencies: {missing_deps}"})
                suggested_fixes.extend([f"Add {dep} before {step_id}." for dep in missing_deps])
            elif missing_deps and idx > 0:
                warnings.append({
                    "step_id": step_id,
                    "message": f"Registry dependencies not selected: {missing_deps}. Allowed because the step can run individually with a FASTQ manifest.",
                })

        return ValidationResult(
            valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
            suggested_fixes=suggested_fixes,
        )

    def resolve_effective_params(self, selected_steps: list[str], user_params: dict) -> dict:
        user_params = user_params or {}
        resolved = {}
        for step_id in selected_steps:
            step = self.steps[step_id]
            for name, spec in step.parameters.items():
                resolved.setdefault(name, spec.default)
            step_params = user_params.get(step_id)
            if isinstance(step_params, dict):
                resolved.update(step_params)
        resolved.update({key: value for key, value in user_params.items() if not isinstance(value, dict)})
        if "cutadapt_minimum_overlap" in resolved:
            resolved["cutadapt_min_overlap"] = resolved["cutadapt_minimum_overlap"]
        elif "cutadapt_min_overlap" in resolved:
            resolved["cutadapt_minimum_overlap"] = resolved["cutadapt_min_overlap"]
        return resolved

    def validate_params(
        self,
        selected_steps: list[str],
        params: dict[str, Any],
        input_path: str | Path | None = None,
    ) -> ValidationResult:
        errors = []
        warnings = []
        suggested_fixes = []
        selected_param_specs = self._selected_param_specs(selected_steps)

        for name, spec in selected_param_specs.items():
            if name not in params:
                continue
            value = params.get(name)
            if value is None or value == "":
                continue
            if spec.choices and value not in spec.choices and str(value).lower() not in {str(choice).lower() for choice in spec.choices}:
                errors.append({"param": name, "message": f"Value must be one of {spec.choices}"})
            if spec.min is not None or spec.max is not None:
                numeric = self._as_number(value)
                if numeric is None:
                    errors.append({"param": name, "message": "Value must be numeric"})
                else:
                    if spec.min is not None and numeric < spec.min:
                        errors.append({"param": name, "message": f"Value must be >= {spec.min}"})
                    if spec.max is not None and numeric > spec.max:
                        errors.append({"param": name, "message": f"Value must be <= {spec.max}"})

        if "rnaseq_04_adapter_quality_trimming" in selected_steps:
            trim_poly_g = params.get("trim_poly_g")
            if trim_poly_g not in {"auto", True, False, "true", "false", "on", "off"}:
                errors.append({"param": "trim_poly_g", "message": "Value must be auto, true, or false"})

            if params.get("trimming_tool") == "cutadapt":
                if not params.get("adapter_sequence_r1"):
                    errors.append({
                        "param": "adapter_sequence_r1",
                        "message": "Cutadapt requires an R1 adapter sequence.",
                    })
                if self._input_may_be_paired(input_path) and not params.get("adapter_sequence_r2"):
                    errors.append({
                        "param": "adapter_sequence_r2",
                        "message": "Cutadapt requires an R2 adapter sequence for paired-end data.",
                    })
                if errors:
                    suggested_fixes.append(
                        "Use fastp, or provide explicit adapter sequences before selecting Cutadapt."
                    )

        if "rnaseq_06_strandedness_inference" in selected_steps:
            method = params.get("strandedness_method", "salmon_auto")
            if method not in {"salmon_auto", "rseqc_validation"}:
                errors.append({
                    "param": "strandedness_method",
                    "message": "Strandedness method must be salmon_auto or rseqc_validation.",
                })
            can_build_salmon_index = (
                "rnaseq_06a_reference_build_validation" in selected_steps
                and bool(params.get("transcriptome_fasta"))
            )
            if method == "salmon_auto" and not params.get("salmon_index") and not can_build_salmon_index:
                errors.append({
                    "param": "salmon_index",
                    "message": "Salmon strandedness inference requires a Salmon transcriptome index path, or transcriptome_fasta so the reference step can build one.",
                })
                suggested_fixes.append("Provide salmon_index, or include reference build with transcriptome_fasta before running strandedness inference.")
            if method == "rseqc_validation":
                if not params.get("rseqc_ref_bed"):
                    errors.append({
                        "param": "rseqc_ref_bed",
                        "message": "RSeQC infer_experiment.py requires a BED annotation file.",
                    })
                if not params.get("existing_bam") and not params.get("star_index"):
                    errors.append({
                        "param": "star_index",
                        "message": "RSeQC fallback requires either an existing sorted BAM or a STAR genome index to create one.",
                    })
                suggested_fixes.append("Use Salmon auto by default, or provide existing_bam/rseqc_ref_bed, or star_index/rseqc_ref_bed for alignment-based validation.")

        if "rnaseq_06a_reference_build_validation" in selected_steps:
            route = params.get("selected_route", "salmon")
            mode = params.get("reference_mode", "demo_reference")
            if route not in {"salmon", "star", "hisat2", "custom", "salmon_de", "star_de", "quant_only", "alignment_only", "hisat2_align"}:
                errors.append({"param": "selected_route", "message": "Route must be salmon, star, hisat2, custom, salmon_de, star_de, quant_only, alignment_only, or hisat2_align."})
            if not params.get("organism"):
                errors.append({"param": "organism", "message": "Reference organism metadata is required."})
            if not params.get("genome_build"):
                errors.append({"param": "genome_build", "message": "Reference genome_build metadata is required."})
            salmon_route = route in {"salmon", "salmon_de", "quant_only"} or "rnaseq_07_salmon_quantification" in selected_steps
            star_route = route in {"star", "star_de", "alignment_only"} or "rnaseq_09_star_alignment" in selected_steps
            hisat2_route = route in {"hisat2", "hisat2_align"} or "rnaseq_09d_hisat2_alignment" in selected_steps
            if mode != "demo_reference":
                if salmon_route and not (params.get("transcriptome_fasta") and params.get("tx2gene")):
                    errors.append({"param": "transcriptome_fasta", "message": "Salmon route requires transcriptome_fasta and tx2gene; salmon_index is optional if it can be built."})
                if star_route and not (params.get("genome_fasta") and params.get("gtf")):
                    errors.append({"param": "genome_fasta", "message": "STAR route requires genome_fasta and gtf; star_index is optional if it can be built."})
                if hisat2_route and not params.get("genome_fasta"):
                    errors.append({"param": "genome_fasta", "message": "HISAT2 route requires genome_fasta; hisat2_index is optional if it can be built."})

        return ValidationResult(
            valid=len(errors) == 0,
            errors=errors,
            warnings=warnings,
            suggested_fixes=suggested_fixes,
        )

    def _selected_param_specs(self, selected_steps: list[str]) -> dict[str, ParameterSpec]:
        specs = {}
        for step_id in selected_steps:
            step = self.steps.get(step_id)
            if step:
                specs.update(step.parameters)
        return specs

    @staticmethod
    def _as_number(value: Any) -> float | None:
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _input_may_be_paired(input_path: str | Path | None) -> bool:
        if not input_path:
            return True
        path = Path(input_path)
        if not path.exists():
            return True
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                layout = (row.get("library_layout") or "").lower()
                if layout == "paired" or row.get("fastq_2"):
                    return True
        return False

    @staticmethod
    def _suggestion(step: StepSpec, recommended: bool, reason: str) -> StepSuggestion:
        return StepSuggestion(
            step_id=step.step_id,
            label=step.step_name,
            reason=reason,
            recommended=recommended,
        )


def get_available_first_steps(omics: str, steps: dict[str, StepSpec]) -> list[StepSpec]:
    return WorkflowPlanner(steps).get_available_first_steps(omics)


def suggest_next_steps(selected_steps: list[str], omics: str, steps: dict[str, StepSpec]) -> dict[str, list[StepSuggestion]]:
    return WorkflowPlanner(steps).suggest_next_steps(selected_steps, omics)


def validate_selected_steps(selected_steps: list[str], omics: str, steps: dict[str, StepSpec]) -> ValidationResult:
    return WorkflowPlanner(steps).validate_selected_steps(selected_steps, omics)


def resolve_effective_params(selected_steps: list[str], user_params: dict, steps: dict[str, StepSpec]) -> dict:
    return WorkflowPlanner(steps).resolve_effective_params(selected_steps, user_params)
