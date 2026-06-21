#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from builder.compiler import NextflowCompiler
from builder.models import RunRequest
from builder.registry import load_steps


@dataclass(frozen=True)
class Case:
    name: str
    selected_steps: list[str]


CASES = [
    Case("genomics_input_validation", ["genomics_01_input_validation"]),
    Case("genomics_reference_prepare", ["genomics_02_reference_prepare"]),
    Case("genomics_raw_qc", ["genomics_01_input_validation", "genomics_03_raw_qc"]),
    Case("genomics_trim_fastq", ["genomics_01_input_validation", "genomics_04_trim_fastq"]),
    Case("genomics_bwa_alignment", ["genomics_05_bwa_alignment"]),
    Case("genomics_bam_processing", ["genomics_06_bam_processing"]),
    Case("genomics_mark_duplicates", ["genomics_07_mark_duplicates"]),
    Case("genomics_coverage_qc", ["genomics_08_coverage_qc"]),
    Case("genomics_bqsr", ["genomics_09_bqsr"]),
    Case("genomics_germline_variants", ["genomics_10_germline_variants"]),
    Case("genomics_somatic_variants", ["genomics_11_somatic_variants"]),
    Case("genomics_variant_filtering", ["genomics_12_variant_filtering"]),
    Case("genomics_structural_variants", ["genomics_13_structural_variants"]),
    Case("genomics_cnv_calling", ["genomics_14_cnv_calling"]),
    Case("genomics_variant_annotation", ["genomics_15_variant_annotation"]),
    Case("genomics_multiqc_report", ["genomics_16_multiqc_report"]),
    Case("genomics_tool_inventory", ["genomics_17_tool_inventory"]),
    Case(
        "genomics_all_tool_inventory_chain",
        [
            "genomics_01_input_validation",
            "genomics_03_raw_qc",
            "genomics_04_trim_fastq",
            "genomics_17_tool_inventory",
            "genomics_99_provenance_manifest",
        ],
    ),
    Case(
        "genomics_full_demo_chain",
        [
            "genomics_01_input_validation",
            "genomics_02_reference_prepare",
            "genomics_03_raw_qc",
            "genomics_04_trim_fastq",
            "genomics_05_bwa_alignment",
            "genomics_06_bam_processing",
            "genomics_07_mark_duplicates",
            "genomics_08_coverage_qc",
            "genomics_09_bqsr",
            "genomics_10_germline_variants",
            "genomics_11_somatic_variants",
            "genomics_12_variant_filtering",
            "genomics_13_structural_variants",
            "genomics_14_cnv_calling",
            "genomics_15_variant_annotation",
            "genomics_16_multiqc_report",
            "genomics_99_provenance_manifest",
        ],
    ),
    Case(
        "genomics_fastq_foundation_chain",
        [
            "genomics_01_input_validation",
            "genomics_02_reference_prepare",
            "genomics_03_raw_qc",
            "genomics_04_trim_fastq",
            "genomics_99_provenance_manifest",
        ],
    ),
]


def write_fastq(path: Path) -> None:
    sequence = "ACGT" * 25
    quality = "I" * len(sequence)
    with gzip.open(path, "wt") as handle:
        handle.write(f"@{path.stem}\n{sequence}\n+\n{quality}\n")


def prepare_inputs(work_dir: Path) -> dict[str, Path]:
    work_dir.mkdir(parents=True, exist_ok=True)
    r1 = work_dir / "genomics_demo_R1.fastq.gz"
    r2 = work_dir / "genomics_demo_R2.fastq.gz"
    write_fastq(r1)
    write_fastq(r2)
    samplesheet = work_dir / "genomics_samplesheet.csv"
    samplesheet.write_text(
        "sample_id,fastq_1,fastq_2,single_end\n"
        f"genomics_demo,{r1},{r2},false\n",
        encoding="utf-8",
    )
    fasta = work_dir / "genome.fa"
    fasta.write_text(">chr_demo\n" + ("ACGT" * 100) + "\n", encoding="utf-8")
    return {"samplesheet": samplesheet, "fasta": fasta}


def compile_case(compiler: NextflowCompiler, case: Case, inputs: dict[str, Path], output_root: Path, results_root: Path) -> Path:
    request = RunRequest(
        run_id=case.name,
        omics="genomics",
        execution_profile="local",
        input=str(inputs["samplesheet"]),
        outdir=str(results_root / case.name / "results"),
        selected_steps=case.selected_steps,
        params={
            "genome_fasta": str(inputs["fasta"]),
            "reference_name": "genomics_demo",
            "workflow_name": case.name,
            "threads": 2,
            "quality_cutoff": 20,
            "min_length": 20,
        },
    )
    return compiler.compile(request, output_root / case.name)


def run_selected(cases: list[Case], keep: Path | None) -> list[str]:
    compiler = NextflowCompiler(load_steps(ROOT / "registry" / "steps.yaml"), ROOT / "workflows" / "templates")
    failures = []
    if keep:
        keep.mkdir(parents=True, exist_ok=True)
        inputs = prepare_inputs(keep / "inputs")
        output_root = keep / "compiled"
        results_root = keep / "results"
        output_root.mkdir(parents=True, exist_ok=True)
        results_root.mkdir(parents=True, exist_ok=True)
        return compile_cases(compiler, cases, inputs, output_root, results_root)
    with TemporaryDirectory(prefix="survom_genomics_smoke_") as tmp:
        root = Path(tmp)
        inputs = prepare_inputs(root / "inputs")
        return compile_cases(compiler, cases, inputs, root / "compiled", root / "results")


def compile_cases(compiler: NextflowCompiler, cases: list[Case], inputs: dict[str, Path], output_root: Path, results_root: Path) -> list[str]:
    failures = []
    for case in cases:
        try:
            compiled = compile_case(compiler, case, inputs, output_root, results_root)
            for required in ("main.nf", "params.yaml", "workflow_plan.json", "ui_run_summary.json"):
                if not (compiled / required).exists():
                    raise RuntimeError(f"missing generated file: {required}")
            print(f"[compile ok] {case.name}: {compiled}")
        except Exception as exc:
            failures.append(case.name)
            print(f"[failed] {case.name}: {exc}", file=sys.stderr)
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Compile SurvOm genomics atomic workflow smoke cases.")
    parser.add_argument("--case", action="append", choices=[case.name for case in CASES])
    parser.add_argument("--keep", type=Path)
    args = parser.parse_args()
    selected = [case for case in CASES if not args.case or case.name in set(args.case)]
    failures = run_selected(selected, args.keep)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
