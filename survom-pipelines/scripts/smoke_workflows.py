#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from builder.compiler import NextflowCompiler
from builder.models import RunRequest
from builder.outputs import collect_run_outputs
from builder.registry import load_steps


@dataclass(frozen=True)
class Case:
    name: str
    selected_steps: list[str]
    input_kind: str
    route: str = "salmon"
    extra_params: dict | None = None


CASES = [
    Case("reference_only", ["rnaseq_06a_reference_build_validation"], "nofastq"),
    Case("raw_qc_only", ["rnaseq_03_raw_read_qc"], "raw"),
    Case("trim_only", ["rnaseq_04_adapter_quality_trimming"], "raw", extra_params={"trimming_tool": "fastp"}),
    Case("post_trim_qc_only", ["rnaseq_05_post_trim_quality_control"], "trimmed"),
    Case("strandedness_only", ["rnaseq_06a_reference_build_validation", "rnaseq_06_strandedness_inference"], "trimmed"),
    Case("salmon_quant_only", ["rnaseq_06a_reference_build_validation", "rnaseq_07_salmon_quantification"], "trimmed"),
    Case("star_align_only", ["rnaseq_06a_reference_build_validation", "rnaseq_09_star_alignment"], "trimmed", route="star"),
    Case("hisat2_align_only", ["rnaseq_06a_reference_build_validation", "rnaseq_09d_hisat2_alignment"], "trimmed", route="hisat2"),
    Case(
        "salmon_count_matrix",
        [
            "rnaseq_06a_reference_build_validation",
            "rnaseq_07_salmon_quantification",
            "rnaseq_08_tximport_gene_summarization",
            "rnaseq_10_count_matrix_qc",
        ],
        "trimmed",
    ),
    Case(
        "star_count_matrix",
        [
            "rnaseq_06a_reference_build_validation",
            "rnaseq_09_star_alignment",
            "rnaseq_09b_bam_sort_index",
            "rnaseq_09c_featurecounts_gene_counting",
            "rnaseq_10_count_matrix_qc",
        ],
        "trimmed",
        route="star",
    ),
    Case(
        "qc_trim_strandedness",
        [
            "rnaseq_03_raw_read_qc",
            "rnaseq_04_adapter_quality_trimming",
            "rnaseq_05_post_trim_quality_control",
            "rnaseq_06a_reference_build_validation",
            "rnaseq_06_strandedness_inference",
        ],
        "raw",
        extra_params={"trimming_tool": "fastp"},
    ),
]


def write_fastq(path: Path, sequence: str) -> None:
    quality = "I" * len(sequence)
    text = f"@{path.stem}\n{sequence}\n+\n{quality}\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def prepare_inputs(work_dir: Path) -> dict[str, Path]:
    work_dir.mkdir(parents=True, exist_ok=True)
    r1 = work_dir / "demo_R1.fastq.gz"
    r2 = work_dir / "demo_R2.fastq.gz"
    write_fastq(r1, "ACGT" * 30)
    write_fastq(r2, "ACGT" * 30)

    raw = work_dir / "raw_manifest.csv"
    raw.write_text(
        "sample_id,fastq_1,fastq_2,single_end,strandedness,platform\n"
        f"demo,{r1},{r2},false,unknown,unknown\n",
        encoding="utf-8",
    )
    trimmed = work_dir / "trim_manifest.tsv"
    trimmed.write_text(
        "sample_id\tfastq_1\tfastq_2\tsingle_end\tstrandedness\n"
        f"demo\t{r1}\t{r2}\tfalse\tunknown\n",
        encoding="utf-8",
    )
    nofastq = work_dir / "no_fastq.csv"
    nofastq.write_text("sample_id,fastq_1,fastq_2,single_end,strandedness,platform\n", encoding="utf-8")
    ref_dir = work_dir / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    reference_sequence = "ACGT" * 80
    (ref_dir / "genome.fa").write_text(f">demo_chr1\n{reference_sequence}\n", encoding="utf-8")
    (ref_dir / "transcriptome.fa").write_text(f">demo_tx1\n{reference_sequence[:240]}\n", encoding="utf-8")
    (ref_dir / "genes.gtf").write_text(
        'demo_chr1\tsmoke\texon\t1\t240\t.\t+\t.\tgene_id "demo_gene1"; transcript_id "demo_tx1";\n',
        encoding="utf-8",
    )
    (ref_dir / "tx2gene.tsv").write_text("tx_id\tgene_id\ndemo_tx1\tdemo_gene1\n", encoding="utf-8")
    return {"raw": raw, "trimmed": trimmed, "nofastq": nofastq, "reference": ref_dir}


def reference_params(route: str, ref_dir: Path) -> dict:
    return {
        "reference_mode": "demo_reference",
        "selected_route": route,
        "organism": "demo",
        "genome_build": "demo_build",
        "genome_fasta": str(ref_dir / "genome.fa"),
        "gtf": str(ref_dir / "genes.gtf"),
        "transcriptome_fasta": str(ref_dir / "transcriptome.fa"),
        "tx2gene": str(ref_dir / "tx2gene.tsv"),
        "salmon_index": None,
        "star_index": None,
        "hisat2_index": None,
        "strandedness_method": "salmon_auto",
        "approved_library_type": "A",
        "quality_threshold": 20,
        "minimum_read_length": 20,
        "trim_poly_g": "auto",
    }


def case_params(case: Case, input_path: Path, ref_dir: Path) -> dict:
    params = reference_params(case.route, ref_dir)
    if case.input_kind == "trimmed":
        params.update({"trim_manifest": str(input_path), "trim_input_mode": "previous_manifest"})
    if case.extra_params:
        params.update(case.extra_params)
    return params


def compile_case(
    compiler: NextflowCompiler,
    case: Case,
    input_path: Path,
    ref_dir: Path,
    output_root: Path,
    results_root: Path,
) -> Path:
    outdir = results_root / case.name / "results"
    request = RunRequest(
        run_id=case.name,
        omics="rnaseq",
        execution_profile="local_docker",
        input=str(input_path),
        outdir=str(outdir),
        selected_steps=case.selected_steps,
        params=case_params(case, input_path, ref_dir),
    )
    return compiler.compile(request, output_root / case.name)


def run_case(compiled: Path, profile: str) -> int:
    command = [
        "nextflow",
        "run",
        "main.nf",
        "-params-file",
        "params.yaml",
        "-profile",
        profile,
        "-work-dir",
        "work",
        "-resume",
    ]
    proc = subprocess.run(command, cwd=compiled, text=True, check=False)
    status = "succeeded" if proc.returncode == 0 else "failed"
    collect_run_outputs(compiled, status=status)
    return proc.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description="Compile or run SurvOm RNA-seq atomic/chained workflow smoke tests.")
    parser.add_argument("--case", action="append", choices=[case.name for case in CASES], help="Run only one case; repeatable.")
    parser.add_argument("--run", action="store_true", help="Execute Nextflow after compiling. Default only compiles.")
    parser.add_argument("--profile", default="local,docker", help="Nextflow profile for --run.")
    parser.add_argument("--keep", type=Path, help="Keep generated files under this directory instead of a temporary directory.")
    args = parser.parse_args()

    selected = [case for case in CASES if not args.case or case.name in set(args.case)]
    compiler = NextflowCompiler(load_steps(ROOT / "registry" / "steps.yaml"), ROOT / "workflows" / "templates")

    if args.keep:
        args.keep.mkdir(parents=True, exist_ok=True)
        inputs = prepare_inputs(args.keep / "inputs")
        output_root = ROOT / "workflows" / "generated" if args.run else args.keep / "compiled"
        results_root = args.keep / "results"
        output_root.mkdir(parents=True, exist_ok=True)
        results_root.mkdir(parents=True, exist_ok=True)
        failures = run_selected(compiler, selected, inputs, output_root, results_root, args.run, args.profile)
        return 1 if failures else 0

    with TemporaryDirectory(prefix="survom_smoke_") as tmp:
        root = Path(tmp)
        inputs = prepare_inputs(root / "inputs")
        output_root = ROOT / "workflows" / "generated" if args.run else root / "compiled"
        results_root = root / "results"
        failures = run_selected(compiler, selected, inputs, output_root, results_root, args.run, args.profile)
        return 1 if failures else 0


def run_selected(
    compiler: NextflowCompiler,
    selected: list[Case],
    inputs: dict[str, Path],
    output_root: Path,
    results_root: Path,
    run_nextflow: bool,
    profile: str,
) -> list[str]:
    failures = []
    for case in selected:
        try:
            compiled = compile_case(compiler, case, inputs[case.input_kind], inputs["reference"], output_root, results_root)
            main_nf = compiled / "main.nf"
            plan = compiled / "workflow_plan.json"
            if not main_nf.exists() or not plan.exists():
                raise RuntimeError("compile completed but main.nf/workflow_plan.json is missing")
            print(f"[compile ok] {case.name}: {compiled}")
            if run_nextflow:
                code = run_case(compiled, profile)
                if code:
                    raise RuntimeError(f"nextflow exited {code}")
                print(f"[run ok] {case.name}")
        except Exception as exc:
            failures.append(case.name)
            print(f"[failed] {case.name}: {exc}", file=sys.stderr)
    if failures:
        print(json.dumps({"failed": failures}, indent=2), file=sys.stderr)
    return failures


if __name__ == "__main__":
    raise SystemExit(main())
