#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parent
sys.path.insert(0, str(ROOT))

from builder.compiler import NextflowCompiler
from builder.models import RunRequest
from builder.registry import load_steps


FULL_DEMO_STEPS = [
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
]


def default_samplesheet(output_dir: Path) -> Path:
    fastq_dir = REPO_ROOT / "testdatasets" / "test_data" / "human_chr22_genomics" / "fastq"
    r1 = fastq_dir / "genomics_test_R1.fastq.gz"
    r2 = fastq_dir / "genomics_test_R2.fastq.gz"
    if not r1.exists() or not r2.exists():
        raise FileNotFoundError(f"Demo FASTQs are missing under {fastq_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    samplesheet = output_dir / "genomics_samplesheet.csv"
    samplesheet.write_text(
        "sample_id,fastq_1,fastq_2,single_end\n"
        f"genomics_test,{r1},{r2},false\n",
        encoding="utf-8",
    )
    return samplesheet


def main() -> int:
    parser = argparse.ArgumentParser(description="Compile the full Genomics Step 1-16 demo chain.")
    parser.add_argument("--outdir", type=Path, default=ROOT / "workflows" / "generated" / "genomics_full_demo_chain_cli")
    parser.add_argument("--samplesheet", type=Path)
    parser.add_argument("--genome-fasta", type=Path, default=REPO_ROOT / "testdatasets" / "test_data" / "human_chr22_rnaseq" / "refs" / "chr22_with_ERCC92.fa")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()

    outdir = args.outdir.resolve()
    input_dir = outdir / "inputs"
    samplesheet = args.samplesheet.resolve() if args.samplesheet else default_samplesheet(input_dir)
    genome_fasta = args.genome_fasta.resolve()
    if not genome_fasta.exists():
        raise FileNotFoundError(f"Genome FASTA does not exist: {genome_fasta}")

    request = RunRequest(
        run_id=outdir.name,
        omics="genomics",
        execution_profile="local",
        input=str(samplesheet),
        outdir=str(outdir / "results"),
        selected_steps=FULL_DEMO_STEPS,
        params={
            "genome_fasta": str(genome_fasta),
            "reference_name": "human_chr22_demo",
            "workflow_name": "genomics_full_demo_chain",
            "threads": args.threads,
            "quality_cutoff": 20,
            "min_length": 20,
        },
    )
    compiled = NextflowCompiler(load_steps(ROOT / "registry" / "steps.yaml"), ROOT / "workflows" / "templates").compile(request, outdir)
    print(compiled)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
