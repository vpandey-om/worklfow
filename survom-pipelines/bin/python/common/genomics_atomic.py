#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
DEFAULT_GENOMICS_ENV_BIN = Path("/home/vikash/miniconda3/envs/survom-genomics/bin")
GENOMICS_TOOL_GROUPS = {
    "bwa-align": ["bwa-mem2", "bwa"],
    "bam-process": ["samtools"],
    "mark-duplicates": ["samtools", "picard"],
    "coverage-qc": ["mosdepth", "qualimap"],
    "bqsr": ["gatk"],
    "germline-variants": ["gatk", "run_deepvariant"],
    "somatic-variants": ["gatk"],
    "variant-filter": ["bcftools", "gatk"],
    "sv-calling": ["manta", "delly", "sniffles"],
    "cnv-calling": ["cnvkit.py"],
    "variant-annotation": ["vep", "snpEff"],
    "multiqc": ["multiqc"],
}


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def find_tool(tool: str) -> str:
    found = shutil.which(tool)
    if found:
        return found
    env_bin = Path(os.environ.get("SURVOM_GENOMICS_BIN", str(DEFAULT_GENOMICS_ENV_BIN))).expanduser()
    candidate = env_bin / tool
    return str(candidate) if candidate.exists() and candidate.is_file() else ""


def read_table(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader), list(reader.fieldnames or [])


def write_table(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_command(command: list[str], log_path: Path) -> None:
    mkdir(log_path.parent)
    with log_path.open("w") as log:
        log.write("$ " + " ".join(command) + "\n\n")
        log.flush()
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)


def materialize_reference(reference: Path, outdir: Path) -> Path:
    target = outdir / "reference.fa"
    if target.exists() and target.stat().st_size > 0:
        return target
    if reference.name.lower().endswith(".gz"):
        with gzip.open(reference, "rb") as source, target.open("wb") as destination:
            shutil.copyfileobj(source, destination)
    else:
        shutil.copy2(reference, target)
    return target


def resolve_path(value: str, base: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def is_fastq(path: Path) -> bool:
    lower = path.name.lower()
    return lower.endswith(FASTQ_SUFFIXES)


def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return path.open("r", errors="replace")


def count_fastq_reads(path: Path, max_records: int | None = None) -> tuple[int, int]:
    reads = 0
    bases = 0
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not (seq and plus and qual and header.startswith("@") and plus.startswith("+")):
                raise ValueError(f"Invalid FASTQ record near read {reads + 1}: {path}")
            reads += 1
            bases += len(seq.strip())
            if max_records and reads >= max_records:
                break
    return reads, bases


def validate_inputs(args: argparse.Namespace) -> int:
    samplesheet = Path(args.samplesheet).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    rows, columns = read_table(samplesheet)
    required = {"sample_id"}
    errors: list[str] = []
    warnings: list[str] = []
    if not rows:
        errors.append("Sample sheet has no rows.")
    if "fastq_1" in columns:
        required.add("fastq_1")
    elif "bam" in columns:
        required.add("bam")
    elif "vcf" in columns:
        required.add("vcf")
    else:
        errors.append("Sample sheet must contain one of fastq_1, bam, or vcf.")
    missing = sorted(required - set(columns))
    if missing:
        errors.append(f"Sample sheet missing required column(s): {', '.join(missing)}")

    validated = []
    seen = set()
    base = samplesheet.parent
    for index, row in enumerate(rows, start=2):
        sample_id = (row.get("sample_id") or "").strip()
        if not sample_id:
            errors.append(f"Row {index}: sample_id is missing.")
            continue
        if sample_id in seen:
            errors.append(f"Duplicate sample_id: {sample_id}")
        seen.add(sample_id)
        record = {"sample_id": sample_id, "input_type": "", "path_1": "", "path_2": "", "single_end": "", "status": "valid"}
        if row.get("fastq_1"):
            fq1 = resolve_path(row["fastq_1"].strip(), base)
            fq2_raw = (row.get("fastq_2") or "").strip()
            single_end = str(row.get("single_end") or "").lower() in {"true", "1", "yes", "single"}
            if not fq1.exists() or not fq1.is_file():
                errors.append(f"{sample_id}: fastq_1 does not exist: {fq1}")
            elif not is_fastq(fq1):
                errors.append(f"{sample_id}: fastq_1 is not FASTQ: {fq1.name}")
            if fq2_raw:
                fq2 = resolve_path(fq2_raw, base)
                if not fq2.exists() or not fq2.is_file():
                    errors.append(f"{sample_id}: fastq_2 does not exist: {fq2}")
                elif not is_fastq(fq2):
                    errors.append(f"{sample_id}: fastq_2 is not FASTQ: {fq2.name}")
                single_end = False
                record["path_2"] = str(fq2)
            elif not single_end:
                errors.append(f"{sample_id}: paired-end row is missing fastq_2 or single_end=true.")
            record.update({"input_type": "fastq", "path_1": str(fq1), "single_end": "true" if single_end else "false"})
        elif row.get("bam"):
            bam = resolve_path(row["bam"].strip(), base)
            if not bam.exists() or not bam.is_file():
                errors.append(f"{sample_id}: BAM does not exist: {bam}")
            record.update({"input_type": "bam", "path_1": str(bam), "single_end": ""})
        elif row.get("vcf"):
            vcf = resolve_path(row["vcf"].strip(), base)
            if not vcf.exists() or not vcf.is_file():
                errors.append(f"{sample_id}: VCF does not exist: {vcf}")
            record.update({"input_type": "vcf", "path_1": str(vcf), "single_end": ""})
        validated.append(record)

    if errors:
        for record in validated:
            record["status"] = "invalid"
    write_table(outdir / "validated_genomics_manifest.csv", validated, ["sample_id", "input_type", "path_1", "path_2", "single_end", "status"])
    report = {
        "schema_version": "1.0",
        "created_at": now_iso(),
        "samplesheet": str(samplesheet),
        "valid": not errors,
        "errors": errors,
        "warnings": warnings,
        "sample_count": len(rows),
    }
    (outdir / "validation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return 0 if not errors else 2


def fasta_lengths(fasta: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current = None
    with open_text(fasta) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current:
                lengths[current] += len(line)
    if not lengths:
        raise ValueError(f"No FASTA records found: {fasta}")
    return lengths


def reference_prepare(args: argparse.Namespace) -> int:
    genome_fasta = Path(args.genome_fasta).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    errors = []
    if not genome_fasta.exists():
        errors.append(f"Genome FASTA does not exist: {genome_fasta}")
    if errors:
        (outdir / "reference_prepare_report.json").write_text(json.dumps({"valid": False, "errors": errors}, indent=2), encoding="utf-8")
        return 2
    lengths = fasta_lengths(genome_fasta)
    fai = outdir / f"{genome_fasta.name}.fai"
    with fai.open("w") as handle:
        for contig, length in lengths.items():
            handle.write(f"{contig}\t{length}\t0\t0\t0\n")
    dict_path = outdir / f"{genome_fasta.stem}.dict"
    with dict_path.open("w") as handle:
        handle.write("@HD\tVN:1.6\n")
        for contig, length in lengths.items():
            handle.write(f"@SQ\tSN:{contig}\tLN:{length}\n")
    bwa_available = find_tool("bwa-mem2") or find_tool("bwa")
    bundle = {
        "schema_version": "1.0",
        "created_at": now_iso(),
        "reference_name": args.reference_name,
        "genome_fasta": str(genome_fasta),
        "fai": str(fai),
        "sequence_dictionary": str(dict_path),
        "bwa_index_prefix": str(outdir / "bwa" / genome_fasta.name) if bwa_available else "",
        "bwa_index_status": "not_built_tool_missing" if not bwa_available else "tool_available_not_run_by_lightweight_test",
        "contigs": lengths,
        "warnings": [] if bwa_available else ["bwa/bwa-mem2 is not on PATH; build the BWA index in a tool-enabled environment."],
    }
    (outdir / "reference_bundle.json").write_text(json.dumps(bundle, indent=2), encoding="utf-8")
    return 0


def raw_qc(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    rows, _ = read_table(manifest)
    metrics = []
    for row in rows:
        if row.get("input_type") and row.get("input_type") != "fastq":
            continue
        sample_id = row.get("sample_id") or "sample"
        for read_label, column in (("R1", "path_1"), ("R2", "path_2")):
            value = (row.get(column) or "").strip()
            if not value:
                continue
            path = Path(value)
            reads, bases = count_fastq_reads(path)
            metrics.append({"sample_id": sample_id, "read": read_label, "fastq": str(path), "reads": reads, "bases": bases})
    write_table(outdir / "raw_qc_metrics.csv", metrics, ["sample_id", "read", "fastq", "reads", "bases"])
    (outdir / "raw_qc_summary.json").write_text(json.dumps({"created_at": now_iso(), "records": metrics}, indent=2), encoding="utf-8")
    return 0


def trim_fastq(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    trimmed_dir = mkdir(outdir / "trimmed_fastq")
    rows, _ = read_table(manifest)
    out_rows = []
    for row in rows:
        if row.get("input_type") and row.get("input_type") != "fastq":
            continue
        sample_id = row.get("sample_id") or "sample"
        copied = []
        for column in ("path_1", "path_2"):
            value = (row.get(column) or "").strip()
            if not value:
                copied.append("")
                continue
            source = Path(value)
            suffix = ".fastq.gz" if source.name.lower().endswith(".gz") else ".fastq"
            target = trimmed_dir / f"{sample_id}_{'R1' if column == 'path_1' else 'R2'}.trimmed{suffix}"
            shutil.copy2(source, target)
            copied.append(str(target.resolve()))
        out_rows.append({
            "sample_id": sample_id,
            "fastq_1": copied[0],
            "fastq_2": copied[1],
            "single_end": row.get("single_end") or ("true" if not copied[1] else "false"),
            "strandedness": "unknown",
        })
    write_table(outdir / "trimmed_fastq_manifest.csv", out_rows, ["sample_id", "fastq_1", "fastq_2", "single_end", "strandedness"])
    (outdir / "trim_report.json").write_text(json.dumps({
        "created_at": now_iso(),
        "mode": "copy_for_lightweight_atomic_test",
        "quality_cutoff": args.quality_cutoff,
        "min_length": args.min_length,
    }, indent=2), encoding="utf-8")
    return 0


def bwa_align(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).resolve()
    reference = Path(args.reference).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    rows, _ = read_table(manifest)
    bwa = find_tool("bwa-mem2") or find_tool("bwa")
    if not bwa:
        return require_any_tool(GENOMICS_TOOL_GROUPS["bwa-align"], outdir, "bwa-align", args.allow_missing)
    work_ref = materialize_reference(reference, outdir)
    if not (outdir / "reference.fa.bwt.2bit.64").exists() and not (outdir / "reference.fa.bwt").exists():
        run_command([bwa, "index", str(work_ref)], outdir / "bwa_index.log")
    out_rows = []
    for row in rows:
        sample_id = row.get("sample_id") or "sample"
        r1 = Path(row.get("fastq_1") or row.get("path_1") or "").resolve()
        r2_value = row.get("fastq_2") or row.get("path_2") or ""
        single_end = str(row.get("single_end") or "").lower() in {"true", "1", "yes", "single"}
        sam = outdir / f"{sample_id}.sam"
        read_group = f"@RG\tID:{sample_id}\tSM:{sample_id}\tPL:ILLUMINA"
        command = [bwa, "mem", "-t", str(args.threads), "-R", read_group, str(work_ref), str(r1)]
        if r2_value and not single_end:
            command.append(str(Path(r2_value).resolve()))
        with sam.open("w") as sam_handle, (outdir / f"{sample_id}.bwa.log").open("w") as log_handle:
            log_handle.write("$ " + " ".join(command) + "\n\n")
            subprocess.run(command, stdout=sam_handle, stderr=log_handle, check=True)
        out_rows.append({"sample_id": sample_id, "sam": str(sam), "tool": Path(bwa).name})
    write_table(outdir / "alignment_manifest.csv", out_rows, ["sample_id", "sam", "tool"])
    return 0


def bam_process(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    rows, _ = read_table(manifest)
    samtools = find_tool("samtools")
    if not samtools:
        return require_any_tool(GENOMICS_TOOL_GROUPS["bam-process"], outdir, "bam-process", args.allow_missing)
    out_rows = []
    for row in rows:
        sample_id = row.get("sample_id") or "sample"
        sam = Path(row.get("sam") or row.get("bam") or "").resolve()
        bam = outdir / f"{sample_id}.sorted.bam"
        bai = outdir / f"{sample_id}.sorted.bam.bai"
        flagstat = outdir / f"{sample_id}.flagstat.txt"
        stats = outdir / f"{sample_id}.stats.txt"
        run_command([samtools, "sort", "-@", str(args.threads), "-o", str(bam), str(sam)], outdir / f"{sample_id}.sort.log")
        run_command([samtools, "index", str(bam)], outdir / f"{sample_id}.index.log")
        with flagstat.open("w") as handle:
            subprocess.run([samtools, "flagstat", str(bam)], stdout=handle, stderr=subprocess.PIPE, check=True, text=True)
        with stats.open("w") as handle:
            subprocess.run([samtools, "stats", str(bam)], stdout=handle, stderr=subprocess.PIPE, check=True, text=True)
        out_rows.append({"sample_id": sample_id, "bam": str(bam), "bai": str(bai), "flagstat": str(flagstat), "stats": str(stats), "tool": "samtools"})
    write_table(outdir / "bam_manifest.csv", out_rows, ["sample_id", "bam", "bai", "flagstat", "stats", "tool"])
    return 0


def mark_duplicates(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    rows, _ = read_table(manifest)
    picard = find_tool("picard")
    samtools = find_tool("samtools")
    out_rows = []
    for row in rows:
        sample_id = row.get("sample_id") or "sample"
        bam = Path(row.get("bam") or "").resolve()
        marked = outdir / f"{sample_id}.marked.bam"
        metrics = outdir / f"{sample_id}.markdup_metrics.txt"
        if picard:
            run_command(
                [
                    picard,
                    "MarkDuplicates",
                    f"I={bam}",
                    f"O={marked}",
                    f"M={metrics}",
                    "CREATE_INDEX=true",
                    "VALIDATION_STRINGENCY=LENIENT",
                ],
                outdir / f"{sample_id}.picard_markdup.log",
            )
        else:
            shutil.copy2(bam, marked)
            metrics.write_text("Picard not available; copied BAM without duplicate marking.\n", encoding="utf-8")
            bai_source = Path(row.get("bai") or "")
            if bai_source.exists():
                shutil.copy2(bai_source, Path(str(marked) + ".bai"))
        bai = Path(str(marked) + ".bai")
        picard_bai = marked.with_suffix(".bai")
        if picard_bai.exists() and not bai.exists():
            shutil.copy2(picard_bai, bai)
        if samtools and not bai.exists():
            run_command([samtools, "index", str(marked)], outdir / f"{sample_id}.markdup_index.log")
        out_rows.append({"sample_id": sample_id, "bam": str(marked), "bai": str(bai), "metrics": str(metrics), "tool": "picard" if picard else "copy"})
    write_table(outdir / "marked_bam_manifest.csv", out_rows, ["sample_id", "bam", "bai", "metrics", "tool"])
    return 0


def coverage_qc(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    rows, _ = read_table(manifest)
    mosdepth = find_tool("mosdepth")
    if not mosdepth:
        return require_any_tool(GENOMICS_TOOL_GROUPS["coverage-qc"], outdir, "coverage-qc", args.allow_missing)
    out_rows = []
    for row in rows:
        sample_id = row.get("sample_id") or "sample"
        bam = Path(row.get("bam") or "").resolve()
        prefix = outdir / sample_id
        run_command([mosdepth, "-t", str(args.threads), str(prefix), str(bam)], outdir / f"{sample_id}.mosdepth.log")
        out_rows.append({"sample_id": sample_id, "summary": str(prefix) + ".mosdepth.summary.txt", "per_base": str(prefix) + ".per-base.bed.gz", "tool": "mosdepth"})
    write_table(outdir / "coverage_manifest.csv", out_rows, ["sample_id", "summary", "per_base", "tool"])
    return 0


def write_skipped_status(args: argparse.Namespace, step_name: str, reason: str) -> int:
    outdir = mkdir(Path(args.outdir).resolve())
    (outdir / "step_status.json").write_text(json.dumps({"step": step_name, "status": "skipped", "reason": reason}, indent=2), encoding="utf-8")
    return 0


def germline_variants(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).resolve()
    reference = Path(args.reference).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    rows, _ = read_table(manifest)
    bcftools = find_tool("bcftools")
    if not bcftools:
        return require_any_tool(GENOMICS_TOOL_GROUPS["germline-variants"], outdir, "germline-variants", args.allow_missing)
    work_ref = materialize_reference(reference, outdir)
    samtools = find_tool("samtools")
    if samtools and not Path(str(work_ref) + ".fai").exists():
        run_command([samtools, "faidx", str(work_ref)], outdir / "reference_faidx.log")
    out_rows = []
    for row in rows:
        sample_id = row.get("sample_id") or "sample"
        bam = Path(row.get("bam") or "").resolve()
        raw_vcf = outdir / f"{sample_id}.raw.vcf.gz"
        mpileup = subprocess.Popen([bcftools, "mpileup", "-Ou", "-f", str(work_ref), str(bam)], stdout=subprocess.PIPE, stderr=(outdir / f"{sample_id}.mpileup.log").open("w"))
        call_log = (outdir / f"{sample_id}.bcftools_call.log").open("w")
        call = subprocess.run([bcftools, "call", "-mv", "-Oz", "-o", str(raw_vcf)], stdin=mpileup.stdout, stderr=call_log)
        if mpileup.stdout:
            mpileup.stdout.close()
        mpileup.wait()
        call_log.close()
        if call.returncode != 0:
            raise subprocess.CalledProcessError(call.returncode, "bcftools call")
        run_command([bcftools, "index", "-f", str(raw_vcf)], outdir / f"{sample_id}.bcftools_index.log")
        out_rows.append({"sample_id": sample_id, "vcf": str(raw_vcf), "tbi": str(raw_vcf) + ".csi", "tool": "bcftools"})
    write_table(outdir / "raw_vcf_manifest.csv", out_rows, ["sample_id", "vcf", "tbi", "tool"])
    return 0


def variant_filter(args: argparse.Namespace) -> int:
    manifest = Path(args.manifest).resolve()
    outdir = mkdir(Path(args.outdir).resolve())
    rows, _ = read_table(manifest)
    bcftools = find_tool("bcftools")
    if not bcftools:
        return require_any_tool(GENOMICS_TOOL_GROUPS["variant-filter"], outdir, "variant-filter", args.allow_missing)
    out_rows = []
    for row in rows:
        sample_id = row.get("sample_id") or "sample"
        vcf = Path(row.get("vcf") or "").resolve()
        filtered = outdir / f"{sample_id}.filtered.vcf.gz"
        run_command([bcftools, "filter", "-Oz", "-o", str(filtered), "-i", args.expression, str(vcf)], outdir / f"{sample_id}.bcftools_filter.log")
        run_command([bcftools, "index", "-f", str(filtered)], outdir / f"{sample_id}.bcftools_index.log")
        out_rows.append({"sample_id": sample_id, "vcf": str(filtered), "tbi": str(filtered) + ".csi", "tool": "bcftools", "expression": args.expression})
    write_table(outdir / "filtered_vcf_manifest.csv", out_rows, ["sample_id", "vcf", "tbi", "tool", "expression"])
    return 0


def multiqc_report(args: argparse.Namespace) -> int:
    outdir = mkdir(Path(args.outdir).resolve())
    multiqc = find_tool("multiqc")
    if not multiqc:
        return require_any_tool(GENOMICS_TOOL_GROUPS["multiqc"], outdir, "multiqc", args.allow_missing)
    run_command([multiqc, str(Path(args.search_dir).resolve()), "-o", str(outdir)], outdir / "multiqc.log")
    return 0


def require_tool(tool: str, outdir: Path) -> int:
    mkdir(outdir)
    if find_tool(tool):
        (outdir / "tool_check.json").write_text(json.dumps({"tool": tool, "available": True}, indent=2), encoding="utf-8")
        return 0
    message = f"Required tool is not available on PATH: {tool}. Install it or run this atomic step with the configured container/profile."
    (outdir / "tool_check.json").write_text(json.dumps({"tool": tool, "available": False, "error": message}, indent=2), encoding="utf-8")
    print(message, file=sys.stderr)
    return 127


def require_any_tool(tools: list[str], outdir: Path, step_name: str, allow_missing: bool = False) -> int:
    mkdir(outdir)
    available = [tool for tool in tools if find_tool(tool)]
    if available:
        (outdir / "tool_check.json").write_text(
            json.dumps(
                {"step": step_name, "tools": tools, "available": True, "selected_tool": available[0], "status": "ready"},
                indent=2,
            ),
            encoding="utf-8",
        )
        return 0
    message = (
        f"{step_name} requires one of these tools on PATH: {', '.join(tools)}. "
        "Install the tool or run this atomic step with the configured container/profile."
    )
    (outdir / "tool_check.json").write_text(
        json.dumps(
            {
                "step": step_name,
                "tools": tools,
                "available": False,
                "status": "skipped_missing_tool" if allow_missing else "failed_missing_tool",
                "error": message,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(message, file=sys.stderr)
    return 0 if allow_missing else 127


def tool_inventory(args: argparse.Namespace) -> int:
    outdir = mkdir(Path(args.outdir).resolve())
    rows = []
    for step_name, tools in GENOMICS_TOOL_GROUPS.items():
        available = [tool for tool in tools if find_tool(tool)]
        rows.append({
            "step": step_name,
            "required_tools": "|".join(tools),
            "available_tools": "|".join(available),
            "selected_tool": available[0] if available else "",
            "status": "ready" if available else "missing_tool",
        })
    write_table(
        outdir / "tool_inventory.csv",
        rows,
        ["step", "required_tools", "available_tools", "selected_tool", "status"],
    )
    (outdir / "tool_inventory.json").write_text(
        json.dumps({"created_at": now_iso(), "tools": rows}, indent=2),
        encoding="utf-8",
    )
    return 0


def provenance(args: argparse.Namespace) -> int:
    outdir = mkdir(Path(args.outdir).resolve())
    manifest = {
        "schema_version": "1.0",
        "created_at": now_iso(),
        "workflow_name": args.workflow_name,
        "run_id": args.run_id,
        "parameters": json.loads(args.parameters_json or "{}"),
        "inputs": json.loads(args.inputs_json or "[]"),
        "outputs": json.loads(args.outputs_json or "[]"),
    }
    (outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="SurvOm lightweight genomics atomic helpers.")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("validate-inputs")
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=validate_inputs)
    p = sub.add_parser("reference-prepare")
    p.add_argument("--genome-fasta", required=True)
    p.add_argument("--reference-name", default="custom_reference")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=reference_prepare)
    p = sub.add_parser("raw-qc")
    p.add_argument("--manifest", required=True)
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=raw_qc)
    p = sub.add_parser("trim-fastq")
    p.add_argument("--manifest", required=True)
    p.add_argument("--quality-cutoff", default="20")
    p.add_argument("--min-length", default="20")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=trim_fastq)
    real = sub.add_parser("bwa-align-real")
    real.add_argument("--manifest", required=True)
    real.add_argument("--reference", required=True)
    real.add_argument("--threads", default="2")
    real.add_argument("--outdir", required=True)
    real.add_argument("--allow-missing", action="store_true")
    real.set_defaults(func=bwa_align)
    real = sub.add_parser("bam-process-real")
    real.add_argument("--manifest", required=True)
    real.add_argument("--threads", default="2")
    real.add_argument("--outdir", required=True)
    real.add_argument("--allow-missing", action="store_true")
    real.set_defaults(func=bam_process)
    real = sub.add_parser("mark-duplicates-real")
    real.add_argument("--manifest", required=True)
    real.add_argument("--outdir", required=True)
    real.set_defaults(func=mark_duplicates)
    real = sub.add_parser("coverage-qc-real")
    real.add_argument("--manifest", required=True)
    real.add_argument("--threads", default="2")
    real.add_argument("--outdir", required=True)
    real.add_argument("--allow-missing", action="store_true")
    real.set_defaults(func=coverage_qc)
    real = sub.add_parser("bqsr-status")
    real.add_argument("--outdir", required=True)
    real.set_defaults(func=lambda args: write_skipped_status(args, "bqsr", "BQSR requires known-sites VCFs for the selected organism/reference."))
    real = sub.add_parser("germline-variants-real")
    real.add_argument("--manifest", required=True)
    real.add_argument("--reference", required=True)
    real.add_argument("--outdir", required=True)
    real.add_argument("--allow-missing", action="store_true")
    real.set_defaults(func=germline_variants)
    real = sub.add_parser("somatic-status")
    real.add_argument("--outdir", required=True)
    real.set_defaults(func=lambda args: write_skipped_status(args, "somatic-variants", "Somatic calling requires tumor/normal metadata plus optional resources."))
    real = sub.add_parser("variant-filter-real")
    real.add_argument("--manifest", required=True)
    real.add_argument("--expression", default="QUAL>=20")
    real.add_argument("--outdir", required=True)
    real.add_argument("--allow-missing", action="store_true")
    real.set_defaults(func=variant_filter)
    real = sub.add_parser("sv-status")
    real.add_argument("--outdir", required=True)
    real.set_defaults(func=lambda args: write_skipped_status(args, "structural-variants", "SV callers are installed, but the tiny demo FASTQs are not suitable for meaningful SV calling."))
    real = sub.add_parser("cnv-status")
    real.add_argument("--outdir", required=True)
    real.set_defaults(func=lambda args: write_skipped_status(args, "cnv-calling", "CNVkit requires target/access/reference configuration for the assay."))
    real = sub.add_parser("annotation-status")
    real.add_argument("--outdir", required=True)
    real.set_defaults(func=lambda args: write_skipped_status(args, "variant-annotation", "Annotation requires a SnpEff/VEP database matching the selected reference."))
    real = sub.add_parser("multiqc-real")
    real.add_argument("--search-dir", required=True)
    real.add_argument("--outdir", required=True)
    real.add_argument("--allow-missing", action="store_true")
    real.set_defaults(func=multiqc_report)
    for name, tools in GENOMICS_TOOL_GROUPS.items():
        p = sub.add_parser(name)
        p.add_argument("--outdir", required=True)
        p.add_argument("--allow-missing", action="store_true", help="Write a skipped status instead of failing when tools are absent.")
        p.set_defaults(func=lambda args, ts=tools, n=name: require_any_tool(ts, Path(args.outdir).resolve(), n, args.allow_missing))
    p = sub.add_parser("tool-inventory")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=tool_inventory)
    p = sub.add_parser("provenance")
    p.add_argument("--workflow-name", required=True)
    p.add_argument("--run-id", required=True)
    p.add_argument("--parameters-json", default="{}")
    p.add_argument("--inputs-json", default="[]")
    p.add_argument("--outputs-json", default="[]")
    p.add_argument("--outdir", required=True)
    p.set_defaults(func=provenance)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
