process REFERENCE_VALIDATE {
    tag "$selected_route"
    label 'process_low'
    publishDir "${params.outdir}/06A_reference", mode: 'copy'

    input:
    val reference_mode
    val selected_route
    val organism
    val genome_build
    val genome_fasta
    val gtf
    val transcriptome_fasta
    val tx2gene
    val salmon_index
    val star_index
    val hisat2_index

    output:
    path "validated_reference_bundle.json", emit: bundle
    path "reference_manifest.tsv", emit: manifest
    path "reference_validation_report.txt", emit: report
    path "versions.yml", emit: versions

    script:
    """
    python3 - <<'PY'
    import json
    from pathlib import Path

    reference_mode = "${reference_mode}"
    selected_route = "${selected_route}"
    fields = {
        "organism": "${organism}",
        "genome_build": "${genome_build}",
        "genome_fasta": "${genome_fasta}",
        "gtf": "${gtf}",
        "transcriptome_fasta": "${transcriptome_fasta}",
        "tx2gene": "${tx2gene}",
        "salmon_index": "${salmon_index}",
        "star_index": "${star_index}",
        "hisat2_index": "${hisat2_index}",
    }
    missing = []
    warnings = []
    required = ["organism", "genome_build"]
    if selected_route.startswith("salmon") or selected_route in {"quant_only", "salmon_de"}:
        required += ["transcriptome_fasta", "tx2gene", "salmon_index"]
    if selected_route.startswith("star") or selected_route in {"alignment_only", "star_de"}:
        required += ["genome_fasta", "gtf", "star_index"]
    if selected_route.startswith("hisat2"):
        required += ["genome_fasta", "hisat2_index"]
    for key in required:
        value = fields.get(key)
        if not value or value in {"None", "null"}:
            missing.append(f"{key} is required for route {selected_route}.")
        elif key not in {"organism", "genome_build"} and not Path(value).exists():
            missing.append(f"{key} does not exist: {value}")
        elif key == "salmon_index" and not (Path(value) / "versionInfo.json").exists():
            missing.append(f"salmon_index is not a valid Salmon index; missing versionInfo.json: {value}")
        elif key == "star_index" and not (Path(value) / "Genome").exists():
            missing.append(f"star_index is not a valid STAR index; missing Genome file: {value}")
        elif key == "hisat2_index" and not (
            (Path(value) / "genome.1.ht2").exists() or (Path(value) / "genome.1.ht2l").exists()
        ):
            missing.append(f"hisat2_index is not a valid HISAT2 index; missing genome.1.ht2 or genome.1.ht2l: {value}")

    if fields.get("genome_fasta") and fields.get("gtf"):
        warnings.append("Chromosome-name compatibility check is recorded as pending lightweight validation.")
    if fields.get("transcriptome_fasta") and fields.get("tx2gene"):
        warnings.append("Transcript ID version compatibility check is recorded as pending lightweight validation.")

    status = "failed" if missing else "validated"
    bundle = {
        "reference_mode": reference_mode,
        "selected_route": selected_route,
        "status": status,
        "organism": fields["organism"],
        "genome_build": fields["genome_build"],
        "genome_fasta": fields["genome_fasta"],
        "gtf": fields["gtf"],
        "transcriptome_fasta": fields["transcriptome_fasta"],
        "tx2gene": fields["tx2gene"],
        "salmon_index": fields["salmon_index"],
        "star_index": fields["star_index"],
        "hisat2_index": fields["hisat2_index"],
        "missing": missing,
        "warnings": warnings,
    }
    Path("validated_reference_bundle.json").write_text(json.dumps(bundle, indent=2) + "\\n")
    with Path("reference_validation_report.txt").open("w") as handle:
        handle.write(f"Reference validation status: {status}\\n")
        handle.write(f"Mode: {reference_mode}\\nRoute: {selected_route}\\n")
        for msg in missing:
            handle.write(f"ERROR: {msg}\\n")
        for msg in warnings:
            handle.write(f"WARN: {msg}\\n")
    with Path("reference_manifest.tsv").open("w") as handle:
        handle.write("reference_name\\tspecies\\tassembly\\tgenome_fasta\\tannotation_gtf\\ttranscript_fasta\\tsalmon_index\\tstar_index\\thisat2_index\\tcreated_at\\n")
        handle.write(
            f"{fields['organism']}_{fields['genome_build']}\\t"
            f"{fields['organism']}\\t"
            f"{fields['genome_build']}\\t"
            f"{fields['genome_fasta']}\\t"
            f"{fields['gtf']}\\t"
            f"{fields['transcriptome_fasta']}\\t"
            f"{fields['salmon_index']}\\t"
            f"{fields['star_index']}\\t"
            f"{fields['hisat2_index']}\\t"
            f"{__import__('datetime').datetime.utcnow().isoformat()}Z\\n"
        )
    if missing:
        raise SystemExit("\\n".join(missing))
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
