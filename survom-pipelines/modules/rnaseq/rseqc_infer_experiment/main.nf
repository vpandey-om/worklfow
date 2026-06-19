process RSEQC_INFER_EXPERIMENT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/06_strandedness", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.rseqc.strandedness_call.json"), emit: call
    tuple val(meta), path("*.rseqc.strandedness_evidence.tsv"), emit: evidence
    path "*.rseqc_infer_experiment.log", emit: log
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    infer_experiment.py \
      -i ${bam} \
      -r ${params.rseqc_ref_bed} \
      > ${prefix}.rseqc_infer_experiment.log 2>&1

    python3 - <<'PY'
    import json
    import re
    from pathlib import Path

    prefix = "${prefix}"
    log_path = Path(f"{prefix}.rseqc_infer_experiment.log")
    text = log_path.read_text(errors="replace") if log_path.exists() else ""

    frac1 = frac2 = None
    for line in text.splitlines():
        if "1++,1--,2+-,2-+" in line:
            match = re.search(r":\\s*([0-9.]+)", line)
            if match:
                frac1 = float(match.group(1))
        if "1+-,1-+,2++,2--" in line:
            match = re.search(r":\\s*([0-9.]+)", line)
            if match:
                frac2 = float(match.group(1))

    threshold = float("${params.rseqc_stranded_threshold}")
    status = "needs_review"
    inferred = None
    reason = "RSeQC evidence was ambiguous or incomplete."
    if frac1 is not None and frac2 is not None:
        if max(frac1, frac2) < threshold:
            inferred = "unstranded"
            status = "completed"
            reason = "Both stranded fractions were below threshold."
        elif abs(frac1 - frac2) < 0.2:
            inferred = "ambiguous"
            reason = "RSeQC stranded fractions are too close."
        elif frac1 > frac2:
            inferred = "forward"
            status = "completed"
            reason = "RSeQC first-read sense fraction was dominant."
        else:
            inferred = "reverse"
            status = "completed"
            reason = "RSeQC second-read sense fraction was dominant."

    call = {
        "sample_id": prefix,
        "method": "rseqc_infer_experiment",
        "status": status,
        "inferred_library_type": inferred,
        "needs_review": status != "completed",
        "approved": False,
        "approved_library_type": None,
        "reason": reason,
        "fraction_1++,1--,2+-,2-+": frac1,
        "fraction_1+-,1-+,2++,2--": frac2,
        "downstream_param": "inferred_library_type",
    }
    Path(f"{prefix}.rseqc.strandedness_call.json").write_text(json.dumps(call, indent=2) + "\\n")
    with Path(f"{prefix}.rseqc.strandedness_evidence.tsv").open("w") as handle:
        handle.write("key\\tvalue\\n")
        for key, value in call.items():
            handle.write(f"{key}\\t{value}\\n")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version 2>&1 | tail -1)
    END_VERSIONS
    """
}
