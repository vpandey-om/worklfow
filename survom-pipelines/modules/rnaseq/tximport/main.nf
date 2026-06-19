process TXIMPORT {
    tag "tximport"
    label 'process_low'
    publishDir "${params.outdir}/08_tximport", mode: 'copy'

    input:
    path quant_dirs

    output:
    path "gene_count_matrix.tsv", emit: counts
    path "tximport_summary.json", emit: summary

    script:
    """
    python3 - <<'PY'
    import json
    from pathlib import Path
    quant_dirs = "${quant_dirs}".split()
    Path("gene_count_matrix.tsv").write_text("gene_id\\tcount\\nplaceholder_gene\\t0\\n")
    Path("tximport_summary.json").write_text(json.dumps({
        "status": "completed",
        "tx2gene": "${params.tx2gene}",
        "quant_dirs": quant_dirs,
        "note": "Placeholder tximport summary; replace with R tximport for production."
    }, indent=2) + "\\n")
    PY
    """
}
