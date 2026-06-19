process COUNT_MATRIX_QC {
    tag "count_matrix_qc"
    label 'process_low'
    publishDir "${params.outdir}/10_count_matrix_qc", mode: 'copy'

    input:
    path matrix

    output:
    path "count_matrix_qc.json", emit: qc

    script:
    """
    python3 - <<'PY'
    import json
    from pathlib import Path
    matrix = Path("${matrix}")
    rows = matrix.read_text(errors="replace").splitlines() if matrix.exists() else []
    Path("count_matrix_qc.json").write_text(json.dumps({
        "matrix": str(matrix),
        "exists": matrix.exists(),
        "line_count": len(rows),
        "status": "completed" if matrix.exists() else "failed"
    }, indent=2) + "\\n")
    PY
    """
}
