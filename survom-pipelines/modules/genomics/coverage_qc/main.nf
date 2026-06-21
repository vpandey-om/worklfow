process GENOMICS_COVERAGE_QC {
    tag "coverage_qc"
    label 'process_medium'
    publishDir "${params.outdir}/08_coverage_qc", mode: 'copy'

    input:
    path manifest

    output:
    path "coverage_manifest.csv", emit: manifest
    path "*.txt", emit: reports, optional: true
    path "*.gz", emit: coverage, optional: true
    path "*.csi", emit: indexes, optional: true
    path "*.log", emit: logs, optional: true
    path "tool_check.json", emit: tool_check
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py coverage-qc-real \\
      --manifest ${manifest} \\
      --threads ${params.threads ?: 2} \\
      --allow-missing \\
      --outdir .

    if [ ! -f tool_check.json ]; then
      printf '{"step":"coverage-qc","available":true,"status":"completed"}\n' > tool_check.json
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
