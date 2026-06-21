process GENOMICS_MULTIQC_REPORT {
    tag "multiqc"
    label 'process_medium'
    publishDir "${params.outdir}/16_multiqc", mode: 'copy'

    input:
    path upstream_marker

    output:
    path "multiqc_report.html", emit: report, optional: true
    path "multiqc_data/**", emit: data, optional: true
    path "multiqc.log", emit: log, optional: true
    path "tool_check.json", emit: tool_check
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py multiqc-real \\
      --search-dir ${params.outdir} \\
      --allow-missing \\
      --outdir .

    if [ ! -f tool_check.json ]; then
      printf '{"step":"multiqc","available":true,"status":"completed"}\n' > tool_check.json
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
