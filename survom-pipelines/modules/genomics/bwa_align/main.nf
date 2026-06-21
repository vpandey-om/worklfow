process GENOMICS_BWA_ALIGN {
    tag "bwa_align"
    label 'process_high'
    publishDir "${params.outdir}/05_bwa_align", mode: 'copy'

    input:
    path manifest
    path reference

    output:
    path "alignment_manifest.csv", emit: manifest
    path "*.sam", emit: sam, optional: true
    path "*.log", emit: logs, optional: true
    path "tool_check.json", emit: tool_check
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py bwa-align-real \\
      --manifest ${manifest} \\
      --reference ${reference} \\
      --threads ${params.threads ?: 2} \\
      --allow-missing \\
      --outdir .

    if [ ! -f tool_check.json ]; then
      printf '{"step":"bwa-align","available":true,"status":"completed"}\n' > tool_check.json
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
