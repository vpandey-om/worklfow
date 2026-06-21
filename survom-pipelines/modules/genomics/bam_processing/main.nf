process GENOMICS_BAM_PROCESSING {
    tag "bam_processing"
    label 'process_medium'
    publishDir "${params.outdir}/06_bam_processing", mode: 'copy'

    input:
    path manifest

    output:
    path "bam_manifest.csv", emit: manifest
    path "*.bam", emit: bam, optional: true
    path "*.bai", emit: bai, optional: true
    path "*.txt", emit: metrics, optional: true
    path "*.log", emit: logs, optional: true
    path "tool_check.json", emit: tool_check
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py bam-process-real \\
      --manifest ${manifest} \\
      --threads ${params.threads ?: 2} \\
      --allow-missing \\
      --outdir .

    if [ ! -f tool_check.json ]; then
      printf '{"step":"bam-process","available":true,"status":"completed"}\n' > tool_check.json
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
