process GENOMICS_MARK_DUPLICATES {
    tag "mark_duplicates"
    label 'process_medium'
    publishDir "${params.outdir}/07_mark_duplicates", mode: 'copy'

    input:
    path manifest

    output:
    path "marked_bam_manifest.csv", emit: manifest
    path "*.bam", emit: bam, optional: true
    path "*.bai", emit: bai, optional: true
    path "*metrics.txt", emit: metrics, optional: true
    path "*.log", emit: logs, optional: true
    path "tool_check.json", emit: tool_check
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py mark-duplicates-real \\
      --manifest ${manifest} \\
      --outdir .

    printf '{"step":"mark-duplicates","available":true,"status":"completed"}\n' > tool_check.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
