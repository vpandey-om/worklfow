process GENOMICS_STRUCTURAL_VARIANTS {
    tag "structural_variants"
    label 'process_high'
    publishDir "${params.outdir}/13_structural_variants", mode: 'copy'

    output:
    path "step_status.json", emit: status
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py sv-status --outdir .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
