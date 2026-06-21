process GENOMICS_SOMATIC_VARIANTS {
    tag "somatic_variants"
    label 'process_high'
    publishDir "${params.outdir}/11_somatic_variants", mode: 'copy'

    output:
    path "step_status.json", emit: status
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py somatic-status --outdir .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
