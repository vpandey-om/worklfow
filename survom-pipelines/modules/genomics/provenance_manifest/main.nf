process GENOMICS_PROVENANCE_MANIFEST {
    tag "provenance"
    label 'process_low'
    publishDir "${params.outdir}/99_provenance", mode: 'copy'

    output:
    path "run_manifest.json", emit: manifest
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py provenance \\
      --workflow-name ${params.workflow_name ?: 'genomics_workflow'} \\
      --run-id ${params.run_id} \\
      --parameters-json '{}' \\
      --inputs-json '[]' \\
      --outputs-json '[]' \\
      --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
