process GENOMICS_INPUT_VALIDATION {
    tag "genomics_inputs"
    label 'process_low'
    publishDir "${params.outdir}/01_input_validation", mode: 'copy'

    input:
    path samplesheet

    output:
    path "validated_genomics_manifest.csv", emit: manifest
    path "validation_report.json", emit: report
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py validate-inputs \\
      --samplesheet ${samplesheet} \\
      --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
