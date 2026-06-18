process VALIDATE_SAMPLESHEET {
    tag "samplesheet"
    label 'process_low'

    input:
    path samplesheet

    output:
    path "validated_samplesheet.csv", emit: samplesheet
    path "versions.yml", emit: versions

    script:
    """
    cp ${samplesheet} validated_samplesheet.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1)
    END_VERSIONS
    """
}
