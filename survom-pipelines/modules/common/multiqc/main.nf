process MULTIQC {
    tag "multiqc"
    label 'process_low'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path inputs

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data", emit: data
    path "versions.yml", emit: versions

    script:
    """
    multiqc . --filename multiqc_report.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}
