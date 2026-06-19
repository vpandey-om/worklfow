process BUILD_SALMON_INDEX {
    tag "salmon_index"
    label 'process_high'
    publishDir "${params.outdir}/06A_reference", mode: 'copy'

    input:
    path transcriptome_fasta

    output:
    path "salmon_index", emit: index
    path "versions.yml", emit: versions

    script:
    """
    salmon index -t ${transcriptome_fasta} -i salmon_index -p ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}
