process MAKE_PREPROCESS_MANIFEST {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    path "*.trimmed_fastq_manifest.tsv", emit: fragments
    path "versions.yml", emit: versions

    script:
    def read_list = reads instanceof List ? reads : [reads]
    def fastq1 = read_list[0]?.name ?: ''
    def fastq2 = read_list.size() > 1 ? read_list[1].name : ''
    """
    printf "${meta.id}\t%s\t%s\t%s\n" "${meta.single_end ? 'single' : 'paired'}" "${fastq1}" "${fastq2}" > ${meta.id}.trimmed_fastq_manifest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1)
    END_VERSIONS
    """
}

process MERGE_PREPROCESS_MANIFEST {
    tag "trimmed_fastq_manifest"
    label 'process_low'
    publishDir "${params.outdir}/04_trimmed", mode: 'copy'

    input:
    path fragments

    output:
    path "trimmed_fastq_manifest.tsv", emit: manifest
    path "versions.yml", emit: versions

    script:
    """
    printf "sample_id\tlibrary_layout\tfastq_1\tfastq_2\n" > trimmed_fastq_manifest.tsv
    cat ${fragments} >> trimmed_fastq_manifest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1)
    END_VERSIONS
    """
}
