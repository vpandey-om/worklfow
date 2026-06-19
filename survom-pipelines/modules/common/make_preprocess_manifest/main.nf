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
    def single_end = meta.single_end ? 'true' : 'false'
    def strandedness = params.approved_library_type ?: params.inferred_library_type ?: 'unknown'
    """
    printf "${meta.id}\t%s\t%s\t%s\t%s\n" "${fastq1}" "${fastq2}" "${single_end}" "${strandedness}" > ${meta.id}.trimmed_fastq_manifest.tsv

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
    path "trim_manifest.tsv", emit: trim_manifest
    path "versions.yml", emit: versions

    script:
    """
    printf "sample_id\tfastq_1\tfastq_2\tsingle_end\tstrandedness\n" > trimmed_fastq_manifest.tsv
    cat ${fragments} >> trimmed_fastq_manifest.tsv
    cp trimmed_fastq_manifest.tsv trim_manifest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1)
    END_VERSIONS
    """
}
