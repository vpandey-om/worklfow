process HISAT2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/09_hisat2", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path hisat2_index

    output:
    tuple val(meta), path("*.hisat2.sam"), emit: sam
    path "*.hisat2.log", emit: log
    path "versions.yml", emit: versions

    script:
    def readArgs = meta.single_end ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    hisat2 \\
      -x ${hisat2_index}/genome \\
      ${readArgs} \\
      -p ${task.cpus} \\
      -S ${meta.id}.hisat2.sam \\
      2> ${meta.id}.hisat2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | head -1 | sed 's/.*version //')
    END_VERSIONS
    """
}
