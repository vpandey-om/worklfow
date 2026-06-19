process BAM_SORT_INDEX {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/09_star/bam", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    tuple val(meta), path("*.sorted.bam.bai"), emit: bai
    path "*.alignment_manifest.tsv", emit: manifest
    path "versions.yml", emit: versions

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam ${bam}
    samtools index ${meta.id}.sorted.bam
    printf "sample_id\\tbam\\tbai\\ttool\\n" > ${meta.id}.alignment_manifest.tsv
    printf "%s\\t%s\\t%s\\tstar\\n" "${meta.id}" "${meta.id}.sorted.bam" "${meta.id}.sorted.bam.bai" >> ${meta.id}.alignment_manifest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | awk '{print \$2}')
    END_VERSIONS
    """
}
