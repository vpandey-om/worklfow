process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    """
    echo "STAR placeholder module; provide genome index before production execution." > ${meta.id}.bam
    """
}
