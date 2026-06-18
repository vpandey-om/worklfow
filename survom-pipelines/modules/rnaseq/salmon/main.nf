process SALMON {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}"), emit: quant

    script:
    """
    mkdir -p ${meta.id}
    echo "SALMON placeholder module; provide transcriptome index before production execution." > ${meta.id}/README.txt
    """
}
