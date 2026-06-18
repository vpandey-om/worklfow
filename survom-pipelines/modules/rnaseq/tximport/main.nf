process TXIMPORT {
    tag "tximport"
    label 'process_low'

    input:
    path quant_dirs

    output:
    path "tximport_counts.tsv", emit: counts

    script:
    """
    echo "gene_id\tcount" > tximport_counts.tsv
    """
}
