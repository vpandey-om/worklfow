process DESEQ2 {
    tag "deseq2"
    label 'process_medium'

    input:
    path counts

    output:
    path "deseq2_results.tsv", emit: results

    script:
    """
    echo "gene_id\tlog2FoldChange\tpadj" > deseq2_results.tsv
    """
}
