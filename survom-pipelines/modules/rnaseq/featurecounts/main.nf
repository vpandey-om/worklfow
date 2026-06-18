process FEATURECOUNTS {
    tag "featurecounts"
    label 'process_medium'

    input:
    path bams

    output:
    path "counts.tsv", emit: counts

    script:
    """
    echo "feature_id\tcount" > counts.tsv
    """
}
