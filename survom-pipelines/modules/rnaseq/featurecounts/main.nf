process FEATURECOUNTS {
    tag "featurecounts"
    label 'process_medium'
    publishDir "${params.outdir}/09_star/featurecounts", mode: 'copy'

    input:
    path bams
    path annotation_gtf

    output:
    path "gene_count_matrix.tsv", emit: counts
    path "featurecounts_summary.txt", emit: summary
    path "featurecounts_status.tsv", emit: status
    path "featurecounts.log", emit: log

    script:
    def strand = params.featurecounts_strand ?: 0
    """
    case "${annotation_gtf}" in
      *.gz) gzip -cd "${annotation_gtf}" > genes.gtf ;;
      *) cp "${annotation_gtf}" genes.gtf ;;
    esac

    printf "status\\treason\\n" > featurecounts_status.tsv
    set +e
    featureCounts \\
      -T ${task.cpus} \\
      -s ${strand} \\
      -a genes.gtf \\
      -o featurecounts_raw.txt \\
      ${bams} \\
      > featurecounts.log 2>&1
    exit_code=\$?
    if [ "\${exit_code}" -ne 0 ] && grep -qi "Paired-end reads were detected" featurecounts.log; then
      featureCounts \\
        -T ${task.cpus} \\
        -s ${strand} \\
        -p \\
        --countReadPairs \\
        -a genes.gtf \\
        -o featurecounts_raw.txt \\
        ${bams} \\
        >> featurecounts.log 2>&1
      exit_code=\$?
    fi
    set -e

    if [ "\${exit_code}" -eq 0 ]; then
      cp featurecounts_raw.txt gene_count_matrix.tsv
      cp featurecounts_raw.txt.summary featurecounts_summary.txt
      printf "completed\\tfeatureCounts completed\\n" >> featurecounts_status.tsv
    else
      printf "Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\n" > gene_count_matrix.tsv
      printf "Status\\t%s\\n" "${bams}" > featurecounts_summary.txt
      printf "needs_review\\tfeatureCounts exited non-zero; inspect featurecounts.log\\n" >> featurecounts_status.tsv
    fi
    """
}
