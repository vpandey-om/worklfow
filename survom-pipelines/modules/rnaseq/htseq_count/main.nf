process HTSEQ_COUNT {
    tag "htseq_count"
    label 'process_medium'
    publishDir "${params.outdir}/09_counts/htseq", mode: 'copy'

    input:
    path alignments
    path annotation_gtf

    output:
    path "htseq_count_matrix.tsv", emit: counts
    path "*.htseq_count.txt", emit: per_sample_counts
    path "*.htseq_count.log", emit: logs
    path "htseq_count_status.tsv", emit: status
    path "versions.yml", emit: versions

    script:
    def stranded = params.htseq_stranded ?: "no"
    """
    case "${annotation_gtf}" in
      *.gz) gzip -cd "${annotation_gtf}" > genes.gtf ;;
      *) cp "${annotation_gtf}" genes.gtf ;;
    esac

    printf "sample_id\\tstatus\\texit_code\\treason\\n" > htseq_count_status.tsv
    printf "gene_id" > htseq_count_matrix.tsv

    first_sample=""
    for aln in ${alignments}; do
    sample=\$(basename "\${aln}")
    sample=\${sample%.sorted.bam}
    sample=\${sample%.bam}
    sample=\${sample%.sam}
    sample=\${sample%.Aligned.sortedByCoord.out}
    sample=\${sample%.hisat2}
      printf "\\t%s" "\${sample}" >> htseq_count_matrix.tsv

      format="sam"
      case "\${aln}" in
        *.bam) format="bam" ;;
      esac

      set +e
      htseq-count \\
        -f "\${format}" \\
        -r pos \\
        -s "${stranded}" \\
        -t exon \\
        -i gene_id \\
        "\${aln}" genes.gtf \\
        > "\${sample}.htseq_count.txt" \\
        2> "\${sample}.htseq_count.log"
      exit_code=\$?
      set -e
      if [ "\${exit_code}" -eq 0 ]; then
        printf "%s\\tcompleted\\t%s\\tHTSeq count completed\\n" "\${sample}" "\${exit_code}" >> htseq_count_status.tsv
      else
        printf "%s\\tneeds_review\\t%s\\tHTSeq count exited non-zero; inspect log\\n" "\${sample}" "\${exit_code}" >> htseq_count_status.tsv
        touch "\${sample}.htseq_count.txt"
      fi

      if [ -z "\${first_sample}" ]; then
        first_sample="\${sample}"
        awk 'BEGIN{FS=OFS="\\t"} \$1 !~ /^__/ {print \$1}' "\${sample}.htseq_count.txt" > .genes.tmp
        awk 'BEGIN{FS=OFS="\\t"} \$1 !~ /^__/ {print \$2+0}' "\${sample}.htseq_count.txt" > ".\${sample}.counts.tmp"
      else
        awk 'BEGIN{FS=OFS="\\t"} \$1 !~ /^__/ {print \$2+0}' "\${sample}.htseq_count.txt" > ".\${sample}.counts.tmp"
      fi
    done
    printf "\\n" >> htseq_count_matrix.tsv

    if [ -s .genes.tmp ]; then
      paste .genes.tmp .*.counts.tmp >> htseq_count_matrix.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq: \$(htseq-count --version 2>&1 | head -1)
    END_VERSIONS
    """
}
