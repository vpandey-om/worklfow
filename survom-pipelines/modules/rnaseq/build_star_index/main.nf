process BUILD_STAR_INDEX {
    tag "star_index"
    label 'process_high'
    publishDir "${params.outdir}/06A_reference", mode: 'copy'

    input:
    path genome_fasta
    path gtf

    output:
    path "star_index", emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p star_index
    case "${genome_fasta}" in
      *.gz) gzip -cd "${genome_fasta}" > genome.fa ;;
      *) cp "${genome_fasta}" genome.fa ;;
    esac
    case "${gtf}" in
      *.gz) gzip -cd "${gtf}" > genes.gtf ;;
      *) cp "${gtf}" genes.gtf ;;
    esac
    genome_bases=\$(
      awk '/^>/ { next } { total += length(\$0) } END { print total + 0 }' genome.fa
    )
    genome_sa_index_nbases=\$(
      awk -v n="\${genome_bases}" 'BEGIN {
        if (n < 1) { print 3; exit }
        value = int(log(n) / log(2) / 2 - 1)
        if (value < 3) value = 3
        if (value > 14) value = 14
        print value
      }'
    )
    STAR \\
      --runMode genomeGenerate \\
      --genomeDir star_index \\
      --genomeFastaFiles genome.fa \\
      --sjdbGTFfile genes.gtf \\
      --genomeSAindexNbases \${genome_sa_index_nbases} \\
      --runThreadN ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSIONS
    """
}
