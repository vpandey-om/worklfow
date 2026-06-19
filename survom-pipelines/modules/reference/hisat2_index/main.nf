process BUILD_HISAT2_INDEX {
    tag "hisat2_index"
    label 'process_high'
    publishDir "${params.outdir}/06A_reference", mode: 'copy'

    input:
    path genome_fasta

    output:
    path "hisat2_index", emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p hisat2_index
    case "${genome_fasta}" in
      *.gz) gzip -cd "${genome_fasta}" > genome.fa ;;
      *) cp "${genome_fasta}" genome.fa ;;
    esac
    hisat2-build -p ${task.cpus} genome.fa hisat2_index/genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | head -1 | sed 's/.*version //')
    END_VERSIONS
    """
}
