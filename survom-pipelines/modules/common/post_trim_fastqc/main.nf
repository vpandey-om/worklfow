process POST_TRIM_FASTQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/05_post_trim_qc/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"),  emit: zip
    path "*.post_trim_qc_manifest.tsv", emit: manifest
    path "versions.yml", emit: versions

    script:
    def read_list = reads instanceof List ? reads : [reads]
    def fastq1 = read_list[0]?.name ?: ''
    def fastq2 = read_list.size() > 1 ? read_list[1].name : ''
    def single_end = meta.single_end ? 'true' : 'false'
    def strandedness = meta.strandedness ?: 'unknown'
    """
    fastqc \
      --threads ${task.cpus} \
      --outdir . \
      --extract \
      ${reads}

    printf "sample_id\tfastq_1\tfastq_2\tsingle_end\tstrandedness\tfastqc_html\tfastqc_zip\n" > ${meta.id}.post_trim_qc_manifest.tsv
    for zip_file in *_fastqc.zip; do
      sample_id=\${zip_file%_fastqc.zip}
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${meta.id}" "${fastq1}" "${fastq2}" "${single_end}" "${strandedness}" "\${sample_id}_fastqc.html" "\${zip_file}" >> ${meta.id}.post_trim_qc_manifest.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
