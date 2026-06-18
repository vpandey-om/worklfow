process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/03_raw_qc/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"),  emit: zip
    path "versions.yml", emit: versions

    script:
    def adapterArg = params.custom_adapters ? "--adapters ${params.custom_adapters}" : ''
    def contaminantArg = params.custom_contaminants ? "--contaminants ${params.custom_contaminants}" : ''
    """
    fastqc \
      --threads ${task.cpus} \
      --outdir . \
      --extract \
      ${adapterArg} \
      ${contaminantArg} \
      ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
