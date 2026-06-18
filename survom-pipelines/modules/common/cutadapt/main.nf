process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/04_trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.cutadapt.json"),    emit: json
    path "*.cutadapt.log", emit: log
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    def minOverlap = params.cutadapt_minimum_overlap ?: params.cutadapt_min_overlap
    if (meta.single_end) {
        """
        cutadapt \
          -j ${task.cpus} \
          -a ${params.adapter_sequence_r1} \
          -q ${params.quality_threshold} \
          -m ${params.minimum_read_length} \
          -e ${params.cutadapt_error_rate} \
          -O ${minOverlap} \
          --json ${prefix}.cutadapt.json \
          -o ${prefix}.trimmed.fastq.gz \
          ${reads[0]} \
          > ${prefix}.cutadapt.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        cutadapt \
          -j ${task.cpus} \
          -a ${params.adapter_sequence_r1} \
          -A ${params.adapter_sequence_r2} \
          -q ${params.quality_threshold} \
          -m ${params.minimum_read_length} \
          -e ${params.cutadapt_error_rate} \
          -O ${minOverlap} \
          --json ${prefix}.cutadapt.json \
          -o ${prefix}_R1.trimmed.fastq.gz \
          -p ${prefix}_R2.trimmed.fastq.gz \
          ${reads[0]} ${reads[1]} \
          > ${prefix}.cutadapt.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
}
