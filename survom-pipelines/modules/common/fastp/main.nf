process FASTP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/04_trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.fastp.json"),       emit: json
    tuple val(meta), path("*.fastp.html"),       emit: html
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    def platform = (meta.platform ?: '').toString().toLowerCase()
    def trimPolyG = params.trim_poly_g == true || params.trim_poly_g == 'true' || params.trim_poly_g == 'on'
    def autoPolyG = params.trim_poly_g == 'auto' && (platform.contains('nextseq') || platform.contains('novaseq'))
    def polyG = (trimPolyG || autoPolyG) ? '--trim_poly_g' : ''
    def polyX = params.trim_poly_x ? '--trim_poly_x' : ''
    def adapterR1 = params.adapter_sequence_r1 ? "--adapter_sequence ${params.adapter_sequence_r1}" : ''
    def adapterR2 = params.adapter_sequence_r2 ? "--adapter_sequence_r2 ${params.adapter_sequence_r2}" : ''
    def adapterFasta = params.adapter_fasta ? "--adapter_fasta ${params.adapter_fasta}" : ''
    def trimFrontR1Value = params.trim_front_r1 ? params.trim_front_r1 as Integer : 0
    def trimFrontR2Value = params.trim_front_r2 ? params.trim_front_r2 as Integer : 0
    def trimFrontR1 = trimFrontR1Value > 0 ? "--trim_front1 ${trimFrontR1Value}" : ''
    def trimFrontR2 = trimFrontR2Value > 0 ? "--trim_front2 ${trimFrontR2Value}" : ''
    def detectPE = !meta.single_end && params.detect_adapter_for_pe && !params.adapter_sequence_r1 && !params.adapter_sequence_r2 && !params.adapter_fasta ? '--detect_adapter_for_pe' : ''

    if (meta.single_end) {
        """
        fastp \
          --in1 ${reads[0]} \
          --out1 ${prefix}.trimmed.fastq.gz \
          --qualified_quality_phred ${params.quality_threshold} \
          --length_required ${params.minimum_read_length} \
          ${polyG} ${polyX} ${adapterR1} ${adapterFasta} ${trimFrontR1} \
          --thread ${task.cpus} \
          --html ${prefix}.fastp.html \
          --json ${prefix}.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    } else {
        """
        fastp \
          --in1 ${reads[0]} \
          --in2 ${reads[1]} \
          --out1 ${prefix}_R1.trimmed.fastq.gz \
          --out2 ${prefix}_R2.trimmed.fastq.gz \
          --qualified_quality_phred ${params.quality_threshold} \
          --length_required ${params.minimum_read_length} \
          ${detectPE} ${polyG} ${polyX} ${adapterR1} ${adapterR2} ${adapterFasta} ${trimFrontR1} ${trimFrontR2} \
          --thread ${task.cpus} \
          --html ${prefix}.fastp.html \
          --json ${prefix}.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    }
}
