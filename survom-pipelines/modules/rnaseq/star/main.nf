process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/06_strandedness/star", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path star_index

    output:
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: bam
    path "*.alignment_manifest.tsv", emit: manifest
    path "*.Log.final.out", emit: log
    path "versions.yml", emit: versions

    script:
    def readFilesCommand = reads[0].toString().endsWith(".gz") ? "--readFilesCommand zcat" : ""
    """
    STAR \
      --genomeDir ${star_index} \
      --readFilesIn ${reads} \
      ${readFilesCommand} \
      --runThreadN ${task.cpus} \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix ${meta.id}.

    printf "sample_id\\tbam\\tbai\\ttool\\n" > ${meta.id}.alignment_manifest.tsv
    printf "%s\\t%s\\t%s\\tstar\\n" "${meta.id}" "${meta.id}.Aligned.sortedByCoord.out.bam" "" >> ${meta.id}.alignment_manifest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSIONS
    """
}
