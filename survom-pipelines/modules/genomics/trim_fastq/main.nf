process GENOMICS_TRIM_FASTQ {
    tag "genomics_trim"
    label 'process_low'
    publishDir "${params.outdir}/04_trim", mode: 'copy'

    input:
    path manifest

    output:
    path "trimmed_fastq_manifest.csv", emit: manifest
    path "trim_report.json", emit: report
    path "trimmed_fastq/*", emit: reads
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py trim-fastq \\
      --manifest ${manifest} \\
      --quality-cutoff ${params.quality_cutoff ?: 20} \\
      --min-length ${params.min_length ?: 20} \\
      --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
