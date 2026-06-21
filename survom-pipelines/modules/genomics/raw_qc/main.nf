process GENOMICS_RAW_QC {
    tag "genomics_raw_qc"
    label 'process_low'
    publishDir "${params.outdir}/03_raw_qc", mode: 'copy'

    input:
    path manifest

    output:
    path "raw_qc_metrics.csv", emit: metrics
    path "raw_qc_summary.json", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py raw-qc \\
      --manifest ${manifest} \\
      --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
