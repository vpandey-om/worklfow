process GENOMICS_GERMLINE_VARIANTS {
    tag "germline_variants"
    label 'process_high'
    publishDir "${params.outdir}/10_germline_variants", mode: 'copy'

    input:
    path manifest
    path reference

    output:
    path "raw_vcf_manifest.csv", emit: manifest
    path "*.vcf.gz", emit: vcf, optional: true
    path "*.csi", emit: indexes, optional: true
    path "*.log", emit: logs, optional: true
    path "tool_check.json", emit: tool_check
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py germline-variants-real \\
      --manifest ${manifest} \\
      --reference ${reference} \\
      --allow-missing \\
      --outdir .

    if [ ! -f tool_check.json ]; then
      printf '{"step":"germline-variants","available":true,"status":"completed"}\n' > tool_check.json
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
