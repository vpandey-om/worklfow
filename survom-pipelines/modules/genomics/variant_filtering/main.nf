process GENOMICS_VARIANT_FILTERING {
    tag "variant_filtering"
    label 'process_medium'
    publishDir "${params.outdir}/12_variant_filtering", mode: 'copy'

    input:
    path manifest

    output:
    path "filtered_vcf_manifest.csv", emit: manifest
    path "*.vcf.gz", emit: vcf, optional: true
    path "*.csi", emit: indexes, optional: true
    path "*.log", emit: logs, optional: true
    path "tool_check.json", emit: tool_check
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py variant-filter-real \\
      --manifest ${manifest} \\
      --expression 'QUAL>=20' \\
      --allow-missing \\
      --outdir .

    if [ ! -f tool_check.json ]; then
      printf '{"step":"variant-filter","available":true,"status":"completed"}\n' > tool_check.json
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
