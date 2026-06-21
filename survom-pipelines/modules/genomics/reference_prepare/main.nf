process GENOMICS_REFERENCE_PREPARE {
    tag "reference"
    label 'process_low'
    publishDir "${params.outdir}/02_reference_prepare", mode: 'copy'

    input:
    path genome_fasta

    output:
    path "reference_bundle.json", emit: bundle
    path "*.fai", emit: fai
    path "*.dict", emit: dict
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py reference-prepare \\
      --genome-fasta ${genome_fasta} \\
      --reference-name ${params.reference_name ?: 'custom_reference'} \\
      --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
