process FALCO {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/03_raw_qc/falco", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.txt"),  emit: txt
    path "versions.yml", emit: versions

    script:
    """
    falco \
      --outdir . \
      ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        falco: \$(falco --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
