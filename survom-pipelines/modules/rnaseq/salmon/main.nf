process SALMON {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/07_salmon", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path salmon_index

    output:
    tuple val(meta), path("${meta.id}"), emit: quant
    path "*.quant_manifest.tsv", emit: manifest
    path "*.salmon_quant_status.tsv", emit: status
    path "*.salmon_quant.log", emit: log
    path "versions.yml", emit: versions

    script:
    def readArgs = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def libType = params.approved_library_type ?: params.inferred_library_type ?: "A"
    """
    set +e
    salmon quant \\
      -i ${salmon_index} \\
      -l ${libType} \\
      ${readArgs} \\
      --validateMappings \\
      --minAssignedFrags 1 \\
      -p ${task.cpus} \\
      -o ${meta.id} \\
      > ${meta.id}.salmon_quant.log 2>&1
    salmon_exit=\$?
    set -e

    status="completed"
    reason="Salmon quantification completed."
    if [ "\${salmon_exit}" -ne 0 ]; then
      status="needs_review"
      reason="Salmon quantification exited non-zero; inspect salmon_quant.log. This often means the reference/index does not match the reads or too few fragments were assigned."
      mkdir -p ${meta.id}
      if [ ! -s ${meta.id}/quant.sf ]; then
        printf "Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n" > ${meta.id}/quant.sf
      fi
    fi

    printf "sample_id\\tstatus\\tsalmon_exit_code\\treason\\tquant_sf\\tlog\\n" > ${meta.id}.salmon_quant_status.tsv
    printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "${meta.id}" "\${status}" "\${salmon_exit}" "\${reason}" "${meta.id}/quant.sf" "${meta.id}.salmon_quant.log" >> ${meta.id}.salmon_quant_status.tsv
    printf "sample_id\\tquant_dir\\tquant_sf\\ttool\\n" > ${meta.id}.quant_manifest.tsv
    printf "%s\\t%s\\t%s\\tsalmon\\n" "${meta.id}" "${meta.id}" "${meta.id}/quant.sf" >> ${meta.id}.quant_manifest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}
