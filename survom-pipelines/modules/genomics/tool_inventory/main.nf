process GENOMICS_TOOL_INVENTORY {
    tag "tool_inventory"
    label 'process_low'
    publishDir "${params.outdir}/17_tool_inventory", mode: 'copy'

    output:
    path "tool_inventory.csv", emit: inventory_csv
    path "tool_inventory.json", emit: inventory_json
    path "versions.yml", emit: versions

    script:
    """
    python3 ${params.survom_pipeline_root}/bin/python/common/genomics_atomic.py tool-inventory \\
      --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
