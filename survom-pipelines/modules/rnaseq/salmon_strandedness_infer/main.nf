process SALMON_STRANDEDNESS_INFERENCE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/06_strandedness", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path salmon_index

    output:
    tuple val(meta), path("*.strandedness_call.json"), emit: call
    tuple val(meta), path("*.strandedness_evidence.tsv"), emit: evidence
    path "*.salmon_inference.log", emit: log
    path "*.strandness_manifest.tsv", emit: manifest
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    def read_list = reads instanceof List ? reads : [reads]
    def fastq1 = read_list[0]?.name ?: ''
    def fastq2 = read_list.size() > 1 ? read_list[1].name : ''
    def single_end = meta.single_end ? 'true' : 'false'
    def strandedness = meta.strandedness ?: 'unknown'
    def readArgs = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    set +e
    salmon quant \\
      -i ${salmon_index} \\
      -l A \\
      ${readArgs} \\
      --validateMappings \\
      --numPreAuxModelSamples ${params.strandedness_inference_reads} \\
      --minAssignedFrags 1 \\
      -p ${task.cpus} \\
      -o ${prefix}.salmon_inference \\
      > ${prefix}.salmon_inference.log 2>&1
    salmon_exit=\$?
    set -e

    meta_path="${prefix}.salmon_inference/aux_info/meta_info.json"
    inferred=""
    status="needs_review"
    reason="Salmon did not report a single confident inferred library type."
    if [ "\${salmon_exit}" -ne 0 ]; then
      status="needs_review"
      reason="Salmon inference exited non-zero; inspect salmon_inference_log. This often means the reference/index does not match the reads or too few fragments were assigned."
    elif [ -s "\${meta_path}" ]; then
      inferred=\$(grep -m 1 -E '"expected_format"|"expectedFormat"' "\${meta_path}" | sed 's/.*: *"//; s/".*//' || true)
      if [ -n "\${inferred}" ]; then
        status="completed"
        reason="Salmon reported expected_format."
      fi
    fi
    needs_review="true"
    if [ "\${status}" = "completed" ]; then
      needs_review="false"
    fi

    {
      printf "{\\n"
      printf "  \\"sample_id\\": \\"%s\\",\\n" "${prefix}"
      printf "  \\"method\\": \\"salmon_auto\\",\\n"
      printf "  \\"status\\": \\"%s\\",\\n" "\${status}"
      printf "  \\"inferred_library_type\\": \\"%s\\",\\n" "\${inferred}"
      printf "  \\"needs_review\\": %s,\\n" "\${needs_review}"
      printf "  \\"approved\\": false,\\n"
      printf "  \\"approved_library_type\\": null,\\n"
      printf "  \\"reason\\": \\"%s\\",\\n" "\${reason}"
      printf "  \\"salmon_library_type_argument\\": \\"A\\",\\n"
      printf "  \\"downstream_param\\": \\"inferred_library_type\\"\\n"
      printf "}\\n"
    } > ${prefix}.strandedness_call.json

    printf "key\tvalue\n" > ${prefix}.strandedness_evidence.tsv
    printf "sample_id\t%s\n" "${prefix}" >> ${prefix}.strandedness_evidence.tsv
    printf "method\tsalmon_auto\n" >> ${prefix}.strandedness_evidence.tsv
    printf "status\t%s\n" "\${status}" >> ${prefix}.strandedness_evidence.tsv
    printf "inferred_library_type\t%s\n" "\${inferred}" >> ${prefix}.strandedness_evidence.tsv
    printf "reason\t%s\n" "\${reason}" >> ${prefix}.strandedness_evidence.tsv
    printf "salmon_exit_code\t%s\n" "\${salmon_exit}" >> ${prefix}.strandedness_evidence.tsv
    printf "meta_info_json\t%s\n" "\${meta_path}" >> ${prefix}.strandedness_evidence.tsv

    printf "sample_id\tfastq_1\tfastq_2\tsingle_end\tstrandedness\tmethod\tstatus\tinferred_library_type\tcall_json\tevidence_tsv\tlog\n" > ${prefix}.strandness_manifest.tsv
    printf "%s\t%s\t%s\t%s\t%s\tsalmon_auto\t%s\t%s\t%s\t%s\t%s\n" "${prefix}" "${fastq1}" "${fastq2}" "${single_end}" "${strandedness}" "\${status}" "\${inferred}" "${prefix}.strandedness_call.json" "${prefix}.strandedness_evidence.tsv" "${prefix}.salmon_inference.log" >> ${prefix}.strandness_manifest.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}
