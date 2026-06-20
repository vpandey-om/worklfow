nextflow.enable.dsl=2

params.featurecounts_input = null
params.metadata = null
params.contrasts = null
params.gene_mapping = null
params.gmt = null
params.run_id = "downstream_demo"
params.runs_dir = "runs"
params.group_column = "condition"
params.cpm = 1.0
params.min_samples = null
params.pca_components = 5
params.seed = 42

process MERGE_FEATURECOUNTS {
    tag "merge_featurecounts"
    publishDir "${params.runs_dir}/${params.run_id}/results/01_merge_featurecounts", mode: 'copy'

    input:
    path featurecounts_input

    output:
    path "merged_raw_counts.csv", emit: counts
    path "merge_summary.json", emit: summary
    path "sample_column_mapping.csv", emit: mapping
    path "metadata.json"
    path "command.txt"
    path "stdout.log"
    path "stderr.log"
    path "checksums.sha256"

    script:
    def mode = featurecounts_input.isDirectory() ? "--input-dir ${featurecounts_input}" : "--counts ${featurecounts_input}"
    """
    python3 ${projectDir}/bin/downstream_cli.py merge-featurecounts \\
      ${mode} \\
      --outdir .
    """
}

process VALIDATE_DOWNSTREAM_INPUTS {
    tag "validate_downstream_inputs"
    publishDir "${params.runs_dir}/${params.run_id}/results/02_validate_inputs", mode: 'copy'

    input:
    path counts
    path metadata
    val contrasts
    val gene_mapping

    output:
    path "validated_counts.csv", emit: counts
    path "validated_metadata.csv", emit: metadata
    path "validated_contrasts.csv", emit: contrasts
    path "validation_report.json", emit: report
    path "sample_id_mapping.csv"
    path "gene_mapping_validation.csv"
    path "metadata.json"
    path "command.txt"
    path "stdout.log"
    path "stderr.log"
    path "checksums.sha256"

    script:
    def contrast_arg = contrasts ? "--contrasts ${contrasts}" : ""
    def mapping_arg = gene_mapping ? "--gene-mapping ${gene_mapping}" : ""
    """
    python3 ${projectDir}/bin/downstream_cli.py validate \\
      --counts ${counts} \\
      --metadata ${metadata} \\
      ${contrast_arg} \\
      ${mapping_arg} \\
      --group-column ${params.group_column} \\
      --outdir .
    """
}

process FILTER_LOW_EXPRESSION {
    tag "filter_low_expression"
    publishDir "${params.runs_dir}/${params.run_id}/results/03_filter_low_expression", mode: 'copy'

    input:
    path counts
    path metadata

    output:
    path "filtered_counts.csv", emit: counts
    path "removed_genes.csv", emit: removed
    path "gene_filter_summary.csv", emit: summary
    path "filter_parameters.json", emit: params_json
    path "metadata.json"
    path "command.txt"
    path "stdout.log"
    path "stderr.log"
    path "checksums.sha256"

    script:
    def min_samples_arg = params.min_samples ? "--min-samples ${params.min_samples}" : ""
    """
    python3 ${projectDir}/bin/downstream_cli.py filter \\
      --counts ${counts} \\
      --metadata ${metadata} \\
      --group-column ${params.group_column} \\
      --cpm ${params.cpm} \\
      ${min_samples_arg} \\
      --outdir .
    """
}

process NORMALIZE_AND_TRANSFORM {
    tag "normalize_transform"
    publishDir "${params.runs_dir}/${params.run_id}/results/04_normalize_transform", mode: 'copy'

    input:
    path counts
    path metadata

    output:
    path "normalized_counts.csv", emit: normalized
    path "vst_expression.csv", emit: vst
    path "size_factors.csv", emit: size_factors
    path "transformation_summary.json", emit: summary
    path "metadata.json"
    path "command.txt"
    path "stdout.log"
    path "stderr.log"
    path "checksums.sha256"

    script:
    """
    python3 ${projectDir}/bin/downstream_cli.py normalize \\
      --counts ${counts} \\
      --metadata ${metadata} \\
      --outdir .
    """
}

process PCA {
    tag "pca"
    publishDir "${params.runs_dir}/${params.run_id}/results/05_pca", mode: 'copy'

    input:
    path expression
    path metadata

    output:
    path "pca_scores.csv", emit: scores
    path "pca_loadings.csv", emit: loadings
    path "pca_explained_variance.csv", emit: variance
    path "pca_parameters.json", emit: params_json
    path "metadata.json"
    path "command.txt"
    path "stdout.log"
    path "stderr.log"
    path "checksums.sha256"

    script:
    """
    python3 ${projectDir}/bin/downstream_cli.py pca \\
      --expression ${expression} \\
      --metadata ${metadata} \\
      --components ${params.pca_components} \\
      --seed ${params.seed} \\
      --outdir .
    """
}

process UMAP {
    tag "umap"
    publishDir "${params.runs_dir}/${params.run_id}/results/06_umap", mode: 'copy'

    input:
    path expression
    path metadata

    output:
    path "umap_coordinates.csv", emit: coordinates
    path "umap_parameters.json", emit: params_json
    path "metadata.json"
    path "command.txt"
    path "stdout.log"
    path "stderr.log"
    path "checksums.sha256"

    script:
    """
    python3 ${projectDir}/bin/downstream_cli.py umap \\
      --expression ${expression} \\
      --metadata ${metadata} \\
      --seed ${params.seed} \\
      --outdir .
    """
}

process PLSDA {
    tag "plsda"
    publishDir "${params.runs_dir}/${params.run_id}/results/08_plsda", mode: 'copy'

    input:
    path expression
    path metadata

    output:
    path "plsda_scores.csv", emit: scores
    path "plsda_loadings.csv", emit: loadings
    path "plsda_vip_scores.csv", emit: vip
    path "plsda_cross_validation.csv", emit: cv
    path "plsda_permutation_test.csv", emit: permutation
    path "plsda_summary.json", emit: summary
    path "metadata.json"
    path "command.txt"
    path "stdout.log"
    path "stderr.log"
    path "checksums.sha256"

    script:
    """
    python3 ${projectDir}/bin/downstream_cli.py plsda \\
      --expression ${expression} \\
      --metadata ${metadata} \\
      --group-column ${params.group_column} \\
      --seed ${params.seed} \\
      --outdir .
    """
}

workflow {
    if (!params.featurecounts_input) {
        error "Missing --featurecounts_input"
    }
    if (!params.metadata) {
        error "Missing --metadata"
    }
    featurecounts_input = file(params.featurecounts_input, checkIfExists: true)
    metadata = file(params.metadata, checkIfExists: true)
    contrasts = params.contrasts ?: ''
    gene_mapping = params.gene_mapping ?: ''

    MERGE_FEATURECOUNTS(featurecounts_input)
    VALIDATE_DOWNSTREAM_INPUTS(
        MERGE_FEATURECOUNTS.out.counts,
        metadata,
        contrasts,
        gene_mapping
    )
    FILTER_LOW_EXPRESSION(VALIDATE_DOWNSTREAM_INPUTS.out.counts, VALIDATE_DOWNSTREAM_INPUTS.out.metadata)
    NORMALIZE_AND_TRANSFORM(FILTER_LOW_EXPRESSION.out.counts, VALIDATE_DOWNSTREAM_INPUTS.out.metadata)
    PCA(NORMALIZE_AND_TRANSFORM.out.vst, VALIDATE_DOWNSTREAM_INPUTS.out.metadata)
    UMAP(NORMALIZE_AND_TRANSFORM.out.vst, VALIDATE_DOWNSTREAM_INPUTS.out.metadata)
    PLSDA(NORMALIZE_AND_TRANSFORM.out.vst, VALIDATE_DOWNSTREAM_INPUTS.out.metadata)
}
