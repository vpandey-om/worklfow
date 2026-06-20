# Downstream Output Contracts

Every atomic module writes:

```text
metadata.json
command.txt
stdout.log
stderr.log
checksums.sha256
```

Main outputs:

```text
01_merge_featurecounts/
  merged_raw_counts.csv
  merge_summary.json
  sample_column_mapping.csv

02_validate_inputs/
  validated_counts.csv
  validated_metadata.csv
  validated_contrasts.csv
  validation_report.json
  sample_id_mapping.csv
  gene_mapping_validation.csv

03_filter_low_expression/
  filtered_counts.csv
  removed_genes.csv
  gene_filter_summary.csv
  filter_parameters.json

04_normalize_transform/
  normalized_counts.csv
  vst_expression.csv
  size_factors.csv
  transformation_summary.json

05_pca/
  pca_scores.csv
  pca_loadings.csv
  pca_explained_variance.csv
  pca_parameters.json

06_umap/
  umap_coordinates.csv
  umap_parameters.json

08_plsda/
  plsda_scores.csv
  plsda_loadings.csv
  plsda_vip_scores.csv
  plsda_cross_validation.csv
  plsda_permutation_test.csv
  plsda_summary.json
```

No plots, HTML reports, PDF reports, Dash layouts, or Plotly figures are produced in this phase.

DESeq2 and OPLS-DA entry points are reserved and fail clearly until their R dependencies are added.
