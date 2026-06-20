# Airway Count-Data Downstream Example

This example starts after featureCounts. It uses the airway count matrix already saved in this repo.

Data folder:

```text
testdatasets/countdata/test_data/airway/
  airway_raw_counts.csv
  airway_sample_metadata.csv
  airway_contrasts.csv
```

Start from repo root:

```bash
cd /data/shared/vikash/mult-omics
```

## Run All Downstream Atomic Steps

```bash
./scripts/run_airway_downstream_atomic.sh
```

This runs:

```text
merge_featurecounts
validate_downstream_inputs
filter_low_expression
normalize_transform
pca
umap
differential_expression
plsda
go_enrichment
pathway_enrichment
ranked_gsea
write_state
status
```

Expected final line:

```text
Overall: PASS
```

Outputs are saved under:

```text
.test-runs/airway_downstream_<timestamp>/
```

## Run In Dash

Restart the Dash app after code changes, then open:

```text
http://127.0.0.1:8056/upload/
```

Select:

```text
Omics type: Bulk RNA-seq
Workflow: Downstream demo - Airway count matrix
```

Then click:

```text
Run analysis
```

This uses the bundled airway files automatically and writes CSV/JSON outputs. No plots are generated.

For student-uploaded count data, choose:

```text
Workflow: Downstream custom - Count matrix CSV outputs
```

Upload the count matrix, metadata, and contrasts CSV/TSV files using:

```text
Upload sample sheet / metadata
```

Select those files in the Uploaded files checklist. Dash will try to auto-fill the downstream fields by filename:

```text
counts / matrix / raw_counts
metadata / sample
contrast / contrasts
```

Students can also paste absolute server paths into the downstream input fields.

## Run One Atomic Step At A Time In Dash

The Workflow dropdown also has one-by-one downstream atomic workflows:

```text
Atomic DS 01 - Merge count matrix
Atomic DS 02 - Validate downstream inputs
Atomic DS 03 - Filter low-expression genes
Atomic DS 04 - Normalize and transform
Atomic DS 05 - PCA CSV outputs
Atomic DS 06 - UMAP CSV outputs
Atomic DS 07 - Differential expression CSV outputs
Atomic DS 08 - PLS-DA CSV outputs
Atomic DS 10 - GO enrichment CSV outputs
Atomic DS 11 - Pathway enrichment CSV outputs
Atomic DS 12 - Ranked GSEA CSV outputs
```

Recommended test order:

```text
1. Merge count matrix
2. Validate downstream inputs
3. Filter low-expression genes
4. Normalize and transform
5. PCA or UMAP or PLS-DA
6. Differential expression
7. GO/pathway enrichment or ranked GSEA
```

Dash tries to reuse previous outputs from the same student/project. For example:

```text
Validate -> uses merged_raw_counts.csv from Merge
Filter -> uses validated_counts.csv and validated_metadata.csv
Normalize -> uses filtered_counts.csv
PCA/UMAP/PLS-DA -> use vst_expression.csv
Enrichment/GSEA -> use ranked_genes.csv and significant_genes.csv from differential expression
```

You can still paste paths manually if you want to test a single step with external files.

## GO And Pathway Enrichment: What Algorithm Is Used?

The current teaching implementation uses local GMT files only.

```text
significant genes + tested gene universe + GMT gene sets
-> one-sided hypergeometric over-representation test
-> Benjamini-Hochberg adjusted p-values
-> CSV outputs
```

For ranked GSEA, the current teaching implementation is lightweight:

```text
ranked genes + GMT gene sets
-> mean rank statistic for genes present in each set
-> leading-edge gene CSV
```

This is good for testing and learning the workflow shape. For publication-grade enrichment later, plug in real MSigDB/Reactome/GO gene-set files and, if needed, a production GSEA method.

## Ensembl IDs, Symbols, Organism, And Mapping

The airway example uses human Ensembl gene IDs:

```text
ENSG00000003402
ENSG00000003096
```

If your GMT file also contains Ensembl IDs, no mapping is needed. Use:

```bash
survom-pipelines/downstream/bin/downstream_cli.py go-enrichment \
  --significant-genes /path/to/significant_genes.csv \
  --ranked-genes /path/to/ranked_genes.csv \
  --universe /path/to/filtered_counts.csv \
  --gmt testdatasets/countdata/test_data/airway/airway_human_ensembl_demo_sets.gmt \
  --organism human \
  --input-id-type ensembl_id \
  --gene-set-id-type ensembl_id \
  --outdir /tmp/airway_go_ensembl
```

If your counts use Ensembl IDs but your pathway GMT uses symbols, provide a mapping CSV with columns like:

```text
gene_id,gene_symbol,ensembl_id,entrez_id
```

Example:

```bash
survom-pipelines/downstream/bin/downstream_cli.py pathway-enrichment \
  --significant-genes /path/to/significant_genes.csv \
  --ranked-genes /path/to/ranked_genes.csv \
  --universe /path/to/filtered_counts.csv \
  --gmt testdatasets/countdata/test_data/airway/airway_human_symbol_demo_sets.gmt \
  --gene-mapping testdatasets/countdata/test_data/airway/airway_human_symbol_mapping_demo.csv \
  --organism human \
  --input-id-type ensembl_id \
  --gene-set-id-type gene_symbol \
  --outdir /tmp/airway_pathway_symbol
```

Dash has the same idea: for GO/pathway/GSEA atomics, paste or upload the GMT file. If the GMT uses a different ID type from the count matrix, also provide the optional gene mapping CSV.

## Choose Your Own Output Folder

```bash
./scripts/run_airway_downstream_atomic.sh \
  --outdir /tmp/airway_downstream_test
```

## Use A Different Count Dataset

Your folder should contain:

```text
my_counts.csv
my_metadata.csv
my_contrasts.csv
```

Counts format:

```csv
gene_id,sample1,sample2,sample3,sample4
GENE1,10,20,30,40
GENE2,0,5,8,12
```

Metadata format:

```csv
sample_id,condition
sample1,control
sample2,control
sample3,treated
sample4,treated
```

Contrasts format:

```csv
contrast_id,numerator,denominator,design_column
treated_vs_control,treated,control,condition
```

Run atomic commands manually:

```bash
python survom-pipelines/downstream/bin/downstream_cli.py merge-featurecounts \
  --counts /path/to/my_counts.csv \
  --outdir /tmp/my_downstream/01_merge_featurecounts

python survom-pipelines/downstream/bin/downstream_cli.py validate \
  --counts /tmp/my_downstream/01_merge_featurecounts/merged_raw_counts.csv \
  --metadata /path/to/my_metadata.csv \
  --contrasts /path/to/my_contrasts.csv \
  --outdir /tmp/my_downstream/02_validate_inputs
```

Continue with:

```bash
python survom-pipelines/downstream/bin/downstream_cli.py filter \
  --counts /tmp/my_downstream/02_validate_inputs/validated_counts.csv \
  --metadata /tmp/my_downstream/02_validate_inputs/validated_metadata.csv \
  --outdir /tmp/my_downstream/03_filter_low_expression

python survom-pipelines/downstream/bin/downstream_cli.py normalize \
  --counts /tmp/my_downstream/03_filter_low_expression/filtered_counts.csv \
  --metadata /tmp/my_downstream/02_validate_inputs/validated_metadata.csv \
  --outdir /tmp/my_downstream/04_normalize_transform

python survom-pipelines/downstream/bin/downstream_cli.py pca \
  --expression /tmp/my_downstream/04_normalize_transform/vst_expression.csv \
  --metadata /tmp/my_downstream/02_validate_inputs/validated_metadata.csv \
  --outdir /tmp/my_downstream/05_pca
```

## Current Dash Status

The downstream layer is backend-ready and command-line ready.

Dash can call the backend service in:

```text
demo_dash_app/services/downstream_service.py
```

No downstream Dash forms or plots are added yet. That is intentional for this implementation phase, because the downstream request asked for CSV/TSV/JSON outputs only and no plotting.
