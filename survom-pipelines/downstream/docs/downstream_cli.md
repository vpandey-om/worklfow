# Downstream RNA-seq CLI

This downstream layer starts after featureCounts or a merged gene-by-sample count matrix.

It does not run FASTQ QC, trimming, alignment, reference indexing, BAM processing, or featureCounts.

## Atomic Commands

From repo root:

```bash
cd /data/shared/vikash/mult-omics
export PYTHONPATH=/data/shared/vikash/mult-omics/survom-pipelines/downstream/bin
```

Merge featureCounts:

```bash
python survom-pipelines/downstream/bin/downstream_cli.py merge-featurecounts \
  --counts survom-pipelines/downstream/tests/fixtures/featurecounts_raw.txt \
  --outdir /tmp/downstream_demo/01_merge_featurecounts
```

Validate:

```bash
python survom-pipelines/downstream/bin/downstream_cli.py validate \
  --counts /tmp/downstream_demo/01_merge_featurecounts/merged_raw_counts.csv \
  --metadata survom-pipelines/downstream/tests/fixtures/sample_metadata.csv \
  --contrasts survom-pipelines/downstream/tests/fixtures/contrasts.csv \
  --outdir /tmp/downstream_demo/02_validate_inputs
```

Filter:

```bash
python survom-pipelines/downstream/bin/downstream_cli.py filter \
  --counts /tmp/downstream_demo/02_validate_inputs/validated_counts.csv \
  --metadata /tmp/downstream_demo/02_validate_inputs/validated_metadata.csv \
  --outdir /tmp/downstream_demo/03_filter_low_expression
```

Normalize/transform:

```bash
python survom-pipelines/downstream/bin/downstream_cli.py normalize \
  --counts /tmp/downstream_demo/03_filter_low_expression/filtered_counts.csv \
  --metadata /tmp/downstream_demo/02_validate_inputs/validated_metadata.csv \
  --outdir /tmp/downstream_demo/04_normalize_transform
```

PCA:

```bash
python survom-pipelines/downstream/bin/downstream_cli.py pca \
  --expression /tmp/downstream_demo/04_normalize_transform/vst_expression.csv \
  --metadata /tmp/downstream_demo/02_validate_inputs/validated_metadata.csv \
  --outdir /tmp/downstream_demo/05_pca
```

UMAP:

```bash
python survom-pipelines/downstream/bin/downstream_cli.py umap \
  --expression /tmp/downstream_demo/04_normalize_transform/vst_expression.csv \
  --metadata /tmp/downstream_demo/02_validate_inputs/validated_metadata.csv \
  --outdir /tmp/downstream_demo/06_umap
```

PLS-DA:

```bash
python survom-pipelines/downstream/bin/downstream_cli.py plsda \
  --expression /tmp/downstream_demo/04_normalize_transform/vst_expression.csv \
  --metadata /tmp/downstream_demo/02_validate_inputs/validated_metadata.csv \
  --outdir /tmp/downstream_demo/08_plsda
```

## Full Downstream Nextflow

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines/downstream

nextflow -C nextflow_downstream.config run main_downstream.nf \
  -profile test \
  --featurecounts_input tests/fixtures/featurecounts_raw.txt \
  --metadata tests/fixtures/sample_metadata.csv \
  --contrasts tests/fixtures/contrasts.csv \
  --run_id downstream_demo_001 \
  --runs_dir /tmp/survom_downstream_runs \
  -resume
```

Resume:

```bash
nextflow -C nextflow_downstream.config run main_downstream.nf \
  -profile test \
  --featurecounts_input tests/fixtures/featurecounts_raw.txt \
  --metadata tests/fixtures/sample_metadata.csv \
  --contrasts tests/fixtures/contrasts.csv \
  --run_id downstream_demo_001 \
  --runs_dir /tmp/survom_downstream_runs \
  -resume
```

Status:

```bash
python bin/downstream_cli.py status --run-id downstream_demo_001 --runs-dir /tmp/survom_downstream_runs
```

## Tests

Dependency-free smoke:

```bash
python survom-pipelines/downstream/tests/run_downstream_smoke.py
```

Pytest, if installed:

```bash
cd survom-pipelines/downstream
pytest -q tests
```

All repo checks:

```bash
cd /data/shared/vikash/mult-omics
./scripts/run_all_tests_and_workflows.sh --quick
```
