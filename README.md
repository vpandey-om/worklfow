# SurvOm Multi-Omics Demo

This repository contains the Dash student app and the modular Nextflow workflow builder used for RNA-seq demo workflows.

## Run All Checks

Run one command from the repository root:

```bash
./scripts/run_all_tests_and_workflows.sh
```

The script detects the checks that exist in this repo, writes logs under `.test-runs/<timestamp>/`, and prints a final PASS/FAIL summary.

Full step-by-step notes are in:

```text
docs/how_to_run_tests_and_demo_workflows.md
```

Command-line workflow examples with custom input/output paths and params are in:

```text
docs/command_line_workflow_examples.md
```

Student/project result tracking notes are in:

```text
docs/student_project_result_tracking.md
```

Downstream RNA-seq CLI and output docs are in:

```text
survom-pipelines/downstream/docs/downstream_cli.md
survom-pipelines/downstream/docs/downstream_outputs.md
survom-pipelines/downstream/docs/downstream_resume_and_versions.md
```

Airway count-data downstream example:

```text
docs/airway_downstream_countdata_example.md
```

## Run a Workflow Locally

Students can run a workflow without opening the Dash app:

```bash
./run_workflow_local.sh --workflow trim_only --input /path/to/input_files --output /path/to/output_dir
```

Override workflow parameters with repeatable `--param key=value` arguments:

```bash
./run_workflow_local.sh \
  --workflow trim_only \
  --input /path/to/input_files \
  --output /tmp/my_trim_output \
  --param quality_threshold=25 \
  --param minimum_read_length=30
```

The input directory can contain either a manifest or FASTQ files.

Useful manifest names:

```text
raw_manifest.csv
samplesheet.csv
trim_manifest.tsv
trimmed_fastq_manifest.tsv
```

If no manifest is found, the runner tries to build one from FASTQ files using common names like `sample_R1.fastq.gz` and `sample_R2.fastq.gz`.

The runner creates:

```text
runs/<run_id>/
  state.json
  inputs/
  outputs/
  logs/
```

Resume a failed or interrupted run later:

```bash
./run_workflow_local.sh --resume <run_id>
```

Completed runs are not rerun if their saved outputs are still present. Failed runs are relaunched with Nextflow `-resume`, so completed Nextflow tasks can be reused.

Common workflow names include:

```text
raw_qc_only
trim_only
post_trim_qc_only
strandedness_only
salmon_quant_only
star_align_only
hisat2_align_only
salmon_count_matrix
star_count_matrix
qc_trim_strandedness
```
