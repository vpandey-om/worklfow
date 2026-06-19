# How To Run Tests And Demo Workflows

This guide is for remembering the exact commands.

Start from the repository root:

```bash
cd /data/shared/vikash/mult-omics
```

## 1. Run Everything

Use this when you want one command that checks code, workflow generation, and real tiny demo Nextflow runs.

```bash
./scripts/run_all_tests_and_workflows.sh
```

This runs:

```text
python compile checks
python unittest local runner
pytest survom-pipelines, if pytest is installed
workflow smoke compile
workflow real Nextflow smoke
```

Expected final result:

```text
Overall: PASS
```

Meaning:

```text
Code compiles.
Local runner save/resume tests pass.
Workflow generator compiles all smoke workflows.
Tiny real FASTQ/reference Nextflow demo workflows run successfully.
```

Logs are saved in:

```text
.test-runs/<timestamp>/
```

Example:

```bash
ls -ltr .test-runs
ls .test-runs/$(ls -1 .test-runs | tail -1)
```

## 2. Run Quick Checks Only

Use this while editing code and you do not want to wait for Docker/Nextflow execution.

```bash
./scripts/run_all_tests_and_workflows.sh --quick
```

This checks code and workflow compilation, but skips real Nextflow runs.

Say this in reports:

```text
Quick code/unit/compile checks passed.
```

Do not call this a full workflow execution test.

## 3. Docker Permission Note

The full test uses Docker containers through Nextflow.

If you see:

```text
permission denied while trying to connect to the docker API at unix:///var/run/docker.sock
```

then the workflow test could not access Docker. That is an environment permission problem.

Try from the normal server shell where your user can run Docker:

```bash
docker run --rm hello-world
```

Then rerun:

```bash
./scripts/run_all_tests_and_workflows.sh
```

## 4. Run One Demo Workflow Directly

The demo smoke script creates tiny FASTQ/reference inputs automatically.

Go to the pipeline folder:

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
```

Run one workflow:

```bash
python scripts/smoke_workflows.py --case qc_trim_strandedness --run --keep /tmp/survom_demo_qc_trim_strandedness
```

Useful demo cases:

```bash
python scripts/smoke_workflows.py --case reference_only --run --keep /tmp/survom_demo_reference
python scripts/smoke_workflows.py --case raw_qc_only --run --keep /tmp/survom_demo_raw_qc
python scripts/smoke_workflows.py --case trim_only --run --keep /tmp/survom_demo_trim
python scripts/smoke_workflows.py --case post_trim_qc_only --run --keep /tmp/survom_demo_post_trim_qc
python scripts/smoke_workflows.py --case strandedness_only --run --keep /tmp/survom_demo_strandedness
python scripts/smoke_workflows.py --case salmon_count_matrix --run --keep /tmp/survom_demo_salmon
python scripts/smoke_workflows.py --case star_count_matrix --run --keep /tmp/survom_demo_star
python scripts/smoke_workflows.py --case qc_trim_strandedness --run --keep /tmp/survom_demo_qc_trim_strandedness
```

Compile only, without running Nextflow:

```bash
python scripts/smoke_workflows.py --case qc_trim_strandedness --keep /tmp/survom_compile_only
```

Run all demo workflow cases:

```bash
python scripts/smoke_workflows.py --run --keep /tmp/survom_all_demo_workflows
```

## 5. Run Student-Style Local Workflow

Use this when a student has their own input folder.

More command-line examples with custom parameters are in:

```text
docs/command_line_workflow_examples.md
```

From repo root:

```bash
cd /data/shared/vikash/mult-omics
```

Example input folder:

```text
/path/to/input_files/
  sample_R1.fastq.gz
  sample_R2.fastq.gz
```

Run trimming:

```bash
./run_workflow_local.sh \
  --workflow trim_only \
  --input /path/to/input_files \
  --output /tmp/my_trim_output
```

The runner creates:

```text
runs/<run_id>/
  state.json
  inputs/
  outputs/
  logs/
```

The command prints the `Run ID`.

## 6. Resume A Local Workflow

If a run failed or the server/app was closed, resume it:

```bash
./run_workflow_local.sh --resume <run_id>
```

Example:

```bash
./run_workflow_local.sh --resume trim_only_20260619_120000
```

Resume uses Nextflow `-resume`, so completed Nextflow tasks can be reused.

## 7. Input Manifest Format

You can provide FASTQ files directly, or provide a manifest.

Raw FASTQ manifest:

```csv
sample_id,fastq_1,fastq_2,single_end,strandedness
sample1,/path/sample1_R1.fastq.gz,/path/sample1_R2.fastq.gz,false,unknown
```

Trimmed FASTQ manifest:

```tsv
sample_id	fastq_1	fastq_2	single_end	strandedness
sample1	/path/sample1_R1.trimmed.fastq.gz	/path/sample1_R2.trimmed.fastq.gz	false	unknown
```

Useful manifest names:

```text
raw_manifest.csv
samplesheet.csv
trim_manifest.tsv
trimmed_fastq_manifest.tsv
```

If no manifest is found, the local runner tries to create one from FASTQ names like:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

## 8. Where To Check Results

For all-test script:

```text
.test-runs/<timestamp>/
```

For direct smoke workflows:

```text
/tmp/survom_demo_*/
```

For local student workflow runs:

```text
runs/<run_id>/
```

Useful files:

```text
state.json
logs/nextflow.stdout.log
logs/nextflow.stderr.log
logs/.nextflow.log
outputs/
```

## 9. Good Report Sentences

After full test passes:

```text
All automated code checks, workflow compile checks, and tiny real Nextflow smoke workflows passed.
```

After quick test passes:

```text
Quick code/unit/compile checks passed. Real Nextflow execution was not run.
```

If Docker permission fails:

```text
Code and compile checks passed, but real Nextflow execution could not run because Docker permission is missing.
```
