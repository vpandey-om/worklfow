# Student Project Result Tracking

This note explains where Dash stores student uploads, projects, jobs, and workflow outputs.

Start from repo root:

```bash
cd /data/shared/vikash/mult-omics
```

## 1. Dash Storage Layout

Uploads:

```text
demo_dash_app/uploads/<omics_type>/<tester_id>/<project_id>/<session_id>/
  fastq/
  metadata/
  reference/
  other/
  uploads_manifest.json
```

Example:

```text
demo_dash_app/uploads/bulk_rnaseq/vikash/testdata/<session_id>/fastq/
```

Runs:

```text
demo_dash_app/runs/<omics_type>/<tester_id>/<project_id>/<run_id>/
  run_request.yaml
  job.json
  compile.stdout.log
  compile.stderr.log
  results/
```

Example:

```text
demo_dash_app/runs/bulk_rnaseq/vikash/testdata/qc_trim_strandedness_20260619_120000_abcd1234/
```

Compiled Nextflow workflows:

```text
survom-pipelines/workflows/generated/<run_id>/
  main.nf
  params.yaml
  workflow_plan.json
  execution_record.json
  logs/
  work/
```

Database:

```text
demo_dash_app/runs/jobs.sqlite
```

This stores:

```text
projects
upload_records
jobs
```

## 2. List Student Uploads

Example for Vikash/testdata:

```bash
find demo_dash_app/uploads/bulk_rnaseq/vikash/testdata -maxdepth 4 -type f | sort
```

List only FASTQ:

```bash
find demo_dash_app/uploads/bulk_rnaseq/vikash/testdata -path '*/fastq/*' -type f | sort
```

List metadata:

```bash
find demo_dash_app/uploads/bulk_rnaseq/vikash/testdata -path '*/metadata/*' -type f | sort
```

List reference files:

```bash
find demo_dash_app/uploads/bulk_rnaseq/vikash/testdata -path '*/reference/*' -type f | sort
```

## 3. List Student Runs

```bash
find demo_dash_app/runs/bulk_rnaseq/vikash/testdata -maxdepth 2 -type f | sort
```

Find job summaries:

```bash
find demo_dash_app/runs/bulk_rnaseq/vikash/testdata -name job.json -print
```

Find run requests:

```bash
find demo_dash_app/runs/bulk_rnaseq/vikash/testdata -name run_request.yaml -print
```

Inspect one job:

```bash
cat demo_dash_app/runs/bulk_rnaseq/vikash/testdata/<run_id>/job.json
```

Inspect selected files and parameters:

```bash
cat demo_dash_app/runs/bulk_rnaseq/vikash/testdata/<run_id>/run_request.yaml
```

## 4. Query SQLite

Open the DB:

```bash
sqlite3 demo_dash_app/runs/jobs.sqlite
```

Inside SQLite:

```sql
.tables
select project_id, project_label, tester_id, omics_type, updated_at from projects order by updated_at desc;
select tester_id, omics_type, project_id, category, filename, stored_path from upload_records order by created_at desc limit 20;
select job_id, status, tester_id, omics_type, project_id, workflow_id, created_at from jobs order by created_at desc limit 20;
```

One-line command:

```bash
sqlite3 demo_dash_app/runs/jobs.sqlite \
  "select job_id,status,tester_id,project_id,workflow_id,created_at from jobs order by created_at desc limit 20;"
```

## 5. Check A Failed Run

Find newest run:

```bash
find demo_dash_app/runs/bulk_rnaseq/vikash/testdata -mindepth 1 -maxdepth 1 -type d | sort | tail -1
```

Check compile logs:

```bash
tail -n 100 demo_dash_app/runs/bulk_rnaseq/vikash/testdata/<run_id>/compile.stderr.log
```

Check Nextflow logs:

```bash
tail -n 100 survom-pipelines/workflows/generated/<run_id>/logs/stderr.log
tail -n 100 survom-pipelines/workflows/generated/<run_id>/.nextflow.log
```

Check collected execution record:

```bash
cat survom-pipelines/workflows/generated/<run_id>/execution_record.json
```

## 6. Backend Command-Line Test Instead Of Dash

Use this when Dash/API has an issue but you want to test the workflow itself.

Raw FASTQ input folder:

```bash
./run_workflow_local.sh \
  --workflow qc_trim_strandedness \
  --input /path/to/fastq_folder \
  --output /tmp/survom_cli_qc_trim_strandedness
```

With params:

```bash
./run_workflow_local.sh \
  --workflow trim_only \
  --input /path/to/fastq_folder \
  --output /tmp/survom_cli_trim_q25 \
  --param quality_threshold=25 \
  --param minimum_read_length=30
```

With a Salmon index:

```bash
./run_workflow_local.sh \
  --workflow salmon_quant_only \
  --input /path/to/trimmed_manifest_folder \
  --output /tmp/survom_cli_salmon \
  --param salmon_index=/path/to/salmon_index
```

Resume:

```bash
./run_workflow_local.sh --resume <run_id>
```

## 7. Practical Rule For Testers

Use one project per class exercise:

```text
tester: vikash
omics: bulk_rnaseq
project: testdata
```

Then all uploads and outputs stay grouped under:

```text
demo_dash_app/uploads/bulk_rnaseq/vikash/testdata/
demo_dash_app/runs/bulk_rnaseq/vikash/testdata/
```

This makes it easy to clean, inspect, or archive one student's project later.
