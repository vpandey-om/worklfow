# SurvOm Demo Dash App

This is a small teaching/demo web app for running the existing SurvOm modular Nextflow workflow in:

```text
/data/shared/vikash/mult-omics/survom-pipelines
```

It does not redesign the pipeline. It wraps the existing builder and job runner.

## Features

- Bootstrap Dash UI with a fixed left workflow panel.
- Analysis setup panel for student/tester, omics type, and workflow selection.
- Native Dash upload boxes for demo/test files.
- Upload multiple FASTQ files: `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`.
- Upload one metadata/sample sheet: `.csv`, `.tsv`, `.xlsx`.
- Dynamic workflow cards from `config/workflow_steps.yaml`.
- Run one step or the two-step workflow.
- SQLite job store survives app restart.
- Uploads and runs are not wiped on startup.
- Captures command, PID, status, logs, and output files.
- Uses one persistent browser session ID per user, so simultaneous uploads do not overwrite each other.
- Uses a lightweight tester workspace model so students only see their own uploads/jobs in the normal UI.
- Polls job status and outputs every 30 seconds while long Nextflow runs continue in the background.

## Directory Layout

```text
demo_dash_app/
  app.py
  server.py
  requirements.txt
  README.md
  config/workflow_steps.yaml
  uploads/
  runs/
  assets/style.css
  services/upload_service.py
  services/workflow_service.py
  services/job_store.py
  systemd/multomics-demo.service
```

## Setup

Use the same environment that runs the workflow tools:

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
source scripts/activate_survom_preprocess.sh
```

Install app dependencies into that environment:

```bash
cd /data/shared/vikash/mult-omics/demo_dash_app
python -m pip install -r requirements.txt
```

The current server did not have `dash`, `dash-bootstrap-components`, or `gunicorn` installed when this app was created.

## Run Locally

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
source scripts/activate_survom_preprocess.sh

cd /data/shared/vikash/mult-omics/demo_dash_app
python app.py
```

Open:

```text
http://SERVER:8056
```

## Run With Gunicorn

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
source scripts/activate_survom_preprocess.sh

cd /data/shared/vikash/mult-omics/demo_dash_app
gunicorn server:server --bind 0.0.0.0:8056 --workers 2 --timeout 300
```

## systemd Service

Edit `<USER>` first:

```bash
cp systemd/multomics-demo.service /tmp/multomics-demo.service
sed -i "s/<USER>/vikash/g" /tmp/multomics-demo.service
sudo cp /tmp/multomics-demo.service /etc/systemd/system/multomics-demo.service
sudo systemctl daemon-reload
sudo systemctl enable --now multomics-demo
```

Check status:

```bash
sudo systemctl status multomics-demo
```

## Upload Storage

Uploaded files are stored under:

```text
/data/shared/vikash/mult-omics/demo_dash_app/uploads/<omics_type>/<tester_id>/<session_id>/
```

For this demo, the browser and Dash callbacks use:

```text
survom_demo_session_id
survom_demo_tester_id
survom_demo_omics_type
```

These values are stored in browser local storage. Each browser/tester workspace gets its own upload folder, for example:

```text
uploads/bulk_rnaseq/tamalika/2f1f9f6a-5d30-4cc5-b78f-6c1c9f7c3b11/
```

This prevents two students uploading files at the same time from sharing one input directory in the demo UI.

## Run Storage

Runs are stored under:

```text
/data/shared/vikash/mult-omics/demo_dash_app/runs/<omics_type>/<tester_id>/<run_id>/
```

Each run contains:

```text
run_request.yaml
job.json
compile.stdout.log
compile.stderr.log
results/
```

Each button click creates a separate `run_id`, `job_id`, run folder, compiled workflow folder, and result folder. Two users can run the same step at the same time because their jobs use different directories.

The SQLite job table stores `session_id`, `tester_id`, `tester_label`, `omics_type`, and `workflow_id`. This is a demo/tester workspace model, not production authentication. For production, replace it with proper login/auth and user authorization.

The compiled Nextflow run is stored in the existing pipeline project:

```text
/data/shared/vikash/mult-omics/survom-pipelines/workflows/generated/<run_id>/
```

That folder includes:

```text
main.nf
params.yaml
nextflow.config
workflow_manifest.json
ui_run_summary.json
methods.md
run_command.sh
job_record.json
execution_record.json
logs/stdout.log
logs/stderr.log
```

## Workflow Configuration

Edit:

```text
config/workflow_steps.yaml
```

Current demo options:

```text
fastq_qc   -> rnaseq_03_raw_read_qc
trim_only  -> rnaseq_04_adapter_quality_trimming
qc_trim    -> rnaseq_03_raw_read_qc + rnaseq_04_adapter_quality_trimming
```

The omics type dropdown currently exposes:

```text
Bulk RNA-seq      active
Genomics          under development
Single-cell RNA-seq under development
Proteomics        under development
Metabolomics      under development
```

Only Bulk RNA-seq can launch workflows today. Other omics modules are visible so the UI/data model can grow without adding fake pipelines.

The app generates run requests and calls the existing pipeline builder:

```bash
python3 -m builder.api compile --request <run_request.yaml> --output <compiled_dir>
python3 -m builder.api run-job --job <job.json>
```

The web request starts the job in the background and returns immediately. The UI then polls the job table every 30 seconds and renders status, stdout/stderr tails, and output links when `execution_record.json` is available.

## Connecting Real Commands

The app currently uses the real SurvOm builder/job runner. If you later add more registry steps in `survom-pipelines`, add a new item to `config/workflow_steps.yaml` and map it to the corresponding `selected_steps`.

Do not hard-code new pipeline commands inside Dash callbacks. Keep workflow definitions in YAML and execution logic in `services/workflow_service.py`.

## Configurable Trimming Parameters

The `trim_only` and `qc_trim` demo workflows expose beginner-safe fastp settings in the Dash UI:

- minimum base quality, default `Q20`
- minimum read length after trimming, default `20`
- polyG trimming, default `auto`

Advanced options are collapsed by default and include tool selection, polyX trimming, explicit adapter sequences, adapter FASTA, fixed 5 prime trimming, and Cutadapt-specific error/overlap settings.

The default trimming tool is `fastp`. `cutadapt` is available for exact adapter handling, but the builder rejects Cutadapt runs early unless the required adapter sequence fields are provided. Dash still launches runs through the existing SurvOm builder and Nextflow path; it does not call fastp or Cutadapt directly.

## Production FASTQ Uploads

The current app uses native Dash upload boxes so the teaching demo can run immediately. For production 5-50 GB FASTQ uploads, use the tus plan in:

```text
docs/fastq_upload_strategy.md
```

The recommended shape is Uppy Dashboard + `@uppy/tus` + Golden Retriever in the browser, with `tusd` running as a standalone upload service. Keep completed files registered into the same `uploads/<session_id>/` layout so the existing workflow execution layer does not need to change.
