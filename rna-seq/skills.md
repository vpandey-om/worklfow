# SurvOm RNA-seq Agent Notes

Use these notes when working on the current RNA-seq demo and Nextflow builder.

## Keep Stable

- Do not rename existing step IDs, process names, output names, or manifest columns.
- Keep `fastp` as the default trimming tool.
- Keep `cutadapt` as an advanced fallback for exact adapter handling.
- Do not rewrite working FastQC, fastp, or Cutadapt modules unless the user asks for a scientific behavior change.

## Known-Good Runtime Images

Use these verified container tags:

```text
FASTQC   quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
FASTP    quay.io/biocontainers/fastp:0.23.4--h5f740d0_0
CUTADAPT quay.io/biocontainers/cutadapt:5.1--py312h0fa9677_0
```

Do not reintroduce these broken tags:

```text
biocontainers/fastqc:v0.12.1_cv4
quay.io/biocontainers/cutadapt:5.1--py313h1a76870_0
```

Before changing any image, verify it with `docker pull`.

## Trimming Rules

- For `trimming_tool: fastp`, generated config should include `withName:FASTP` only.
- For `trimming_tool: cutadapt`, generated config should include `withName:CUTADAPT` only.
- Do not emit stale process selectors for a process that is not included in `main.nf`.
- If local Cutadapt is unavailable, use container execution such as `local_docker`.

## Atomic Workflow Rules

- Treat workflows as chains of atomic step IDs from `demo_dash_app/config/workflow_steps.yaml`.
- Keep step contracts in `survom-pipelines/registry/steps.yaml`.
- Do not hard-code full RNA-seq chains inside Dash callbacks or one Nextflow file.
- Post-trim QC must run after trimming and consume trimmed FASTQ.
- Strandedness inference runs after post-trim QC. Salmon `-l A` is the default; STAR/RSeQC is fallback or validation.
- Salmon count matrix workflow route is Salmon; STAR count matrix workflow route is STAR.
- Reference validation must run before Salmon quantification or STAR/featureCounts counting.
- Do not emit stale Nextflow config selectors for processes absent from generated `main.nf`.

## Minimal Regression Checks

After any code or workflow change, run the repository test wrapper from repo root:

```bash
cd /data/shared/vikash/mult-omics
./scripts/run_all_tests_and_workflows.sh --quick
```

Before claiming the workflow works end-to-end, run the full wrapper with Docker/Nextflow access:

```bash
cd /data/shared/vikash/mult-omics
./scripts/run_all_tests_and_workflows.sh
```

Meaning of the two modes:

```text
--quick = code/unit/compile checks only
default = code/unit/compile checks + tiny real FASTQ/reference Nextflow smoke runs
```

If Docker is unavailable, clearly say:

```text
Code and compile checks passed, but real Nextflow execution could not run because Docker permission is missing.
```

Also compile/check both trimming paths when touching trimming:

```text
qc_trim with fastp
qc_trim with cutadapt
```

Expected outputs for Cutadapt:

```text
*.trimmed.fastq.gz
*.cutadapt.json
*.cutadapt.log
trimmed_fastq_manifest.tsv
```

## Tutorial And CLI Contract

Whenever adding or changing a workflow, step, parameter, input contract, output contract, or resume behavior:

- update tests or smoke cases first, preferably in `survom-pipelines/scripts/smoke_workflows.py` and/or `survom-pipelines/tests/`
- update `docs/how_to_run_tests_and_demo_workflows.md`
- update `docs/command_line_workflow_examples.md` with at least one real command-line example
- keep `README.md` pointing to the current tutorial
- keep `run_workflow_local.sh` working as a fallback when Dash/API is unavailable

Students and instructors must be able to test a workflow from the command line with:

```bash
./run_workflow_local.sh \
  --workflow <workflow_name> \
  --input <input_dir> \
  --output <output_dir>
```

Parameter overrides must be possible with repeatable `--param key=value` arguments:

```bash
./run_workflow_local.sh \
  --workflow trim_only \
  --input /path/to/input_files \
  --output /tmp/survom_trim_test \
  --param quality_threshold=25 \
  --param minimum_read_length=30
```

Reference/index examples should use params too:

```bash
./run_workflow_local.sh \
  --workflow salmon_quant_only \
  --input /path/to/trimmed_inputs \
  --output /tmp/survom_salmon_test \
  --param salmon_index=/path/to/salmon_index
```

Resume must stay available:

```bash
./run_workflow_local.sh --resume <run_id>
```

Do not add a workflow that only works through Dash. Every implemented workflow should have a command-line path for debugging and teaching.
