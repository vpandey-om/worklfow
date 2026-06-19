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

After changes, compile/check both paths:

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
