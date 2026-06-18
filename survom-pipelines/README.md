# SurvOm Modular Nextflow Workflow Builder

This directory contains a registry-driven Nextflow DSL2 workflow builder for SurvOm. It is organized so a UI or chat interface can build workflows step by step without hard-coding one monolithic RNA-seq pipeline.

## Architecture

```text
registry/steps.yaml -> builder/planner.py -> builder/compiler.py -> generated Nextflow -> reusable modules
```

The first implemented route supports RNA-seq FASTQ preprocessing:

```text
rnaseq_03_raw_read_qc
rnaseq_04_adapter_quality_trimming
```

Both steps are implemented as reusable COMMON modules, so they can also be reused by genomics, ATAC-seq, ChIP-seq, and methylation workflows.

## Layout

```text
registry/                 machine-readable step metadata and schemas
modules/common/           reusable omics-agnostic DSL2 modules
modules/rnaseq/           RNA-seq-specific module placeholders
workflows/templates/      Jinja2 templates for generated Nextflow files
workflows/generated/      compiled workflows
builder/                  Python registry, planner, compiler, runner, API/CLI
examples/                 run request YAML files
tests/                    registry, planner, and compiler tests
```

## Compile Examples

From this directory:

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
python3 -m builder.api compile --request examples/rnaseq_qc_only.run.yaml --output workflows/generated/run_qc_only
python3 -m builder.api compile --request examples/rnaseq_qc_trim.run.yaml --output workflows/generated/run_qc_trim
python3 -m builder.api compile --request examples/rnaseq_trim_only.run.yaml --output workflows/generated/run_trim_only
```

Run compiled Nextflow directly:

```bash
nextflow run workflows/generated/run_qc_trim/main.nf \
  -profile local,docker \
  -params-file workflows/generated/run_qc_trim/params.yaml \
  -resume
```

AWS uses the same generated pipeline:

```bash
nextflow run workflows/generated/run_qc_trim/main.nf \
  -profile aws \
  -params-file workflows/generated/run_qc_trim/params.yaml \
  --aws_queue survom-rnaseq-queue \
  --aws_region eu-north-1 \
  --aws_workdir s3://survom-rnaseq/work/run_qc_trim \
  -resume
```

## Planner Behavior

The planner can:

- list first steps for an omics type,
- suggest recommended and allowed next steps,
- validate selected step chains,
- resolve default parameters plus user overrides.

Trimming can run individually when a valid FASTQ manifest is supplied, even though it normally follows raw QC.

## API

`builder/api.py` exposes a `create_app()` function for FastAPI:

```bash
uvicorn builder.api:create_app --factory --reload
```

Endpoints:

```text
GET  /api/omics
GET  /api/steps?omics=rnaseq
GET  /api/steps/{step_id}
POST /api/workflows/plan
POST /api/workflows/compile
POST /api/workflows/run
GET  /api/runs/{run_id}
GET  /api/runs/{run_id}/outputs
```

## Testing

The project uses `pytest` tests:

```bash
python3 -m pytest tests
```

If `pytest` is not installed:

```bash
python3 -m pip install pytest
```

## Notes

The module code is intentionally separated from workflow orchestration. The compiler generates `main.nf`, `nextflow.config`, `params.yaml`, and `workflow_manifest.json` for each selected workflow.
