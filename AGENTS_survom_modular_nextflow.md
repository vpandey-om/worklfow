# AGENTS.md — SurvOm Modular Nextflow Workflow Builder

This file defines how Codex/AI coding agents must work on the SurvOm modular multiomics workflow system.

The goal is to let users build dynamic workflows step by step, where each process can run individually or as part of a selected chain, while preserving all existing working modules.

---

## 1. Prime Directive

Build incrementally. Do not break working steps.

Every change must preserve:

- existing step IDs
- existing registry fields
- existing Nextflow module process names
- existing input/output contracts
- existing generated workflow behavior
- existing tests and examples
- existing local and AWS execution profiles

When adding a new step, add it as a new registry entry and a new reusable module. Do not rewrite unrelated modules.

---

## 2. Non-Negotiable Rules

### Do not delete or rename working contracts

Never rename these without a migration plan:

- `step_id`
- `process_name`
- emitted Nextflow channel names
- output names in the registry
- manifest column names
- public API response fields
- generated params keys

Bad:

```text
rename rnaseq_03_raw_read_qc to raw_qc
rename FASTQC to RUN_FASTQC
rename trimmed_fastq_manifest to cleaned_fastqs
```

Good:

```text
keep rnaseq_03_raw_read_qc
add aliases only if needed
add a schema_version or deprecated flag when replacing behavior
```

### Do not create one monolithic pipeline

Never hard-code the entire RNA-seq flow as one fixed script.

Good architecture:

```text
registry -> planner -> compiler -> generated Nextflow workflow -> reusable modules
```

Bad architecture:

```text
one giant rnaseq_pipeline.nf containing all 20 steps directly
```

### Do not mix UI logic into Nextflow modules

Nextflow modules should only do scientific processing.

UI/chat guidance belongs in:

```text
registry/steps.yaml
builder/planner.py
frontend step cards
```

### Do not make Nextflow modules depend on RNA-seq when they are COMMON

COMMON modules must be reusable by genomics, ATAC-seq, ChIP-seq, methylation, and other omics workflows.

For example:

```text
modules/common/fastqc/main.nf
modules/common/fastp/main.nf
modules/common/multiqc/main.nf
```

These modules must not contain RNA-seq-specific assumptions like condition columns, DESeq2 design, or transcriptome references.

---

## 3. Repository Ownership Map

Use this ownership model when changing files.

| Area | Owner agent | Change policy |
|---|---|---|
| `registry/steps.yaml` | Registry Curator | Additive changes preferred. Do not rename working IDs. |
| `modules/common/*` | Nextflow Module Engineer | Keep modules generic and reusable. |
| `modules/rnaseq/*` | RNA-seq Module Engineer | RNA-seq-specific logic only. |
| `builder/planner.py` | Workflow Planner Engineer | Validate step graph and suggest next steps. |
| `builder/compiler.py` | Workflow Compiler Engineer | Generate Nextflow from selected steps. |
| `builder/runner.py` | Execution Engineer | Local/AWS execution wrapper only. |
| `schemas/*` | Schema Engineer | Version schemas; avoid breaking old requests. |
| `tests/*` | QA Agent | Add tests for every new step and compiler case. |
| `docs/*` | Documentation Agent | Keep user/developer docs synced with registry. |

---

## 4. Recommended Agents / Skills

### Agent 1 — Registry Curator

Purpose: maintain the workflow step registry.

Responsibilities:

- add new step definitions
- validate required inputs, outputs, and params
- define `next_steps.recommended` and `next_steps.allowed`
- mark modules as `COMMON`, `RNA_SEQ_SPECIFIC`, `GENOMICS_COMMON`, `STATISTICS_COMMON`, or `REPORTING_COMMON`
- keep parameter exposure as `user`, `expert`, or `internal`

Must not:

- change existing `step_id` values
- remove working outputs
- add tool flags not implemented in the Nextflow module
- expose advanced/internal parameters to beginner users

Completion checklist:

```text
[ ] step_id is stable and unique
[ ] category is correct
[ ] omics reuse list is correct
[ ] required inputs are declared
[ ] outputs are declared with stable names
[ ] user/expert/internal params are separated
[ ] next_steps are defined
[ ] tests updated
```

---

### Agent 2 — Nextflow Module Engineer

Purpose: implement reusable DSL2 modules.

Responsibilities:

- create one module per tool/process
- keep modules small and testable
- use stable input and output channels
- emit `versions.yml`
- publish expected output files
- support both single-end and paired-end FASTQ when relevant
- use containers from config, not hard-coded runtime logic

Must not:

- embed UI text in modules
- hard-code project-specific paths
- use absolute paths like `/data/shared/vikash/...`
- change emitted output names without updating registry and tests
- combine multiple scientific steps into one process unless the step is intentionally composite

Preferred process style:

```nextflow
process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"),  emit: zip
    path "versions.yml", emit: versions

    script:
    """
    fastqc --threads ${task.cpus} --outdir . --extract ${reads}
    """
}
```

Completion checklist:

```text
[ ] process name matches registry
[ ] input tuple contract is documented
[ ] output emits match registry
[ ] versions.yml emitted
[ ] container configured in nextflow.config
[ ] works for local profile
[ ] works for docker profile
[ ] compatible with AWS Batch profile
```

---

### Agent 3 — Workflow Planner Engineer

Purpose: decide what steps can be added next.

Responsibilities:

- read `registry/steps.yaml`
- validate selected steps
- detect missing dependencies
- allow individual step execution if required inputs are provided
- suggest next steps after each selected step
- prevent invalid graphs
- return user-friendly reasons for recommendations

Must not:

- hard-code RNA-seq-only step order in Python
- assume all workflows are linear forever
- silently skip missing required inputs
- auto-add destructive steps without user confirmation

Core planner functions:

```python
def get_available_first_steps(omics: str) -> list[Step]:
    pass


def suggest_next_steps(selected_steps: list[str], omics: str) -> list[StepSuggestion]:
    pass


def validate_selected_steps(selected_steps: list[str], omics: str, provided_inputs: dict) -> ValidationResult:
    pass


def resolve_effective_params(selected_steps: list[str], user_params: dict) -> dict:
    pass
```

Completion checklist:

```text
[ ] raw QC suggests trimming
[ ] trimming suggests post-trim QC, Salmon, or STAR
[ ] trimming can run alone if FASTQ manifest is supplied
[ ] missing FASTQ input gives clear error
[ ] invalid next step gives clear error
```

---

### Agent 4 — Workflow Compiler Engineer

Purpose: generate runnable Nextflow workflows from selected steps.

Responsibilities:

- convert selected steps into `main.nf`
- include only required modules
- wire channels between modules
- create `params.yaml`
- create `nextflow.config`
- create `workflow_manifest.json`
- keep generated workflows reproducible

Must not:

- manually write one workflow per possible step combination
- duplicate module code into generated workflows
- alter source modules during compilation
- invent parameters not in registry

Compilation strategy:

```text
run_request.yaml
  -> validate selected steps
  -> resolve parameters
  -> identify module includes
  -> generate workflow channel wiring
  -> write generated workflow folder
```

Generated folder example:

```text
workflows/generated/run001/
├── main.nf
├── nextflow.config
├── params.yaml
└── workflow_manifest.json
```

Completion checklist:

```text
[ ] raw QC-only workflow generated
[ ] trimming-only workflow generated
[ ] raw QC + trimming workflow generated
[ ] params.yaml contains defaults + overrides
[ ] source modules are not modified
[ ] generated workflow passes nextflow syntax check
```

---

### Agent 5 — Execution Engineer

Purpose: run generated workflows locally, on HPC, or on AWS.

Responsibilities:

- implement run wrapper
- support `local`, `local_docker`, `singularity`, `aws_batch`
- record run status
- capture logs
- collect outputs
- preserve provenance

Must not:

- create separate pipeline code for AWS
- embed AWS credentials in files
- hard-code bucket names
- write outputs outside declared `outdir`

Execution command examples:

```bash
nextflow run workflows/generated/run001/main.nf \
  -profile local,docker \
  -params-file workflows/generated/run001/params.yaml \
  -resume
```

```bash
nextflow run workflows/generated/run001/main.nf \
  -profile aws \
  -params-file workflows/generated/run001/params.yaml \
  -work-dir s3://survom-work/work \
  -resume
```

Completion checklist:

```text
[ ] local execution works
[ ] docker execution works
[ ] AWS profile validates aws_queue, aws_region, aws_workdir
[ ] run logs are captured
[ ] run manifest is written
[ ] outputs can be listed for UI/chat
```

---

### Agent 6 — QA / Regression Agent

Purpose: prevent breaking existing working steps.

Responsibilities:

- write tests before or alongside changes
- test registry validation
- test planner recommendations
- test generated Nextflow snapshots
- test backward compatibility
- test each module independently

Must not:

- delete tests to make builds pass
- update snapshots without explaining why
- skip tests for changed registry fields

Required tests for every new step:

```text
[ ] registry schema validation
[ ] step ID uniqueness
[ ] required inputs declared
[ ] outputs declared
[ ] params exposure declared
[ ] next step suggestions valid
[ ] generated workflow includes correct module
[ ] generated params include defaults
[ ] individual step run request valid
[ ] chained workflow request valid
```

---

### Agent 7 — Documentation Agent

Purpose: keep user-facing and developer-facing docs accurate.

Responsibilities:

- update tutorial pages
- update parameter tables
- update examples
- update command snippets
- explain beginner vs expert settings
- document outputs for UI cards

Must not:

- document parameters that are not implemented
- expose internal implementation details to beginner users
- claim AWS support for a module that was not tested with container/profile

---

## 5. Safe Change Workflow for Codex

Every Codex task should follow this order:

```text
1. Inspect current registry and modules.
2. Identify the smallest change needed.
3. Add or modify registry entry.
4. Add or modify one module only if needed.
5. Update planner/compiler only if the current logic cannot support the new step.
6. Add tests.
7. Run tests.
8. Show changed files and explain compatibility.
```

Codex must answer these before editing:

```text
What existing step IDs are affected?
What output contracts are affected?
Can existing workflows still run?
Is this change additive or breaking?
Which tests prove backward compatibility?
```

---

## 6. Compatibility Policy

### Additive changes are safe

Examples:

```text
add rnaseq_05_post_trim_quality_control
add modules/common/multiqc/main.nf
add one new optional parameter with default value
add one new allowed next step
```

### Breaking changes require migration

Examples:

```text
rename output trimmed_fastq_manifest
change FASTP output pattern
remove parameter quality_threshold
change sample sheet column fastq_1 to r1
```

For breaking changes, create:

```text
registry/migrations/YYYYMMDD_description.md
```

Include:

```text
old behavior
new behavior
why change is needed
migration steps
affected workflows
tests updated
```

---

## 7. Parameter Exposure Policy

Use exactly three exposure levels.

### user

Visible in guided UI and chat.

Examples:

```text
quality_threshold
minimum_read_length
trim_poly_g
run_trimming
```

### expert

Hidden under advanced settings.

Examples:

```text
adapter_sequence_r1
adapter_sequence_r2
cutadapt_error_rate
custom_contaminants
```

### internal

Controlled by platform.

Examples:

```text
threads
container_image
publish_dir
work_dir
output_prefix
```

Do not expose internal parameters in beginner UI.

---

## 8. Output Contract Policy

Every step must produce machine-readable outputs for UI/chat.

Minimum output contract:

```yaml
outputs:
  - name: primary_output
    type: file_or_manifest
    pattern: path/pattern
    ui_visible: true
  - name: summary_json
    type: json
    pattern: path/summary.json
    ui_visible: true
  - name: versions
    type: yaml
    pattern: versions.yml
    ui_visible: false
```

UI cards should never parse random logs when a JSON/TSV summary can be produced.

---

## 9. Step Addition Template

When adding a new step, Codex should use this checklist.

```text
New step ID:
Category:
Omics reuse:
Default tool:
Fallback tool:
Module path:
Process name:
Required input contract:
Output contract:
User parameters:
Expert parameters:
Internal parameters:
Recommended next steps:
Allowed next steps:
Can run individually:
Can be skipped:
Tests added:
Docs updated:
```

---

## 10. Codex Task Prompt Template

Use this prompt for any new coding task.

```text
You are working in the SurvOm modular Nextflow workflow builder.

Follow AGENTS.md strictly.

Task:
<describe exact task>

Constraints:
- Make the smallest safe change.
- Do not rename existing step IDs, process names, output names, or params.
- Do not remove working code.
- Prefer additive changes.
- Keep COMMON modules reusable across omics.
- If a breaking change is unavoidable, create a migration note and explain it.
- Add or update tests.
- Show changed files and backward-compatibility impact.

Before coding, inspect:
- registry/steps.yaml
- builder/planner.py
- builder/compiler.py
- modules/common
- tests

After coding, run relevant tests and summarize:
- what changed
- what did not change
- how existing workflows remain safe
- how to run the new or changed step individually
- how to run it in a chain
```

---

## 11. Example Task: Add Post-Trim QC Safely

Prompt:

```text
Add rnaseq_05_post_trim_quality_control as a new registry step using the existing common/fastqc module.

Do not modify rnaseq_03_raw_read_qc or rnaseq_04_adapter_quality_trimming behavior.

The new step should accept trimmed_fastq_manifest or trimmed FASTQ files, run FastQC, emit post_trim_fastqc_html, post_trim_fastqc_zip, and post_trim_qc_summary, and recommend Salmon or STAR as next steps.

Update planner tests and compiler tests so this step can run individually and after trimming.
```

Expected safe change:

```text
registry/steps.yaml                       modified additively
builder/planner.py                        maybe unchanged
builder/compiler.py                       maybe unchanged if generic
modules/common/fastqc/main.nf             unchanged
examples/rnaseq_post_trim_qc.run.yaml     added
tests/test_registry.py                    updated
tests/test_planner.py                     updated
tests/test_compiler.py                    updated
```

---

## 12. Example Task: Add Genomics FASTQ QC Reusing Existing Module

Prompt:

```text
Add a genomics FASTQ QC step that reuses modules/common/fastqc/main.nf.

Do not duplicate the FastQC module.
Do not change RNA-seq step IDs or outputs.
Create a new genomics step ID in registry/steps.yaml and mark it as GENOMICS_COMMON or COMMON as appropriate.
Add planner tests showing the module is reusable across rnaseq and genomics.
```

Expected safe change:

```text
registry/steps.yaml       modified additively
modules/common/fastqc     unchanged
tests                     updated
```

---

## 13. Red Flags

Stop and ask for review if a change requires:

```text
renaming a step_id
renaming a process_name
changing sample sheet schema
changing output manifest schema
changing generated workflow structure for all steps
removing an existing parameter
changing AWS profile behavior
moving modules between directories
```

---

## 14. Definition of Done

A task is done only when:

```text
[ ] existing examples still compile
[ ] existing tests pass
[ ] new tests are added for new behavior
[ ] step can run individually if intended
[ ] step can run in a chain if intended
[ ] registry and module contracts match
[ ] documentation or examples are updated
[ ] backward compatibility is explained
```
