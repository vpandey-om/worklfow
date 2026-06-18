# SurvOm Dynamic Modular Nextflow Pipeline Builder — Implementation Prompt

You are implementing a production-grade modular multiomics workflow builder for SurvOm.

## Goal

Build a registry-driven system where a biologist can create a workflow step-by-step from a UI or chat interface. Each process/module must be able to run:

1. Individually as a single step.
2. In a selected chain of steps.
3. As a full predefined workflow such as RNA-seq 20-step analysis.
4. Reused across other omics workflows such as genomics, ATAC-seq, ChIP-seq, scRNA-seq, methylation, microbiome, and proteomics where appropriate.

The first implementation target is RNA-seq preprocessing using:

- `rnaseq_03_raw_read_qc`
- `rnaseq_04_adapter_quality_trimming`

But the architecture must scale to all RNA-seq modules and other omics pipelines.

## Core Architecture Rule

Do not hard-code one monolithic RNA-seq pipeline.

Use three layers:

1. **Step registry layer** — YAML/JSON metadata describing each module, inputs, outputs, parameters, dependencies, and compatible next steps.
2. **Module implementation layer** — reusable Nextflow DSL2 modules, preferably grouped by domain: `common`, `rnaseq`, `genomics`, `statistics`, `reporting`.
3. **Workflow builder layer** — Python or TypeScript code that validates the selected steps and generates a valid Nextflow `main.nf`, `nextflow.config`, and `params.yaml` for that selected workflow.

Nextflow modules should remain static and reusable. The dynamic behavior should happen in the workflow builder, which compiles/generates the selected workflow before running Nextflow.

## Required Repository Structure

Create this structure:

```text
survom-pipelines/
├── registry/
│   ├── steps.yaml
│   ├── routes.yaml
│   └── schemas/
│       ├── step.schema.json
│       ├── samplesheet.schema.json
│       └── run_request.schema.json
├── modules/
│   ├── common/
│   │   ├── fastqc/main.nf
│   │   ├── falco/main.nf
│   │   ├── fastp/main.nf
│   │   ├── cutadapt/main.nf
│   │   ├── multiqc/main.nf
│   │   └── validate_samplesheet/main.nf
│   ├── rnaseq/
│   │   ├── salmon/main.nf
│   │   ├── star/main.nf
│   │   ├── featurecounts/main.nf
│   │   ├── tximport/main.nf
│   │   └── deseq2/main.nf
│   ├── genomics/
│   │   └── README.md
│   └── reporting/
│       └── render_report/main.nf
├── workflows/
│   ├── templates/
│   │   ├── main.nf.j2
│   │   ├── nextflow.config.j2
│   │   └── params.yaml.j2
│   └── generated/
│       └── .gitkeep
├── builder/
│   ├── __init__.py
│   ├── models.py
│   ├── registry.py
│   ├── planner.py
│   ├── compiler.py
│   ├── runner.py
│   └── api.py
├── examples/
│   ├── rnaseq_qc_only.run.yaml
│   ├── rnaseq_qc_trim.run.yaml
│   └── rnaseq_full_preprocess.run.yaml
├── tests/
│   ├── test_registry.py
│   ├── test_planner.py
│   ├── test_compiler.py
│   └── test_nextflow_generation.py
└── README.md
```

## Step Registry Requirements

Implement `registry/steps.yaml` as the source of truth for UI, chat, validation, and workflow generation.

Each step must define:

- `step_id`
- `step_name`
- `category`
- `omics` list
- `purpose`
- `module_path`
- `process_name`
- `default_tool`
- `fallback_tool`
- `depends_on`
- `optional_depends_on`
- `compatible_previous_outputs`
- `required_inputs`
- `outputs`
- `parameters`
- `next_steps`
- `can_run_individually`
- `can_be_skipped`
- `skip_conditions`
- `resources`
- `containers`
- `ui`
- `chat_guidance`

Use this first registry content:

```yaml
version: "1.0"
steps:
  rnaseq_03_raw_read_qc:
    step_id: rnaseq_03_raw_read_qc
    step_name: Raw Read Quality Control
    category: COMMON
    omics: [rnaseq, genomics, atacseq, chipseq, methylation]
    purpose: Assess sequencing quality before modifying reads.
    module_path: modules/common/fastqc/main.nf
    process_name: FASTQC
    default_tool: fastqc
    fallback_tool: falco
    can_run_individually: true
    can_be_skipped: false
    depends_on:
      - rnaseq_02_fastq_staging_integrity
    optional_depends_on: []
    compatible_previous_outputs:
      - staged_fastq_manifest
      - fastq_manifest
      - raw_fastq_files
    required_inputs:
      - name: fastq_manifest
        type: manifest
        formats: [csv, tsv, json]
        required: true
      - name: fastq_1
        type: file
        formats: [fastq, fastq.gz, fq, fq.gz]
        required: true
      - name: fastq_2
        type: file
        formats: [fastq, fastq.gz, fq, fq.gz]
        required_when: library_layout == paired
    outputs:
      - name: fastqc_html
        type: file
        pattern: "03_raw_qc/fastqc/*_fastqc.html"
        ui_visible: true
      - name: fastqc_zip
        type: file
        pattern: "03_raw_qc/fastqc/*_fastqc.zip"
        ui_visible: true
      - name: raw_qc_summary
        type: json
        pattern: "03_raw_qc/raw_qc_summary.json"
        ui_visible: true
      - name: raw_qc_status
        type: tsv
        pattern: "03_raw_qc/raw_qc_status.tsv"
        ui_visible: true
    parameters:
      run_raw_qc:
        type: boolean
        default: true
        exposure: user
        ui_label: Run raw read QC
      qc_tool:
        type: enum
        choices: [fastqc, falco]
        default: fastqc
        exposure: expert
        ui_label: QC tool
      threads:
        type: integer
        default: 4
        min: 1
        max: 32
        exposure: internal
      custom_adapters:
        type: file
        default: null
        exposure: expert
      custom_contaminants:
        type: file
        default: null
        exposure: expert
    next_steps:
      recommended:
        - rnaseq_04_adapter_quality_trimming
      allowed:
        - rnaseq_04_adapter_quality_trimming
        - rnaseq_07_salmon_quantification
        - rnaseq_09_star_alignment
    resources:
      cpus: 4
      memory: 8.GB
      time: 2.h
    containers:
      fastqc: biocontainers/fastqc:v0.12.1_cv4
      falco: quay.io/biocontainers/falco:1.3.2--h077b44d_0
    ui:
      card_title: Check raw FASTQ quality
      card_description: Detect adapter contamination, low-quality cycles, GC bias, duplication, and overrepresented sequences.
      mode: automated
      expose_in_guided: true
    chat_guidance:
      beginner: "I will check the raw FASTQ quality first using FastQC."
      expert: "Use Falco only for large cohorts or speed-sensitive runs."

  rnaseq_04_adapter_quality_trimming:
    step_id: rnaseq_04_adapter_quality_trimming
    step_name: Adapter and Quality Trimming
    category: COMMON
    omics: [rnaseq, genomics, atacseq, chipseq, methylation]
    purpose: Remove adapters, poor-quality read tails, and platform-specific artifacts before mapping or quantification.
    module_path: modules/common/fastp/main.nf
    process_name: FASTP
    default_tool: fastp
    fallback_tool: cutadapt
    can_run_individually: true
    can_be_skipped: true
    skip_conditions:
      - raw QC shows no adapter contamination and high read quality
      - expert user explicitly selects no trimming
    depends_on:
      - rnaseq_03_raw_read_qc
    optional_depends_on:
      - rnaseq_02_fastq_staging_integrity
    compatible_previous_outputs:
      - staged_fastq_manifest
      - raw_qc_summary
      - raw_fastq_files
    required_inputs:
      - name: fastq_manifest
        type: manifest
        formats: [csv, tsv, json]
        required: true
      - name: fastq_1
        type: file
        formats: [fastq, fastq.gz, fq, fq.gz]
        required: true
      - name: fastq_2
        type: file
        formats: [fastq, fastq.gz, fq, fq.gz]
        required_when: library_layout == paired
    outputs:
      - name: trimmed_fastq_manifest
        type: tsv
        pattern: "04_trimmed/trimmed_fastq_manifest.tsv"
        ui_visible: true
      - name: trimmed_fastq
        type: file
        pattern: "04_trimmed/fastq/*.trimmed.fastq.gz"
        ui_visible: true
      - name: trimming_html
        type: file
        pattern: "04_trimmed/reports/*.fastp.html"
        ui_visible: true
      - name: trimming_json
        type: json
        pattern: "04_trimmed/reports/*.fastp.json"
        ui_visible: true
      - name: trimming_summary
        type: json
        pattern: "04_trimmed/trimming_summary.json"
        ui_visible: true
    parameters:
      run_trimming:
        type: boolean
        default: true
        exposure: user
        ui_label: Run adapter and quality trimming
      trimming_tool:
        type: enum
        choices: [fastp, cutadapt]
        default: fastp
        exposure: expert
        ui_label: Trimming tool
      quality_threshold:
        type: integer
        default: 20
        min: 5
        max: 35
        exposure: user
        ui_label: Minimum base quality
      minimum_read_length:
        type: integer
        default: 20
        min: 10
        max: 100
        exposure: user
        ui_label: Minimum read length after trimming
      trim_poly_g:
        type: enum
        choices: [auto, on, off]
        default: auto
        exposure: user
        ui_label: Trim polyG tails
      trim_poly_x:
        type: boolean
        default: false
        exposure: user
        ui_label: Trim polyX tails
      adapter_sequence_r1:
        type: string
        default: null
        exposure: expert
      adapter_sequence_r2:
        type: string
        default: null
        exposure: expert
      cutadapt_error_rate:
        type: float
        default: 0.1
        min: 0
        max: 0.3
        exposure: expert
      cutadapt_min_overlap:
        type: integer
        default: 3
        min: 1
        max: 20
        exposure: expert
      threads:
        type: integer
        default: 4
        exposure: internal
    next_steps:
      recommended:
        - rnaseq_05_post_trim_quality_control
      allowed:
        - rnaseq_05_post_trim_quality_control
        - rnaseq_07_salmon_quantification
        - rnaseq_09_star_alignment
    resources:
      cpus: 4
      memory: 8.GB
      time: 4.h
    containers:
      fastp: quay.io/biocontainers/fastp:0.23.4--h5f740d0_0
      cutadapt: quay.io/biocontainers/cutadapt:5.1--py313h1a76870_0
    ui:
      card_title: Trim adapters and low-quality reads
      card_description: Clean reads before alignment or quantification.
      mode: configurable
      expose_in_guided: true
    chat_guidance:
      beginner: "I recommend fastp with Q20, minimum length 20, and automatic polyG handling."
      expert: "Use Cutadapt only when exact adapter grammar or older-study reproducibility is required."
```

## Workflow Request Format

Implement a run request file that the UI/chat sends to the builder:

```yaml
run_id: run001
omics: rnaseq
execution_profile: local_docker
input: assets/samplesheet.csv
outdir: results/run001
selected_steps:
  - rnaseq_03_raw_read_qc
  - rnaseq_04_adapter_quality_trimming
params:
  qc_tool: fastqc
  trimming_tool: fastp
  quality_threshold: 20
  minimum_read_length: 20
  trim_poly_g: auto
  trim_poly_x: false
  threads: 4
```

The builder must support:

```bash
# Run one step only
survom build-run --request examples/rnaseq_qc_only.run.yaml

# Run selected chain
survom build-run --request examples/rnaseq_qc_trim.run.yaml

# Generate but do not run
survom compile --request examples/rnaseq_qc_trim.run.yaml --output workflows/generated/run001

# Run compiled Nextflow directly
nextflow run workflows/generated/run001/main.nf -profile local,docker -params-file workflows/generated/run001/params.yaml -resume
```

## Planner Requirements

Implement `builder/planner.py` with these functions:

```python
def get_available_first_steps(omics: str) -> list[Step]:
    """Return steps that can start a workflow for this omics type."""


def suggest_next_steps(selected_steps: list[str], omics: str) -> list[StepSuggestion]:
    """Given current selected steps, return allowed next steps with recommendation labels."""


def validate_selected_steps(selected_steps: list[str], omics: str) -> ValidationResult:
    """Check dependencies, input/output compatibility, cycles, duplicate steps, and unsupported omics."""


def resolve_effective_params(selected_steps: list[str], user_params: dict) -> dict:
    """Merge defaults from registry with user overrides and profile-level resources."""
```

The planner must return errors like:

```json
{
  "valid": false,
  "errors": [
    {
      "step_id": "rnaseq_04_adapter_quality_trimming",
      "message": "Trimming requires FASTQ input. Add FASTQ staging or provide raw FASTQ manifest."
    }
  ],
  "suggested_fixes": [
    "Add rnaseq_02_fastq_staging_integrity before trimming."
  ]
}
```

## Compiler Requirements

Implement `builder/compiler.py` that generates:

- `main.nf`
- `nextflow.config`
- `params.yaml`
- `workflow_manifest.json`

Use Jinja2 templates.

The compiler should not generate tool-specific logic in Python. Instead, each registry step points to a Nextflow module and process name. The compiler includes selected modules and wires outputs to inputs.

For the first version, support linear workflows. Later the system can support branching and merging.

## Generated Nextflow Pattern

For selected steps `[rnaseq_03_raw_read_qc, rnaseq_04_adapter_quality_trimming]`, generate a `main.nf` similar to:

```nextflow
nextflow.enable.dsl = 2

include { FASTQC } from '../../modules/common/fastqc/main'
include { FASTP  } from '../../modules/common/fastp/main'
include { MAKE_PREPROCESS_MANIFEST } from '../../modules/common/make_preprocess_manifest/main'

workflow {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample_id,
                single_end: row.library_layout == 'single',
                platform: row.platform ?: 'unknown'
            ]
            def reads = meta.single_end ? [ file(row.fastq_1) ] : [ file(row.fastq_1), file(row.fastq_2) ]
            tuple(meta, reads)
        }
        .set { ch_reads }

    FASTQC(ch_reads)

    FASTP(ch_reads)

    MAKE_PREPROCESS_MANIFEST(FASTP.out.reads)
}
```

For single-step raw QC only, generate:

```nextflow
nextflow.enable.dsl = 2

include { FASTQC } from '../../modules/common/fastqc/main'

workflow {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample_id, single_end: row.library_layout == 'single' ]
            def reads = meta.single_end ? [ file(row.fastq_1) ] : [ file(row.fastq_1), file(row.fastq_2) ]
            tuple(meta, reads)
        }
        .set { ch_reads }

    FASTQC(ch_reads)
}
```

For trimming-only, generate:

```nextflow
nextflow.enable.dsl = 2

include { FASTP } from '../../modules/common/fastp/main'
include { MAKE_PREPROCESS_MANIFEST } from '../../modules/common/make_preprocess_manifest/main'

workflow {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample_id, single_end: row.library_layout == 'single', platform: row.platform ?: 'unknown' ]
            def reads = meta.single_end ? [ file(row.fastq_1) ] : [ file(row.fastq_1), file(row.fastq_2) ]
            tuple(meta, reads)
        }
        .set { ch_reads }

    FASTP(ch_reads)
    MAKE_PREPROCESS_MANIFEST(FASTP.out.reads)
}
```

## Required Nextflow Modules

Create reusable modules with clean input/output contracts.

### `modules/common/fastqc/main.nf`

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
    def adapterArg = params.custom_adapters ? "--adapters ${params.custom_adapters}" : ''
    def contaminantArg = params.custom_contaminants ? "--contaminants ${params.custom_contaminants}" : ''
    """
    fastqc \
      --threads ${task.cpus} \
      --outdir . \
      --extract \
      ${adapterArg} \
      ${contaminantArg} \
      ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """
}
```

### `modules/common/falco/main.nf`

```nextflow
process FALCO {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.txt"),  emit: txt
    path "versions.yml", emit: versions

    script:
    """
    falco \
      --outdir . \
      ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        falco: \$(falco --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}
```

### `modules/common/fastp/main.nf`

```nextflow
process FASTP {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.fastp.json"),       emit: json
    tuple val(meta), path("*.fastp.html"),       emit: html
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    def polyG = params.trim_poly_g == 'on' ? '--trim_poly_g' : ''
    def polyX = params.trim_poly_x ? '--trim_poly_x' : ''
    def adapterR1 = params.adapter_sequence_r1 ? "--adapter_sequence ${params.adapter_sequence_r1}" : ''
    def adapterR2 = params.adapter_sequence_r2 ? "--adapter_sequence_r2 ${params.adapter_sequence_r2}" : ''
    def detectPE = !meta.single_end && !params.adapter_sequence_r1 && !params.adapter_sequence_r2 ? '--detect_adapter_for_pe' : ''

    if (meta.single_end) {
        """
        fastp \
          --in1 ${reads[0]} \
          --out1 ${prefix}.trimmed.fastq.gz \
          --qualified_quality_phred ${params.quality_threshold} \
          --length_required ${params.minimum_read_length} \
          ${polyG} ${polyX} ${adapterR1} \
          --thread ${task.cpus} \
          --html ${prefix}.fastp.html \
          --json ${prefix}.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    } else {
        """
        fastp \
          --in1 ${reads[0]} \
          --in2 ${reads[1]} \
          --out1 ${prefix}_R1.trimmed.fastq.gz \
          --out2 ${prefix}_R2.trimmed.fastq.gz \
          --qualified_quality_phred ${params.quality_threshold} \
          --length_required ${params.minimum_read_length} \
          ${detectPE} ${polyG} ${polyX} ${adapterR1} ${adapterR2} \
          --thread ${task.cpus} \
          --html ${prefix}.fastp.html \
          --json ${prefix}.fastp.json

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    }
}
```

### `modules/common/cutadapt/main.nf`

```nextflow
process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("*.cutadapt.json"),    emit: json
    path "*.cutadapt.log", emit: log
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    if (meta.single_end) {
        """
        cutadapt \
          -j ${task.cpus} \
          -a ${params.adapter_sequence_r1} \
          -q ${params.quality_threshold} \
          -m ${params.minimum_read_length} \
          -e ${params.cutadapt_error_rate} \
          -O ${params.cutadapt_min_overlap} \
          --json ${prefix}.cutadapt.json \
          -o ${prefix}.trimmed.fastq.gz \
          ${reads[0]} \
          > ${prefix}.cutadapt.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        cutadapt \
          -j ${task.cpus} \
          -a ${params.adapter_sequence_r1} \
          -A ${params.adapter_sequence_r2} \
          -q ${params.quality_threshold} \
          -m ${params.minimum_read_length} \
          -e ${params.cutadapt_error_rate} \
          -O ${params.cutadapt_min_overlap} \
          --json ${prefix}.cutadapt.json \
          -o ${prefix}_R1.trimmed.fastq.gz \
          -p ${prefix}_R2.trimmed.fastq.gz \
          ${reads[0]} ${reads[1]} \
          > ${prefix}.cutadapt.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
}
```

## Configuration Requirements

Create `nextflow.config` with local, docker, singularity, and AWS profiles:

```groovy
profiles {
  local {
    process.executor = 'local'
  }

  docker {
    docker.enabled = true
  }

  singularity {
    singularity.enabled = true
    singularity.autoMounts = true
  }

  aws {
    process.executor = 'awsbatch'
    process.queue = params.aws_queue
    aws.region = params.aws_region
    docker.enabled = true
    workDir = params.aws_workdir
  }
}

process {
  withName:FASTQC {
    cpus = params.threads ?: 4
    memory = '8 GB'
    container = 'biocontainers/fastqc:v0.12.1_cv4'
  }

  withName:FALCO {
    cpus = params.threads ?: 4
    memory = '4 GB'
    container = 'quay.io/biocontainers/falco:1.3.2--h077b44d_0'
  }

  withName:FASTP {
    cpus = params.threads ?: 4
    memory = '8 GB'
    container = 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'
  }

  withName:CUTADAPT {
    cpus = params.threads ?: 4
    memory = '8 GB'
    container = 'quay.io/biocontainers/cutadapt:5.1--py313h1a76870_0'
  }
}
```

## UI/API Requirements

Implement API endpoints:

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

`POST /api/workflows/plan` request:

```json
{
  "omics": "rnaseq",
  "selected_steps": ["rnaseq_03_raw_read_qc"],
  "current_outputs": ["raw_qc_summary"]
}
```

Response:

```json
{
  "valid": true,
  "selected_steps": ["rnaseq_03_raw_read_qc"],
  "recommended_next_steps": [
    {
      "step_id": "rnaseq_04_adapter_quality_trimming",
      "label": "Adapter and Quality Trimming",
      "reason": "Raw QC can identify adapter contamination; trimming cleans reads before alignment or quantification."
    }
  ],
  "allowed_next_steps": [
    "rnaseq_04_adapter_quality_trimming",
    "rnaseq_07_salmon_quantification",
    "rnaseq_09_star_alignment"
  ]
}
```

## Builder Code Skeleton

Implement `builder/models.py` using Pydantic:

```python
from pydantic import BaseModel, Field
from typing import Any, Literal

Exposure = Literal['user', 'expert', 'internal']

class ParameterSpec(BaseModel):
    type: str
    default: Any = None
    choices: list[Any] | None = None
    min: float | None = None
    max: float | None = None
    exposure: Exposure
    ui_label: str | None = None

class OutputSpec(BaseModel):
    name: str
    type: str
    pattern: str | None = None
    ui_visible: bool = True

class InputSpec(BaseModel):
    name: str
    type: str
    formats: list[str] = Field(default_factory=list)
    required: bool = True
    required_when: str | None = None

class StepSpec(BaseModel):
    step_id: str
    step_name: str
    category: str
    omics: list[str]
    purpose: str
    module_path: str
    process_name: str
    default_tool: str
    fallback_tool: str | None = None
    can_run_individually: bool = True
    can_be_skipped: bool = False
    depends_on: list[str] = Field(default_factory=list)
    optional_depends_on: list[str] = Field(default_factory=list)
    compatible_previous_outputs: list[str] = Field(default_factory=list)
    required_inputs: list[InputSpec] = Field(default_factory=list)
    outputs: list[OutputSpec] = Field(default_factory=list)
    parameters: dict[str, ParameterSpec] = Field(default_factory=dict)
    next_steps: dict[str, list[str]] = Field(default_factory=dict)
    resources: dict[str, Any] = Field(default_factory=dict)
    containers: dict[str, str] = Field(default_factory=dict)
    ui: dict[str, Any] = Field(default_factory=dict)
    chat_guidance: dict[str, str] = Field(default_factory=dict)

class RunRequest(BaseModel):
    run_id: str
    omics: str
    execution_profile: str
    input: str
    outdir: str
    selected_steps: list[str]
    params: dict[str, Any] = Field(default_factory=dict)

class ValidationResult(BaseModel):
    valid: bool
    errors: list[dict[str, Any]] = Field(default_factory=list)
    warnings: list[dict[str, Any]] = Field(default_factory=list)
    suggested_fixes: list[str] = Field(default_factory=list)
```

Implement `builder/planner.py`:

```python
from .models import StepSpec, ValidationResult

class WorkflowPlanner:
    def __init__(self, steps: dict[str, StepSpec]):
        self.steps = steps

    def get_steps_for_omics(self, omics: str) -> list[StepSpec]:
        return [s for s in self.steps.values() if omics in s.omics]

    def get_available_first_steps(self, omics: str) -> list[StepSpec]:
        candidates = self.get_steps_for_omics(omics)
        return [s for s in candidates if s.can_run_individually]

    def suggest_next_steps(self, selected_steps: list[str], omics: str) -> dict:
        if not selected_steps:
            first = self.get_available_first_steps(omics)
            return {
                'recommended': first[:5],
                'allowed': first,
            }

        last_id = selected_steps[-1]
        last = self.steps[last_id]
        recommended_ids = last.next_steps.get('recommended', [])
        allowed_ids = last.next_steps.get('allowed', [])

        return {
            'recommended': [self.steps[sid] for sid in recommended_ids if sid in self.steps and omics in self.steps[sid].omics],
            'allowed': [self.steps[sid] for sid in allowed_ids if sid in self.steps and omics in self.steps[sid].omics],
        }

    def validate_selected_steps(self, selected_steps: list[str], omics: str) -> ValidationResult:
        errors = []
        warnings = []

        seen = set()
        for idx, step_id in enumerate(selected_steps):
            if step_id not in self.steps:
                errors.append({'step_id': step_id, 'message': 'Unknown step_id'})
                continue

            step = self.steps[step_id]

            if omics not in step.omics:
                errors.append({'step_id': step_id, 'message': f'Step is not compatible with omics={omics}'})

            if step_id in seen:
                errors.append({'step_id': step_id, 'message': 'Duplicate step in workflow'})
            seen.add(step_id)

            if idx > 0:
                prev_id = selected_steps[idx - 1]
                prev = self.steps.get(prev_id)
                if prev and step_id not in prev.next_steps.get('allowed', []) and step_id not in prev.next_steps.get('recommended', []):
                    warnings.append({
                        'step_id': step_id,
                        'message': f'{step_id} is not listed as a normal next step after {prev_id}; verify input/output compatibility.'
                    })

            missing_deps = [d for d in step.depends_on if d not in selected_steps[:idx]]
            if missing_deps and not step.can_run_individually:
                errors.append({'step_id': step_id, 'message': f'Missing dependencies: {missing_deps}'})

        return ValidationResult(valid=len(errors) == 0, errors=errors, warnings=warnings)

    def resolve_effective_params(self, selected_steps: list[str], user_params: dict) -> dict:
        resolved = {}
        for step_id in selected_steps:
            step = self.steps[step_id]
            for name, spec in step.parameters.items():
                if name not in resolved:
                    resolved[name] = spec.default
        resolved.update(user_params)
        return resolved
```

Implement `builder/compiler.py`:

```python
from pathlib import Path
import json
import yaml
from jinja2 import Environment, FileSystemLoader
from .models import RunRequest, StepSpec
from .planner import WorkflowPlanner

class NextflowCompiler:
    def __init__(self, registry_steps: dict[str, StepSpec], template_dir: str):
        self.steps = registry_steps
        self.planner = WorkflowPlanner(registry_steps)
        self.env = Environment(loader=FileSystemLoader(template_dir), trim_blocks=True, lstrip_blocks=True)

    def compile(self, request: RunRequest, output_dir: str) -> Path:
        validation = self.planner.validate_selected_steps(request.selected_steps, request.omics)
        if not validation.valid:
            raise ValueError(validation.model_dump_json(indent=2))

        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        selected = [self.steps[sid] for sid in request.selected_steps]
        params = self.planner.resolve_effective_params(request.selected_steps, request.params)
        params['input'] = request.input
        params['outdir'] = request.outdir
        params['run_id'] = request.run_id

        main_template = self.env.get_template('main.nf.j2')
        config_template = self.env.get_template('nextflow.config.j2')

        (out / 'main.nf').write_text(main_template.render(request=request, steps=selected), encoding='utf-8')
        (out / 'nextflow.config').write_text(config_template.render(request=request, steps=selected, params=params), encoding='utf-8')
        (out / 'params.yaml').write_text(yaml.safe_dump(params, sort_keys=False), encoding='utf-8')

        manifest = {
            'run_id': request.run_id,
            'omics': request.omics,
            'selected_steps': request.selected_steps,
            'params': params,
            'generated_files': ['main.nf', 'nextflow.config', 'params.yaml'],
        }
        (out / 'workflow_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
        return out
```

Create `workflows/templates/main.nf.j2`:

```jinja2
nextflow.enable.dsl = 2

{% for step in steps %}
include { {{ step.process_name }} } from '../../{{ step.module_path | replace('/main.nf', '/main') }}'
{% endfor %}
{% if 'rnaseq_04_adapter_quality_trimming' in request.selected_steps %}
include { MAKE_PREPROCESS_MANIFEST } from '../../modules/common/make_preprocess_manifest/main'
{% endif %}

workflow {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample_id,
                single_end: row.library_layout == 'single',
                platform: row.platform ?: 'unknown'
            ]
            def reads = meta.single_end ? [ file(row.fastq_1) ] : [ file(row.fastq_1), file(row.fastq_2) ]
            tuple(meta, reads)
        }
        .set { ch_reads }

{% for step in steps %}
{% if step.step_id == 'rnaseq_03_raw_read_qc' %}
    {{ step.process_name }}(ch_reads)
{% elif step.step_id == 'rnaseq_04_adapter_quality_trimming' %}
    {{ step.process_name }}(ch_reads)
    MAKE_PREPROCESS_MANIFEST({{ step.process_name }}.out.reads)
{% else %}
    // TODO: wire {{ step.step_id }} using registry input/output contracts
    {{ step.process_name }}(ch_previous)
{% endif %}

{% endfor %}
}
```

Create `workflows/templates/nextflow.config.j2`:

```jinja2
params {
  input = null
  outdir = 'results'
  run_id = '{{ request.run_id }}'
  aws_queue = null
  aws_region = 'eu-north-1'
  aws_workdir = null
}

profiles {
  local {
    process.executor = 'local'
  }

  docker {
    docker.enabled = true
  }

  singularity {
    singularity.enabled = true
    singularity.autoMounts = true
  }

  aws {
    process.executor = 'awsbatch'
    process.queue = params.aws_queue
    aws.region = params.aws_region
    docker.enabled = true
    workDir = params.aws_workdir
  }
}

process {
{% for step in steps %}
  withName:{{ step.process_name }} {
    cpus = {{ step.resources.get('cpus', 4) }}
    memory = '{{ step.resources.get('memory', '8.GB').replace('.', ' ') }}'
{% if step.containers %}
    container = '{{ step.containers.get(step.default_tool) }}'
{% endif %}
  }
{% endfor %}
}
```

## Single-Step and Full-Workflow Run Examples

Single raw QC:

```bash
nextflow run workflows/generated/run_qc_only/main.nf \
  -profile local,docker \
  -params-file workflows/generated/run_qc_only/params.yaml \
  -resume
```

Trimming only:

```bash
nextflow run workflows/generated/run_trim_only/main.nf \
  -profile local,docker \
  -params-file workflows/generated/run_trim_only/params.yaml \
  -resume
```

QC plus trimming:

```bash
nextflow run workflows/generated/run_qc_trim/main.nf \
  -profile local,docker \
  -params-file workflows/generated/run_qc_trim/params.yaml \
  -resume
```

AWS:

```bash
nextflow run workflows/generated/run_qc_trim/main.nf \
  -profile aws \
  -params-file workflows/generated/run_qc_trim/params.yaml \
  --aws_queue survom-rnaseq-queue \
  --aws_region eu-north-1 \
  --aws_workdir s3://survom-rnaseq/work/run_qc_trim \
  -resume
```

## Product Behavior

The UI/chat should work like this:

1. User selects omics type: RNA-seq.
2. System shows possible first steps from registry.
3. User selects first step, e.g. Raw Read QC.
4. System validates required input and parameters.
5. System suggests next steps from `next_steps.recommended` and `next_steps.allowed`.
6. User adds Adapter and Quality Trimming.
7. User can edit workflow-level parameters and step-level parameters.
8. System compiles the selected workflow into Nextflow DSL2.
9. User runs locally or on AWS using the same compiled workflow and different profile.
10. Output contracts are read back into UI/chat result cards.

## Testing Requirements

Write tests for:

- Registry loads and validates against schema.
- Each step has unique `step_id`.
- `COMMON` steps are usable by more than one omics type.
- User/expert/internal parameter exposure works.
- Planner recommends correct next steps after step 03.
- Planner allows step 04 to run individually if FASTQ manifest is supplied.
- Compiler generates valid `main.nf` for:
  - raw QC only
  - trimming only
  - raw QC + trimming
- Generated `params.yaml` includes defaults plus user overrides.
- AWS profile requires `aws_queue`, `aws_region`, and `aws_workdir`.

## Acceptance Criteria

The implementation is complete when:

1. A user can run one module only.
2. A user can build a selected chain of modules.
3. The UI/chat can ask “what next?” and receive valid recommended next steps.
4. Common modules are reusable across RNA-seq and genomics workflows without rewriting them.
5. Nextflow module code is separated from workflow orchestration.
6. AWS vs local execution is controlled by profiles, not different pipeline code.
7. Every module has machine-readable outputs for UI/chat.
8. Every run writes a provenance manifest containing selected steps, parameters, software versions, input paths, output paths, and execution profile.
