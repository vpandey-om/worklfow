# RNA-seq Tool Selection Policy For Survom

## Purpose

Survom is a production multiomics platform, not a tool encyclopedia. Each workflow step should expose one recommended default and, only when justified, one fallback or special-case tool.

The policy below keeps workflows understandable for guided users, predictable for production support, and reusable across multiomics assays.

## Core Rule

Each step may define at most:

- **Default tool / algorithm**
- **Fallback tool / algorithm**

If one tool is enough, no fallback is required.

## Selection Criteria

Prefer tools that have:

- stable, machine-readable outputs,
- MultiQC support where relevant,
- strong community adoption,
- simple parameter surface,
- active maintenance,
- mature Docker/Conda packaging,
- predictable resource use,
- clear Nextflow DSL2 module boundaries,
- strong compatibility with downstream steps.

Avoid exposing tools only because they are popular in the ecosystem. More options create more validation cases, more UI complexity, more documentation burden, and more support risk.

## Supported RNA-seq MVP Routes

### Route A: Alignment-Based Gene-Level Analysis

```text
FASTQ -> QC -> trimming -> STAR -> samtools sort/index -> featureCounts -> DESeq2 -> enrichment
```

Use this when:

- BAM files are required,
- genome alignment QC is important,
- splice junction information matters,
- downstream visualization in IGV is expected,
- users need conventional gene-level counts from aligned reads.

### Route B: Pseudoalignment Transcript/Gene-Level Analysis

```text
FASTQ -> QC -> trimming -> Salmon -> tximport -> DESeq2 -> enrichment
```

Use this when:

- the primary goal is expression quantification,
- speed and lower storage are important,
- BAM output is not required,
- transcript-level estimates may be useful.

Default Survom route for standard bulk RNA-seq MVP:

```text
Route B: Salmon -> tximport -> DESeq2
```

Reason: faster, easier to scale, less storage-heavy, and sufficient for many differential expression workflows.

## Curated Tool Choices

| Task | Default | Fallback | Rationale |
|---|---|---|---|
| Project initialization | custom validator | none | Platform-specific |
| Sample sheet validation | Python/Pydantic or JSON Schema | none | Keep schema enforcement centralized |
| FASTQ integrity | Nextflow staging + md5/gzip checks | seqkit stats | seqkit only when deeper FASTQ stats are needed |
| Raw FASTQ QC | FastQC | Falco | Falco only for speed/large-scale cases |
| Trimming | fastp | Cutadapt | fastp is simple; Cutadapt handles exact custom adapters |
| Post-trim QC | FastQC | Falco | Same as raw QC |
| Strandedness inference | Salmon auto lib type | RSeQC infer_experiment.py | RSeQC requires BAM, so fallback only |
| Reference validation | custom reference bundle validator | none | Platform-specific |
| Alignment | STAR | HISAT2 | STAR default for bulk RNA-seq; HISAT2 for lower memory/special cases |
| BAM sort/index | samtools | none | Stable standard |
| Gene counting | featureCounts | none | Stable, widely used |
| Pseudoalignment | Salmon | Kallisto | Salmon default for richer output and mature workflow use |
| Transcript-to-gene summarization | tximport | tximeta | tximeta when provenance-aware transcriptome metadata is needed |
| Differential expression | DESeq2 | edgeR | DESeq2 default; edgeR for specific designs/count models |
| ORA enrichment | clusterProfiler ORA | g:Profiler API/offline export | ORA is a distinct analytical task |
| GSEA enrichment | fgsea | clusterProfiler GSEA | GSEA is separate from ORA |
| PCA/sample QC | DESeq2/R base PCA | plotly export | Plotly only for interactive product view |
| Volcano plot | EnhancedVolcano | ggplot2 custom | EnhancedVolcano default for standard output |
| Heatmap | pheatmap | ComplexHeatmap | ComplexHeatmap for complex annotations |
| MultiQC aggregation | MultiQC | none | Production standard |
| Report rendering | Quarto/R Markdown | none | Reproducible reporting |
| Provenance capture | custom run manifest | RO-Crate export | RO-Crate for advanced interoperability |

## Parameter Exposure Policy

Use three parameter levels:

### User-Settable Parameters

Visible in guided and dynamic UI.

Examples:

- organism,
- genome build,
- library layout,
- strandedness,
- route choice,
- condition column,
- batch column,
- FDR threshold.

### Expert Parameters

Hidden under advanced settings.

Examples:

- STAR index overhang,
- Salmon library type override,
- featureCounts strand mode,
- DESeq2 design formula,
- minimum count filtering threshold.

### Internal Parameters

Controlled by platform defaults and not normally shown.

Examples:

- Nextflow channel names,
- publish directory structure,
- container paths,
- retry policy,
- temporary work directory.

## Verification Policy

If a command-line flag, default, or tool behavior is uncertain, write:

```text
VERIFY_IN_OFFICIAL_DOCS
```

Do not invent exact flags or defaults.

## Multiomics Reuse Tags

Every step should be tagged as one of:

- `COMMON`
- `RNA_SEQ_SPECIFIC`
- `GENOMICS_COMMON`
- `STATISTICS_COMMON`
- `REPORTING_COMMON`

This supports reusing common steps across RNA-seq, genomics, ATAC-seq, ChIP-seq, scRNA-seq, proteomics where appropriate, metabolomics where applicable, microbiome, and methylation workflows.

