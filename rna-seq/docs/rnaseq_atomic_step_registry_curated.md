# Curated RNA-seq Atomic Step Registry For Survom

This registry defines a clean production-oriented RNA-seq workflow model for Survom.

It intentionally avoids listing every possible tool. Each step has one default tool or algorithm and, only when useful, one fallback or special-case tool.

Supported MVP routes:

- **Route A:** FASTQ -> QC -> trimming -> STAR -> samtools -> featureCounts -> DESeq2 -> enrichment
- **Route B:** FASTQ -> QC -> trimming -> Salmon -> tximport -> DESeq2 -> enrichment

Both routes require sample-sheet validation and reference-bundle validation before quantification or alignment.

Default route for standard bulk RNA-seq:

```text
Route B: Salmon -> tximport -> DESeq2
```

---

## 00 Project Initialization

```yaml
step_id: rnaseq_00_project_initialization
step_name: Project initialization
category: COMMON
purpose: Create a reproducible run context before any data processing starts.
why_this_step_exists: Production runs need a stable run ID, output directory, execution profile, provenance record, and parameter snapshot.
required_inputs:
  - project_name
  - organism
  - genome_build
  - outdir
optional_inputs:
  - profile
  - resume
  - container_engine
input_formats:
  - UI form
  - JSON
  - YAML
required_outputs:
  - run_manifest
  - project_run_directory
  - resolved_params
optional_outputs:
  - audit_event
output_formats:
  - JSON
  - directory
depends_on: []
can_be_skipped: false
skip_conditions: []
default_tool_or_algorithm: custom project config validator
fallback_tool_or_algorithm: none
tool_selection_logic: Always use the platform validator.
user_settable_parameters:
  - project_name
  - organism
  - genome_build
  - outdir
  - profile
  - resume
  - container_engine
  - max_cpus
  - max_memory
  - max_time
expert_parameters:
  - work_dir
  - publish_mode
internal_parameters:
  - run_id
  - audit_actor
  - registry_version
qc_metrics:
  - config_valid
qc_pass_fail_rules:
  - fail if required project metadata is missing
  - fail if output directory is not writable
common_failure_modes:
  - invalid organism/genome build pair
  - unwritable output directory
  - unsupported execution profile
nextflow_module_name: platform/init_project
reusable_for_other_omics: COMMON
downstream_steps:
  - rnaseq_01_sample_sheet_validation
notes_for_guided_flowchart: Show this as project setup, not as a scientific analysis step.
notes_for_dynamic_flowchart: Lock this as the root node.
notes_for_chat_mode: Chat may collect missing project metadata before drafting a workflow.
```

---

## 01 Sample Sheet Validation

```yaml
step_id: rnaseq_01_sample_sheet_validation
step_name: Sample sheet validation
category: COMMON
purpose: Validate sample metadata and map biological samples to input files.
why_this_step_exists: Most RNA-seq failures come from incorrect sample metadata, missing files, duplicated sample IDs, or invalid condition labels.
required_inputs:
  - sample_sheet
optional_inputs:
  - project_metadata
input_formats:
  - CSV
  - TSV
required_outputs:
  - validated_sample_sheet
  - sample_manifest
optional_outputs:
  - validation_report
output_formats:
  - CSV
  - JSON
  - HTML
depends_on:
  - rnaseq_00_project_initialization
can_be_skipped: false
skip_conditions: []
default_tool_or_algorithm: Python/Pydantic or JSON Schema validator
fallback_tool_or_algorithm: none
tool_selection_logic: Use one platform schema implementation to keep validation consistent across UI, API, and Nextflow.
user_settable_parameters:
  - sample_sheet
  - sample_id
  - fastq_1
  - fastq_2
  - condition
  - batch
  - replicate
  - strandedness
  - library_layout
expert_parameters:
  - allow_single_replicate_warning
  - allow_missing_batch_column
internal_parameters:
  - schema_version
  - file_resolver
qc_metrics:
  - number_of_samples
  - number_of_conditions
  - replicates_per_condition
  - missing_fastq_count
  - duplicate_sample_id_count
qc_pass_fail_rules:
  - fail if FASTQ paths do not exist
  - fail if sample IDs are not unique
  - fail if condition column is missing
  - fail if paired-end files are not correctly paired
  - warn if fewer than 3 biological replicates per condition
  - warn if batch column is missing
  - warn or fail if batch is confounded with condition
common_failure_modes:
  - sample IDs contain spaces or unsafe characters
  - condition labels inconsistent by capitalization
  - R1/R2 naming mismatch
  - metadata says paired-end but only one FASTQ is present
nextflow_module_name: common/validate_samplesheet
reusable_for_other_omics: COMMON
downstream_steps:
  - rnaseq_02_fastq_staging_integrity
notes_for_guided_flowchart: Show actionable errors in plain language.
notes_for_dynamic_flowchart: This node must produce a typed sample_manifest before FASTQ steps unlock.
notes_for_chat_mode: Chat can ask for missing condition, batch, or replicate columns.
```

---

## 02 FASTQ Staging And Integrity

```yaml
step_id: rnaseq_02_fastq_staging_integrity
step_name: FASTQ staging and integrity
category: COMMON
purpose: Stage FASTQ files into the workflow and confirm they are readable before compute-heavy steps.
why_this_step_exists: Broken gzip files, incomplete uploads, or path resolution errors should fail early.
required_inputs:
  - sample_manifest
  - raw_fastq_files
optional_inputs:
  - md5_manifest
input_formats:
  - FASTQ.GZ
  - FASTQ
  - MD5 text
required_outputs:
  - staged_fastq_manifest
  - integrity_report
optional_outputs:
  - fastq_stats
output_formats:
  - JSON
  - TSV
depends_on:
  - rnaseq_01_sample_sheet_validation
can_be_skipped: false
skip_conditions: []
default_tool_or_algorithm: Nextflow channel staging plus md5/gzip integrity check
fallback_tool_or_algorithm: seqkit stats
tool_selection_logic: Use built-in staging and lightweight integrity checks by default; use seqkit when deeper read statistics are needed.
user_settable_parameters:
  - verify_md5
expert_parameters:
  - run_seqkit_stats
internal_parameters:
  - staged_file_paths
  - storage_backend
qc_metrics:
  - files_staged
  - gzip_integrity_failures
  - md5_mismatches
  - read_count_estimate
qc_pass_fail_rules:
  - fail if any required FASTQ cannot be staged
  - fail if gzip integrity check fails
  - fail if provided md5 does not match
common_failure_modes:
  - interrupted upload
  - incorrect object storage path
  - uncompressed FASTQ too large for configured policy
nextflow_module_name: common/stage_fastq_integrity
reusable_for_other_omics: COMMON
downstream_steps:
  - rnaseq_03_raw_read_qc
notes_for_guided_flowchart: Hide implementation details; show "files verified".
notes_for_dynamic_flowchart: Emits typed fastq_pairs or fastq_single.
notes_for_chat_mode: Chat should explain missing or corrupt files before workflow generation.
```

---

## 03 Raw Read Quality Control

```yaml
step_id: rnaseq_03_raw_read_qc
step_name: Raw read quality control
category: COMMON
purpose: Assess sequencing quality before modifying reads.
why_this_step_exists: Raw QC identifies adapter contamination, low-quality cycles, GC bias, duplication, and overrepresented sequences.
required_inputs:
  - staged_fastq_manifest
optional_inputs:
  - custom_adapters
input_formats:
  - FASTQ.GZ
  - FASTQ
required_outputs:
  - fastqc_html
  - fastqc_zip
optional_outputs:
  - raw_qc_summary
output_formats:
  - HTML
  - ZIP
  - JSON
depends_on:
  - rnaseq_02_fastq_staging_integrity
can_be_skipped: false
skip_conditions: []
default_tool_or_algorithm: FastQC
fallback_tool_or_algorithm: Falco
tool_selection_logic: Use FastQC by default for mature MultiQC compatibility; use Falco for large cohorts where FastQC runtime is the bottleneck.
user_settable_parameters:
  - threads
expert_parameters:
  - custom_adapters
  - custom_contaminants
internal_parameters:
  - output_prefix
qc_metrics:
  - per_base_quality_status
  - adapter_content_status
  - gc_content_status
  - duplication_status
  - overrepresented_sequences_status
qc_pass_fail_rules:
  - warn if adapter content is high
  - warn if per-base quality degrades strongly at read ends
  - warn if GC distribution is abnormal
  - do not fail solely on FastQC warnings without user review
common_failure_modes:
  - wrong file type
  - corrupted FASTQ
  - unexpected library contamination
nextflow_module_name: common/fastqc
reusable_for_other_omics: COMMON
downstream_steps:
  - rnaseq_04_adapter_quality_trimming
  - rnaseq_07_salmon_quantification
  - rnaseq_09_star_alignment
notes_for_guided_flowchart: Explain QC warnings but do not overload users with all FastQC modules.
notes_for_dynamic_flowchart: This node can feed trimming or direct quantification if trimming is skipped.
notes_for_chat_mode: Chat should summarize QC issues and recommend whether trimming is needed.
```

---

## 04 Adapter And Quality Trimming

```yaml
step_id: rnaseq_04_adapter_quality_trimming
step_name: Adapter and quality trimming
category: COMMON
purpose: Remove adapters and low-quality bases that can reduce quantification or alignment quality.
why_this_step_exists: Adapter contamination and poor-quality read tails can bias downstream mapping and counting.
required_inputs:
  - staged_fastq_manifest
optional_inputs:
  - custom_adapter_sequences
input_formats:
  - FASTQ.GZ
  - FASTQ
required_outputs:
  - trimmed_fastq_manifest
  - trimming_report
optional_outputs:
  - trimming_html
  - trimming_json
output_formats:
  - FASTQ.GZ
  - HTML
  - JSON
depends_on:
  - rnaseq_03_raw_read_qc
can_be_skipped: true
skip_conditions:
  - raw QC shows no adapter contamination and high read quality
  - user selects no trimming with expert approval
default_tool_or_algorithm: fastp
fallback_tool_or_algorithm: Cutadapt
tool_selection_logic: fastp is default because it combines adapter detection, quality filtering, and machine-readable reports. Cutadapt is fallback for exact custom adapters or reproducibility with older studies.
user_settable_parameters:
  - minimum_read_length
  - quality_threshold
  - trim_poly_g
  - trim_poly_x
expert_parameters:
  - explicit_adapter_r1
  - explicit_adapter_r2
  - front_trim_r1
  - front_trim_r2
internal_parameters:
  - threads
  - report_prefix
qc_metrics:
  - reads_before_filtering
  - reads_after_filtering
  - percent_reads_retained
  - bases_trimmed
  - adapter_trimmed_reads
qc_pass_fail_rules:
  - warn if too many reads are discarded
  - warn if adapter content remains high after trimming
  - fail if output FASTQ is empty
common_failure_modes:
  - over-trimming due to incorrect adapter sequence
  - minimum read length too strict
  - mixed library types in one sample sheet
nextflow_module_name: common/fastp
reusable_for_other_omics: COMMON
downstream_steps:
  - rnaseq_05_post_trim_qc
notes_for_guided_flowchart: Default to fastp and expose only simple quality/min-length controls.
notes_for_dynamic_flowchart: Trimming may be skipped if QC and route allow it.
notes_for_chat_mode: Chat should recommend trimming only when QC or protocol supports it.
```

---

## 05 Post-Trim Quality Control

```yaml
step_id: rnaseq_05_post_trim_qc
step_name: Post-trim quality control
category: COMMON
purpose: Confirm trimming improved reads and generated usable FASTQ files.
why_this_step_exists: Users need pre/post trimming comparison before alignment or quantification.
required_inputs:
  - trimmed_fastq_manifest
optional_inputs: []
input_formats:
  - FASTQ.GZ
required_outputs:
  - post_trim_fastqc_html
  - post_trim_fastqc_zip
optional_outputs:
  - post_trim_qc_summary
output_formats:
  - HTML
  - ZIP
  - JSON
depends_on:
  - rnaseq_04_adapter_quality_trimming
can_be_skipped: true
skip_conditions:
  - trimming was skipped
default_tool_or_algorithm: FastQC
fallback_tool_or_algorithm: Falco
tool_selection_logic: Same as raw read QC.
user_settable_parameters:
  - threads
expert_parameters: []
internal_parameters:
  - output_prefix
qc_metrics:
  - adapter_content_status_after_trimming
  - per_base_quality_status_after_trimming
qc_pass_fail_rules:
  - warn if adapter content remains high
  - warn if read quality remains poor
common_failure_modes:
  - trimming did not remove adapter due to unusual adapter sequence
  - trimmed reads too short
nextflow_module_name: common/fastqc_post_trim
reusable_for_other_omics: COMMON
downstream_steps:
  - rnaseq_06_strandedness_inference
  - rnaseq_07_salmon_quantification
  - rnaseq_09_star_alignment
notes_for_guided_flowchart: Show pre/post QC comparison.
notes_for_dynamic_flowchart: Output is QC-only; it does not change reads.
notes_for_chat_mode: Chat should explain whether trimming improved the data.
```

---

## 06 Strandedness Inference

```yaml
step_id: rnaseq_06_strandedness_inference
step_name: Strandedness inference
category: RNA_SEQ_SPECIFIC
purpose: Determine whether the library is unstranded, forward-stranded, or reverse-stranded.
why_this_step_exists: Incorrect strandedness causes systematically wrong quantification and can silently invalidate results.
required_inputs:
  - fastq_manifest
  - transcriptome_index
optional_inputs:
  - aligned_bam
  - bed12_annotation
input_formats:
  - FASTQ.GZ
  - Salmon index
  - BAM
  - BED12
required_outputs:
  - strandedness_call
  - strandedness_evidence
optional_outputs:
  - inference_log
output_formats:
  - JSON
  - TXT
depends_on:
  - rnaseq_05_post_trim_qc
can_be_skipped: true
skip_conditions:
  - strandedness is known and explicitly supplied in sample sheet
default_tool_or_algorithm: Salmon auto library type on representative subsample
fallback_tool_or_algorithm: RSeQC infer_experiment.py
tool_selection_logic: Salmon is default because it can infer library type without a BAM. RSeQC is fallback when BAM exists or transcriptome indexing is unavailable.
user_settable_parameters:
  - strandedness
expert_parameters:
  - subsample_read_count
  - representative_samples_per_batch
internal_parameters:
  - inference_threshold
qc_metrics:
  - compatible_fragment_ratio
  - strandedness_confidence
qc_pass_fail_rules:
  - fail if auto-inference is ambiguous and strandedness is required
  - warn if batches disagree in inferred strandedness
common_failure_modes:
  - wrong transcriptome index
  - low-complexity subsample
  - mixed library prep protocols
nextflow_module_name: rnaseq/infer_strandedness
reusable_for_other_omics: RNA_SEQ_SPECIFIC
downstream_steps:
  - rnaseq_07_salmon_quantification
  - rnaseq_11_featurecounts_gene_counting
notes_for_guided_flowchart: Present as "library orientation check".
notes_for_dynamic_flowchart: Required when strandedness is set to auto.
notes_for_chat_mode: Chat should ask for lab protocol if inference is ambiguous.
```

---

## 06A Reference Bundle Validation

```yaml
step_id: rnaseq_06a_reference_bundle_validation
step_name: Reference bundle validation
category: GENOMICS_COMMON
purpose: Confirm that genome, transcriptome, annotation, and prebuilt indexes are mutually compatible before alignment or quantification.
why_this_step_exists: Wrong genome builds, mismatched GTF files, transcript ID version mismatches, or stale indexes can silently corrupt RNA-seq results.
required_inputs:
  - organism
  - genome_build
  - selected_route
optional_inputs:
  - genome_fasta
  - gtf_annotation
  - star_index
  - transcriptome_fasta
  - salmon_index
  - tx2gene_mapping
input_formats:
  - FASTA
  - GTF
  - STAR index directory
  - Salmon index directory
  - CSV
  - TSV
required_outputs:
  - validated_reference_bundle
  - reference_validation_report
optional_outputs:
  - index_build_recommendation
output_formats:
  - JSON
  - HTML
depends_on:
  - rnaseq_00_project_initialization
can_be_skipped: false
skip_conditions: []
default_tool_or_algorithm: custom Survom reference bundle validator
fallback_tool_or_algorithm: none
tool_selection_logic: Reference compatibility must be enforced by the platform because it spans files, indexes, and metadata rather than one external tool.
user_settable_parameters:
  - organism
  - genome_build
  - selected_route
expert_parameters:
  - allow_custom_reference
  - allow_index_reuse
internal_parameters:
  - reference_registry_version
  - expected_chromosome_style
  - transcript_id_version_policy
qc_metrics:
  - reference_files_present
  - chromosome_name_compatibility
  - transcript_gene_mapping_rate
  - index_matches_reference
qc_pass_fail_rules:
  - fail if selected route lacks its required index
  - fail if genome FASTA and GTF chromosome naming are incompatible
  - fail if Salmon index, transcriptome, and tx2gene mapping are incompatible
  - warn if custom reference lacks provenance
common_failure_modes:
  - hg19/hg38 mismatch
  - Ensembl/UCSC chromosome naming mismatch
  - transcript IDs with version suffixes do not match tx2gene table
  - STAR or Salmon index built from a different reference release
nextflow_module_name: common/validate_reference_bundle
reusable_for_other_omics: GENOMICS_COMMON
downstream_steps:
  - rnaseq_07_salmon_quantification
  - rnaseq_09_star_alignment
notes_for_guided_flowchart: Show this as "reference checked" and surface only actionable mismatches.
notes_for_dynamic_flowchart: Must emit a validated reference bundle before alignment or quantification nodes unlock.
notes_for_chat_mode: Chat should ask for organism/genome build when missing and warn about custom reference provenance.
```

---

## 07 Salmon Quantification

```yaml
step_id: rnaseq_07_salmon_quantification
step_name: Transcript quantification with Salmon
category: RNA_SEQ_SPECIFIC
purpose: Quantify transcript abundance directly from FASTQ reads.
why_this_step_exists: Provides fast transcript-level quantification without generating large BAM files.
required_inputs:
  - analysis_fastq_manifest
  - transcriptome_index
  - strandedness_call
optional_inputs:
  - gene_transcript_map
input_formats:
  - FASTQ.GZ
  - Salmon index
  - JSON
required_outputs:
  - salmon_quant_files
  - salmon_logs
optional_outputs:
  - mapping_metrics
output_formats:
  - quant.sf
  - JSON
  - LOG
depends_on:
  - rnaseq_03_raw_read_qc
  - rnaseq_06_strandedness_inference
  - rnaseq_06a_reference_bundle_validation
can_be_skipped: true
skip_conditions:
  - user chooses alignment-based route
default_tool_or_algorithm: Salmon
fallback_tool_or_algorithm: Kallisto
tool_selection_logic: Salmon is default for mature output, bias correction support, and common tximport workflows. Kallisto is fallback for studies that require compatibility with existing Kallisto outputs.
user_settable_parameters:
  - quant_route
  - strandedness
expert_parameters:
  - bias_correction
  - gc_bias_correction
  - validate_mappings
internal_parameters:
  - threads
  - index_path
qc_metrics:
  - mapping_rate
  - fragments_processed
  - library_type_compatibility
qc_pass_fail_rules:
  - warn if mapping rate is low
  - fail if no quantification files are produced
  - fail if library type is incompatible with inferred strandedness
common_failure_modes:
  - wrong transcriptome index
  - incorrect strandedness
  - poor read quality
nextflow_module_name: rnaseq/salmon_quant
reusable_for_other_omics: RNA_SEQ_SPECIFIC
downstream_steps:
  - rnaseq_08_tximport_gene_summarization
notes_for_guided_flowchart: Default route for standard expression analysis.
notes_for_dynamic_flowchart: Emits transcript_quant and can feed tximport.
notes_for_chat_mode: Chat should choose this route for fast differential expression unless BAM is requested.
```

---

## 08 Transcript-To-Gene Summarization

```yaml
step_id: rnaseq_08_tximport_gene_summarization
step_name: Transcript-to-gene summarization
category: RNA_SEQ_SPECIFIC
purpose: Convert transcript-level quantification to gene-level counts suitable for differential expression.
why_this_step_exists: DESeq2 gene-level analysis needs a consistent gene-by-sample matrix.
required_inputs:
  - salmon_quant_files
  - tx2gene_mapping
optional_inputs:
  - transcriptome_metadata
input_formats:
  - quant.sf
  - CSV
  - TSV
required_outputs:
  - gene_count_matrix
  - abundance_matrix
optional_outputs:
  - transcript_count_matrix
output_formats:
  - CSV
  - TSV
  - RDS
depends_on:
  - rnaseq_07_salmon_quantification
can_be_skipped: true
skip_conditions:
  - user chooses alignment-based featureCounts route
default_tool_or_algorithm: tximport
fallback_tool_or_algorithm: tximeta
tool_selection_logic: tximport is simple and stable. tximeta is fallback when transcriptome provenance metadata is needed.
user_settable_parameters:
  - summarization_level
expert_parameters:
  - counts_from_abundance
internal_parameters:
  - tx2gene_source
qc_metrics:
  - genes_detected
  - samples_imported
  - missing_quant_files
qc_pass_fail_rules:
  - fail if any expected sample quantification is missing
  - fail if transcript IDs cannot map to genes
common_failure_modes:
  - transcript ID version mismatch
  - wrong annotation release
  - missing tx2gene mapping
nextflow_module_name: rnaseq/tximport
reusable_for_other_omics: RNA_SEQ_SPECIFIC
downstream_steps:
  - rnaseq_12_count_matrix_qc_filtering
notes_for_guided_flowchart: Explain this as "make gene count table".
notes_for_dynamic_flowchart: Accepts transcript_quant and emits gene_counts.
notes_for_chat_mode: Chat should include this automatically after Salmon when DESeq2 is requested.
```

---

## 09 STAR Alignment

```yaml
step_id: rnaseq_09_star_alignment
step_name: Splice-aware genome alignment
category: RNA_SEQ_SPECIFIC
purpose: Align RNA-seq reads to the genome and generate alignment files.
why_this_step_exists: Alignment route is needed for BAM outputs, genome-level QC, splice-aware analysis, and conventional counting.
required_inputs:
  - analysis_fastq_manifest
  - star_genome_index
  - strandedness_call
optional_inputs:
  - gtf_annotation
input_formats:
  - FASTQ.GZ
  - STAR index directory
  - GTF
required_outputs:
  - aligned_bam_unsorted_or_sorted
  - alignment_log
optional_outputs:
  - splice_junctions
output_formats:
  - BAM
  - LOG
  - SJ.out.tab
depends_on:
  - rnaseq_03_raw_read_qc
  - rnaseq_06_strandedness_inference
  - rnaseq_06a_reference_bundle_validation
can_be_skipped: true
skip_conditions:
  - user chooses Salmon route
default_tool_or_algorithm: STAR
fallback_tool_or_algorithm: HISAT2
tool_selection_logic: STAR is default for standard bulk RNA-seq because it is fast and mature for splice-aware alignment. HISAT2 is fallback for lower memory environments or specific legacy compatibility.
user_settable_parameters:
  - alignment_route
expert_parameters:
  - read_files_command
  - genome_index_overhang
  - two_pass_mode
internal_parameters:
  - threads
  - genome_index_path
qc_metrics:
  - uniquely_mapped_percent
  - multimapped_percent
  - unmapped_percent
  - chimeric_reads
qc_pass_fail_rules:
  - warn if uniquely mapped percent is low
  - fail if BAM output is missing
  - warn if high multimapping suggests contamination or wrong reference
common_failure_modes:
  - wrong genome build
  - missing or mismatched GTF
  - insufficient memory
nextflow_module_name: rnaseq/star_align
reusable_for_other_omics: RNA_SEQ_SPECIFIC
downstream_steps:
  - rnaseq_10_bam_sort_index
notes_for_guided_flowchart: Choose this route when users request BAM or genome alignment.
notes_for_dynamic_flowchart: Emits BAM and alignment metrics.
notes_for_chat_mode: Chat should not choose alignment route unless user needs BAM, splice, or alignment QC.
```

---

## 10 BAM Sort And Index

```yaml
step_id: rnaseq_10_bam_sort_index
step_name: BAM sorting and indexing
category: GENOMICS_COMMON
purpose: Produce coordinate-sorted and indexed BAM files for counting, QC, and visualization.
why_this_step_exists: Many downstream tools require sorted/indexed BAM files.
required_inputs:
  - aligned_bam
optional_inputs: []
input_formats:
  - BAM
required_outputs:
  - sorted_bam
  - bam_index
optional_outputs:
  - samtools_stats
output_formats:
  - BAM
  - BAI
  - TXT
depends_on:
  - rnaseq_09_star_alignment
can_be_skipped: true
skip_conditions:
  - aligner directly emits sorted and indexed BAM accepted by downstream tools
default_tool_or_algorithm: samtools sort/index
fallback_tool_or_algorithm: none
tool_selection_logic: samtools is the standard supported implementation.
user_settable_parameters: []
expert_parameters:
  - sort_memory
internal_parameters:
  - threads
  - temporary_directory
qc_metrics:
  - bam_records
  - index_created
qc_pass_fail_rules:
  - fail if sorted BAM or index is missing
common_failure_modes:
  - insufficient disk space
  - insufficient memory for sorting
nextflow_module_name: genomics/samtools_sort_index
reusable_for_other_omics: GENOMICS_COMMON
downstream_steps:
  - rnaseq_11_featurecounts_gene_counting
  - rnaseq_18_multiqc_aggregation
notes_for_guided_flowchart: Hide unless alignment route is expanded.
notes_for_dynamic_flowchart: Required edge between alignment and counting.
notes_for_chat_mode: Chat should include it automatically after STAR.
```

---

## 11 Gene Counting With featureCounts

```yaml
step_id: rnaseq_11_featurecounts_gene_counting
step_name: Gene-level read counting
category: RNA_SEQ_SPECIFIC
purpose: Count aligned reads overlapping annotated genes.
why_this_step_exists: Differential expression requires a gene-by-sample count matrix.
required_inputs:
  - sorted_bam
  - gtf_annotation
  - strandedness_call
optional_inputs: []
input_formats:
  - BAM
  - GTF
required_outputs:
  - gene_count_matrix
  - featurecounts_summary
optional_outputs: []
output_formats:
  - TSV
  - TXT
depends_on:
  - rnaseq_10_bam_sort_index
  - rnaseq_06_strandedness_inference
can_be_skipped: true
skip_conditions:
  - user chooses Salmon route
default_tool_or_algorithm: featureCounts
fallback_tool_or_algorithm: none
tool_selection_logic: featureCounts is stable, fast, widely used, and produces MultiQC-compatible summaries.
user_settable_parameters:
  - strandedness
expert_parameters:
  - feature_type
  - attribute_type
  - count_multimapping_reads
  - count_overlapping_features
internal_parameters:
  - threads
qc_metrics:
  - assigned_reads_percent
  - unassigned_no_feature_percent
  - unassigned_ambiguity_percent
qc_pass_fail_rules:
  - warn if assigned reads percent is low
  - fail if count matrix is empty
common_failure_modes:
  - wrong strandedness
  - wrong annotation release
  - chromosome naming mismatch between BAM and GTF
nextflow_module_name: rnaseq/featurecounts
reusable_for_other_omics: RNA_SEQ_SPECIFIC
downstream_steps:
  - rnaseq_12_count_matrix_qc_filtering
notes_for_guided_flowchart: Explain as "count reads per gene".
notes_for_dynamic_flowchart: Accepts sorted_bam and emits gene_counts.
notes_for_chat_mode: Chat should include this only for alignment route.
```

---

## 12 Count Matrix QC And Filtering

```yaml
step_id: rnaseq_12_count_matrix_qc_filtering
step_name: Count matrix QC and low-expression filtering
category: STATISTICS_COMMON
purpose: Check gene count matrix quality and remove genes with insufficient counts for reliable testing.
why_this_step_exists: Low-count genes inflate multiple testing burden and add noise.
required_inputs:
  - gene_count_matrix
  - validated_sample_sheet
optional_inputs: []
input_formats:
  - CSV
  - TSV
required_outputs:
  - filtered_count_matrix
  - count_qc_report
optional_outputs:
  - normalized_preview_matrix
output_formats:
  - CSV
  - TSV
  - HTML
depends_on:
  one_of:
    - rnaseq_08_tximport_gene_summarization
    - rnaseq_11_featurecounts_gene_counting
can_be_skipped: false
skip_conditions: []
default_tool_or_algorithm: DESeq2 prefiltering logic
fallback_tool_or_algorithm: edgeR filterByExpr
tool_selection_logic: Use DESeq2-compatible filtering by default when DESeq2 is downstream; use edgeR filterByExpr when edgeR is selected.
user_settable_parameters:
  - minimum_count
  - minimum_samples
expert_parameters:
  - filter_strategy
internal_parameters:
  - design_formula
qc_metrics:
  - genes_before_filtering
  - genes_after_filtering
  - library_size_per_sample
  - zero_count_fraction
qc_pass_fail_rules:
  - warn if many samples have very low library size
  - warn if most genes are filtered out
  - fail if fewer than minimum samples remain
common_failure_modes:
  - sample sheet and count matrix sample names do not match
  - severe outlier sample
  - wrong count source
nextflow_module_name: statistics/count_matrix_qc
reusable_for_other_omics: STATISTICS_COMMON
downstream_steps:
  - rnaseq_13_sample_level_qc_pca
  - rnaseq_14_differential_expression
notes_for_guided_flowchart: Show how many genes and samples remain.
notes_for_dynamic_flowchart: Emits filtered_gene_counts.
notes_for_chat_mode: Chat should flag sample/count mismatches before DE.
```

---

## 13 Sample-Level QC And PCA

```yaml
step_id: rnaseq_13_sample_level_qc_pca
step_name: Sample-level QC and PCA
category: STATISTICS_COMMON
purpose: Identify sample outliers, batch effects, and broad structure before differential testing.
why_this_step_exists: Differential expression can be invalid if samples cluster by batch, quality, or hidden covariates.
required_inputs:
  - filtered_count_matrix
  - validated_sample_sheet
optional_inputs:
  - normalized_counts
input_formats:
  - CSV
  - TSV
required_outputs:
  - pca_plot
  - sample_distance_plot
  - sample_qc_table
optional_outputs:
  - interactive_pca_html
output_formats:
  - PNG
  - PDF
  - HTML
  - CSV
depends_on:
  - rnaseq_12_count_matrix_qc_filtering
can_be_skipped: false
skip_conditions: []
default_tool_or_algorithm: DESeq2 variance-stabilized PCA
fallback_tool_or_algorithm: plotly export for interactive product view
tool_selection_logic: Use DESeq2/R outputs for reproducibility; use Plotly only when interactive UI output is needed.
user_settable_parameters:
  - color_by
  - shape_by
expert_parameters:
  - transformation_method
internal_parameters:
  - plot_theme
qc_metrics:
  - percent_variance_pc1
  - percent_variance_pc2
  - sample_distance_outliers
qc_pass_fail_rules:
  - warn if samples cluster by batch more strongly than condition
  - warn if obvious outlier sample is detected
common_failure_modes:
  - confounded design
  - mislabeled sample
  - hidden batch effect
nextflow_module_name: statistics/sample_qc_pca
reusable_for_other_omics: STATISTICS_COMMON
downstream_steps:
  - rnaseq_14_differential_expression
notes_for_guided_flowchart: Show this as a review checkpoint.
notes_for_dynamic_flowchart: Can be inserted before most matrix-based statistics.
notes_for_chat_mode: Chat should summarize possible outliers but not remove samples automatically.
```

---

## 14 Differential Expression

```yaml
step_id: rnaseq_14_differential_expression
step_name: Differential expression testing
category: STATISTICS_COMMON
purpose: Identify genes with statistically significant expression differences between conditions.
why_this_step_exists: This is the primary statistical analysis for many bulk RNA-seq studies.
required_inputs:
  - filtered_count_matrix
  - validated_sample_sheet
  - contrast_definition
optional_inputs:
  - batch_column
input_formats:
  - CSV
  - TSV
  - JSON
required_outputs:
  - differential_expression_table
  - normalized_counts
  - model_summary
optional_outputs:
  - shrinkage_results
output_formats:
  - CSV
  - TSV
  - RDS
  - JSON
depends_on:
  - rnaseq_12_count_matrix_qc_filtering
  - rnaseq_13_sample_level_qc_pca
can_be_skipped: true
skip_conditions:
  - user only wants quantification/count matrix
default_tool_or_algorithm: DESeq2
fallback_tool_or_algorithm: edgeR
tool_selection_logic: DESeq2 is default for standard bulk RNA-seq because it is widely adopted and has strong shrinkage/modeling support. edgeR is fallback for specific designs or edgeR-preferred statistical workflows.
user_settable_parameters:
  - condition_column
  - reference_condition
  - contrast
  - batch_column
  - fdr_threshold
expert_parameters:
  - design_formula
  - lfc_shrinkage
  - independent_filtering
internal_parameters:
  - model_output_prefix
qc_metrics:
  - genes_tested
  - significant_genes
  - dispersion_fit_status
  - samples_in_model
qc_pass_fail_rules:
  - fail if condition column is missing
  - fail if contrast cannot be constructed
  - warn if fewer than 3 biological replicates per condition
  - warn if design matrix is not full rank
common_failure_modes:
  - batch confounded with condition
  - too few replicates
  - sample names mismatch
  - using normalized non-count data as input
nextflow_module_name: statistics/deseq2
reusable_for_other_omics: STATISTICS_COMMON
downstream_steps:
  - rnaseq_15_ora_enrichment
  - rnaseq_16_gsea_enrichment
  - rnaseq_17_visualization_volcano_heatmap
notes_for_guided_flowchart: Ask for comparison in biological language.
notes_for_dynamic_flowchart: Requires count matrix and validated sample metadata.
notes_for_chat_mode: Chat may propose design formula but must ask for confirmation before run.
```

---

## 15 ORA Enrichment

```yaml
step_id: rnaseq_15_ora_enrichment
step_name: Over-representation enrichment analysis
category: STATISTICS_COMMON
purpose: Test whether significant genes are enriched for pathways or gene sets.
why_this_step_exists: Converts long DE gene lists into interpretable biological themes.
required_inputs:
  - differential_expression_table
  - gene_set_database
optional_inputs:
  - gene_universe
input_formats:
  - CSV
  - TSV
  - GMT
required_outputs:
  - ora_enrichment_table
  - ora_plot
optional_outputs:
  - leading_gene_list
output_formats:
  - CSV
  - PNG
  - PDF
depends_on:
  - rnaseq_14_differential_expression
can_be_skipped: true
skip_conditions:
  - user only wants DE table
default_tool_or_algorithm: clusterProfiler ORA
fallback_tool_or_algorithm: g:Profiler API/offline export
tool_selection_logic: Use clusterProfiler by default for local reproducible ORA. Use g:Profiler only when the platform supports its API/offline export policy.
user_settable_parameters:
  - gene_set_collection
  - fdr_threshold
  - log2fc_threshold
expert_parameters:
  - universe_definition
  - min_gene_set_size
  - max_gene_set_size
internal_parameters:
  - gene_id_mapping_strategy
qc_metrics:
  - input_gene_count
  - mapped_gene_count
  - enriched_term_count
qc_pass_fail_rules:
  - warn if many genes fail ID mapping
  - warn if significant gene list is too small
common_failure_modes:
  - wrong gene identifier namespace
  - no explicit universe/background
  - outdated pathway database
nextflow_module_name: statistics/clusterprofiler_ora
reusable_for_other_omics: STATISTICS_COMMON
downstream_steps:
  - rnaseq_19_report_rendering
notes_for_guided_flowchart: Explain ORA as "pathways overrepresented among significant genes".
notes_for_dynamic_flowchart: Separate from GSEA because inputs and assumptions differ.
notes_for_chat_mode: Chat should choose ORA when user asks about significant gene lists.
```

---

## 16 GSEA Enrichment

```yaml
step_id: rnaseq_16_gsea_enrichment
step_name: Gene set enrichment analysis
category: STATISTICS_COMMON
purpose: Test whether ranked genes show coordinated pathway-level shifts.
why_this_step_exists: GSEA can detect pathway signal even when individual genes do not pass a hard DE threshold.
required_inputs:
  - ranked_gene_table
  - gene_set_database
optional_inputs: []
input_formats:
  - CSV
  - TSV
  - GMT
required_outputs:
  - gsea_enrichment_table
  - enrichment_plots
optional_outputs:
  - leading_edge_genes
output_formats:
  - CSV
  - PNG
  - PDF
depends_on:
  - rnaseq_14_differential_expression
can_be_skipped: true
skip_conditions:
  - user only wants DE table
default_tool_or_algorithm: fgsea
fallback_tool_or_algorithm: clusterProfiler GSEA
tool_selection_logic: fgsea is default for fast preranked GSEA. clusterProfiler GSEA is fallback for teams already using clusterProfiler outputs.
user_settable_parameters:
  - ranking_metric
  - gene_set_collection
  - fdr_threshold
expert_parameters:
  - min_gene_set_size
  - max_gene_set_size
  - permutation_strategy
internal_parameters:
  - gene_id_mapping_strategy
qc_metrics:
  - ranked_genes
  - mapped_ranked_genes
  - significant_gene_sets
qc_pass_fail_rules:
  - warn if too few ranked genes map to gene sets
  - fail if ranking metric cannot be computed
common_failure_modes:
  - duplicate gene IDs
  - wrong ranking direction
  - mismatched species database
nextflow_module_name: statistics/fgsea
reusable_for_other_omics: STATISTICS_COMMON
downstream_steps:
  - rnaseq_19_report_rendering
notes_for_guided_flowchart: Explain as "pathway shifts across all ranked genes".
notes_for_dynamic_flowchart: Do not merge with ORA; keep as separate analytical route.
notes_for_chat_mode: Chat should choose GSEA when user asks about pathway-level trends.
```

---

## 17 Visualization

```yaml
step_id: rnaseq_17_visualization_volcano_heatmap
step_name: RNA-seq result visualization bundle
category: REPORTING_COMMON
purpose: Generate standard figures for interpreting differential expression results.
why_this_step_exists: Users need interpretable plots, but each plot should be task-specific rather than a generic plotting toolbox.
required_inputs:
  - differential_expression_table
  - normalized_counts
optional_inputs:
  - gene_annotations
input_formats:
  - CSV
  - TSV
required_outputs:
  - volcano_plot
  - ma_plot
  - heatmap_plot
optional_outputs:
  - top_gene_table
output_formats:
  - PNG
  - PDF
  - HTML
  - CSV
depends_on:
  - rnaseq_14_differential_expression
can_be_skipped: true
skip_conditions:
  - user only wants machine-readable tables
default_tool_or_algorithm: task-specific R plotting functions for volcano, MA, and heatmap outputs
fallback_tool_or_algorithm: Plotly for interactive UI exports
tool_selection_logic: Use standard R plotting for reproducible static reports; use Plotly only for interactive application views.
user_settable_parameters:
  - fdr_threshold
  - log2fc_threshold
  - top_n_genes
expert_parameters:
  - gene_label_strategy
  - heatmap_scaling_method
internal_parameters:
  - plot_theme
  - color_palette
qc_metrics:
  - plots_generated
  - genes_available_for_heatmap
qc_pass_fail_rules:
  - warn if too few significant genes for heatmap
  - fail if required DE table columns are missing
common_failure_modes:
  - missing adjusted p-value column
  - gene symbol mapping failure
  - too many labels causing unreadable plots
nextflow_module_name: reporting/rnaseq_visualization_bundle
reusable_for_other_omics: REPORTING_COMMON
downstream_steps:
  - rnaseq_19_report_rendering
notes_for_guided_flowchart: Display figures as interpretation outputs, not analysis decisions.
notes_for_dynamic_flowchart: Production UI should expose volcano, MA plot, and heatmap as child nodes even if the backend runs one reporting bundle.
notes_for_chat_mode: Chat can propose figures but should not alter thresholds without confirmation.
```

---

## 18 MultiQC Aggregation

```yaml
step_id: rnaseq_18_multiqc_aggregation
step_name: MultiQC aggregation
category: REPORTING_COMMON
purpose: Aggregate QC logs and reports into one run-level QC report.
why_this_step_exists: Production users need one place to review QC across all samples and tools.
required_inputs:
  - fastqc_outputs
optional_inputs:
  - trimming_reports
  - alignment_logs
  - quantification_logs
  - featurecounts_summary
  - custom_multiqc_config
input_formats:
  - HTML
  - ZIP
  - TXT
  - JSON
required_outputs:
  - multiqc_html
  - multiqc_data
optional_outputs:
  - multiqc_plots
output_formats:
  - HTML
  - TXT
  - JSON
  - PNG
depends_on:
  - rnaseq_03_raw_read_qc
can_be_skipped: false
skip_conditions: []
default_tool_or_algorithm: MultiQC
fallback_tool_or_algorithm: none
tool_selection_logic: MultiQC is the supported QC aggregation standard.
user_settable_parameters:
  - report_title
expert_parameters:
  - custom_multiqc_config
internal_parameters:
  - search_directory
  - module_order
qc_metrics:
  - modules_detected
  - samples_in_report
  - missing_qc_sections
qc_pass_fail_rules:
  - warn if expected modules are missing
  - fail if no QC files are detected
common_failure_modes:
  - tool output paths not published
  - unsupported output format
  - sample naming mismatch across tools
nextflow_module_name: common/multiqc
reusable_for_other_omics: REPORTING_COMMON
downstream_steps:
  - rnaseq_19_report_rendering
notes_for_guided_flowchart: Show as "combined QC report".
notes_for_dynamic_flowchart: Aggregates any compatible QC artifacts.
notes_for_chat_mode: Chat should summarize MultiQC warnings in plain language.
```

---

## 19 Report Rendering

```yaml
step_id: rnaseq_19_report_rendering
step_name: Final report rendering
category: REPORTING_COMMON
purpose: Produce a human-readable analysis report with methods, QC, figures, and result summaries.
why_this_step_exists: A production workflow needs a reproducible report, not only raw output files.
required_inputs:
  - multiqc_html
  - run_manifest
  - tool_versions
optional_inputs:
  - differential_expression_table
  - visualizations
  - enrichment_results
  - user_notes
input_formats:
  - HTML
  - CSV
  - TSV
  - JSON
  - PNG
  - PDF
required_outputs:
  - final_report
optional_outputs:
  - methods_section
  - downloadable_summary
output_formats:
  - HTML
  - PDF
  - DOCX
depends_on:
  - rnaseq_18_multiqc_aggregation
can_be_skipped: true
skip_conditions:
  - machine-readable run only
default_tool_or_algorithm: Quarto or R Markdown
fallback_tool_or_algorithm: none
tool_selection_logic: Use one reproducible reporting system to avoid template drift.
user_settable_parameters:
  - report_title
  - include_methods
  - include_qc
  - include_enrichment
expert_parameters:
  - report_template
internal_parameters:
  - artifact_paths
  - render_profile
qc_metrics:
  - report_sections_rendered
  - missing_artifacts
qc_pass_fail_rules:
  - fail if required report artifacts are missing
  - warn if optional sections cannot render
common_failure_modes:
  - missing figure paths
  - invalid table schema
  - rendering dependency failure
nextflow_module_name: reporting/render_report
reusable_for_other_omics: REPORTING_COMMON
downstream_steps:
  - rnaseq_20_provenance_export
notes_for_guided_flowchart: Final user-facing deliverable.
notes_for_dynamic_flowchart: Can accept outputs from many assay workflows.
notes_for_chat_mode: Chat may draft interpretation but report generation requires user confirmation.
```

---

## 20 Provenance And Artifact Export

```yaml
step_id: rnaseq_20_provenance_export
step_name: Provenance and artifact export
category: COMMON
purpose: Export final artifacts, run metadata, tool versions, parameters, and audit-ready provenance.
why_this_step_exists: Results must be reproducible, inspectable, and traceable.
required_inputs:
  - run_manifest
  - nextflow_trace
  - tool_versions
  - final_outputs
optional_inputs:
  - audit_context
input_formats:
  - JSON
  - TXT
  - CSV
  - HTML
required_outputs:
  - artifact_manifest
  - provenance_manifest
  - published_results_directory
optional_outputs:
  - ro_crate
output_formats:
  - JSON
  - directory
  - RO-Crate
depends_on:
  - rnaseq_19_report_rendering
can_be_skipped: false
skip_conditions: []
default_tool_or_algorithm: custom Survom run manifest exporter
fallback_tool_or_algorithm: RO-Crate export
tool_selection_logic: Use platform manifest for product needs; offer RO-Crate when interoperability is required.
user_settable_parameters:
  - export_format
expert_parameters:
  - include_intermediate_files
internal_parameters:
  - object_storage_prefix
  - audit_actor
  - retention_policy
qc_metrics:
  - artifacts_published
  - manifest_complete
  - tool_versions_recorded
qc_pass_fail_rules:
  - fail if final required outputs are not published
  - fail if tool versions are missing
common_failure_modes:
  - incomplete run due to failed upstream process
  - object storage permission error
  - missing software version capture
nextflow_module_name: common/provenance_export
reusable_for_other_omics: COMMON
downstream_steps: []
notes_for_guided_flowchart: Usually hidden; shown as "results saved".
notes_for_dynamic_flowchart: Terminal node for all workflows.
notes_for_chat_mode: Chat can explain where outputs are stored and what parameters were used.
```
