# SurvOm Genomics Atomic Steps

This first genomics layer keeps RNA-seq unchanged and adds independent DNA/genomics foundation steps.

Implemented now:

- `genomics_01_input_validation`
- `genomics_02_reference_prepare`
- `genomics_03_raw_qc`
- `genomics_04_trim_fastq`
- `genomics_17_tool_inventory`
- `genomics_99_provenance_manifest`
- Tool-gated placeholders for:
  - `genomics_05_bwa_alignment`
  - `genomics_06_bam_processing`
  - `genomics_07_mark_duplicates`
  - `genomics_08_coverage_qc`
  - `genomics_09_bqsr`
  - `genomics_10_germline_variants`
  - `genomics_11_somatic_variants`
  - `genomics_12_variant_filtering`
  - `genomics_13_structural_variants`
  - `genomics_14_cnv_calling`
  - `genomics_15_variant_annotation`
  - `genomics_16_multiqc_report`

Tool-heavy steps such as BWA alignment, BAM processing, GATK variant calling, CNV, SV, VEP/SnpEff annotation, and MultiQC are visible as atomics in Dash, but they are tool-gated. They check for the required command and fail clearly if the tool/container is not available. They do not fabricate BAM/VCF/report outputs.

Use `genomics_17_tool_inventory` or the Dash workflow `Genomics chain - demo FASTQ plus all tool checks` when you want one safe test that completes and tells you which downstream tools are missing.

## Test Data

FASTQ:

```bash
/data/shared/vikash/mult-omics/testdatasets/test_data/human_chr22_genomics/fastq
```

Reference:

```bash
/data/shared/vikash/mult-omics/testdatasets/test_data/human_chr22_rnaseq/refs/chr22_with_ERCC92.fa
```

## Run All Genomics Atomic Checks

```bash
cd /data/shared/vikash/mult-omics
bash tests/run_all_atomic_tests.sh
```

Logs are written under:

```text
tests/logs/<timestamp>/
```

This script also compiles and runs the Dash-compatible full demo chain:

```text
genomics_01_input_validation -> genomics_02_reference_prepare -> genomics_03_raw_qc -> genomics_04_trim_fastq -> genomics_05_bwa_alignment -> genomics_06_bam_processing -> genomics_07_mark_duplicates -> genomics_08_coverage_qc -> genomics_09_bqsr -> genomics_10_germline_variants -> genomics_11_somatic_variants -> genomics_12_variant_filtering -> genomics_13_structural_variants -> genomics_14_cnv_calling -> genomics_15_variant_annotation -> genomics_16_multiqc_report
```

The generated test chain lives at:

```text
survom-pipelines/workflows/generated/genomics_full_demo_chain_cli_test/
```

The main outputs are under:

```text
survom-pipelines/workflows/generated/genomics_full_demo_chain_cli_test/results/
```

## Genomics Tool Environment

The Genomics tool environment is separate from the Dash/RNA-seq environment:

```bash
mamba activate survom-genomics
```

Installed or checked tools include:

```text
bwa-mem2, bwa, samtools, bcftools, mosdepth, fastp, fastqc, multiqc,
gatk, picard, cnvkit.py, delly, sniffles, snpEff
```

Dash/Nextflow can detect these tools through:

```text
/home/vikash/miniconda3/envs/survom-genomics/bin
```

To override this on another machine:

```bash
export SURVOM_GENOMICS_BIN=/path/to/genomics/env/bin
```

## Run Only Input Validation

Create a sample sheet:

```bash
cat > /tmp/genomics_samplesheet.csv <<'EOF'
sample_id,fastq_1,fastq_2,single_end
genomics_test,/data/shared/vikash/mult-omics/testdatasets/test_data/human_chr22_genomics/fastq/genomics_test_R1.fastq.gz,/data/shared/vikash/mult-omics/testdatasets/test_data/human_chr22_genomics/fastq/genomics_test_R2.fastq.gz,false
EOF
```

Run:

```bash
python3 survom-pipelines/bin/python/common/genomics_atomic.py validate-inputs \
  --samplesheet /tmp/genomics_samplesheet.csv \
  --outdir /tmp/genomics_01_validate
```

Main outputs:

```text
/tmp/genomics_01_validate/validated_genomics_manifest.csv
/tmp/genomics_01_validate/validation_report.json
```

## Run Only Reference Preparation

```bash
python3 survom-pipelines/bin/python/common/genomics_atomic.py reference-prepare \
  --genome-fasta /data/shared/vikash/mult-omics/testdatasets/test_data/human_chr22_rnaseq/refs/chr22_with_ERCC92.fa \
  --reference-name human_chr22_demo \
  --outdir /tmp/genomics_02_reference
```

Main outputs:

```text
/tmp/genomics_02_reference/reference_bundle.json
/tmp/genomics_02_reference/*.fai
/tmp/genomics_02_reference/*.dict
```

## Run Only Raw QC

```bash
python3 survom-pipelines/bin/python/common/genomics_atomic.py raw-qc \
  --manifest /tmp/genomics_01_validate/validated_genomics_manifest.csv \
  --outdir /tmp/genomics_03_raw_qc
```

Main outputs:

```text
/tmp/genomics_03_raw_qc/raw_qc_metrics.csv
/tmp/genomics_03_raw_qc/raw_qc_summary.json
```

## Run Only Trimming Contract

This lightweight implementation copies FASTQs and writes the downstream manifest contract. Full `fastp` trimming can be swapped into this atomic boundary later without changing downstream inputs.

```bash
python3 survom-pipelines/bin/python/common/genomics_atomic.py trim-fastq \
  --manifest /tmp/genomics_01_validate/validated_genomics_manifest.csv \
  --quality-cutoff 20 \
  --min-length 20 \
  --outdir /tmp/genomics_04_trim
```

Main outputs:

```text
/tmp/genomics_04_trim/trimmed_fastq_manifest.csv
/tmp/genomics_04_trim/trim_report.json
/tmp/genomics_04_trim/trimmed_fastq/
```

## Run Tool Inventory

This step checks downstream Genomics tool availability without failing if tools are missing.

```bash
python3 survom-pipelines/bin/python/common/genomics_atomic.py tool-inventory \
  --outdir /tmp/genomics_17_tool_inventory
```

Main outputs:

```text
/tmp/genomics_17_tool_inventory/tool_inventory.csv
/tmp/genomics_17_tool_inventory/tool_inventory.json
```

## Run One Safe Dash Chain

In Dash, choose:

```text
Omics type: Genomics
Workflow: Genomics chain - demo FASTQ plus all tool checks
```

If you do not upload files, the app uses the demo FASTQs from:

```text
/data/shared/vikash/mult-omics/testdatasets/test_data/human_chr22_genomics/fastq
```

Main outputs:

```text
01_input_validation/validated_genomics_manifest.csv
03_raw_qc/raw_qc_metrics.csv
04_trim/trimmed_fastq_manifest.csv
17_tool_inventory/tool_inventory.csv
99_provenance/run_manifest.json
```

## Run Full Step 1-16 Chain From Command Line

For the most practical local demo with real BWA/samtools/bcftools/mosdepth/MultiQC outputs, run:

```bash
cd /data/shared/vikash/mult-omics
bash scripts/run_genomics_demo_pipeline.sh
```

This writes a timestamped folder under:

```text
.test-runs/genomics_demo_<timestamp>/
```

Key real demo outputs include sorted BAM, BAI, flagstat/stats, coverage summaries, VCFs if variants are found, and a MultiQC report when MultiQC succeeds.

The Dash-compatible generated Nextflow chain is also available. It runs every visible atomic step and writes status/download artifacts for each step:

Compile the chain:

```bash
cd /data/shared/vikash/mult-omics
/home/vikash/miniconda3/envs/survom-preprocess/bin/python \
  survom-pipelines/scripts/compile_genomics_demo_chain.py \
  --outdir survom-pipelines/workflows/generated/genomics_full_demo_chain_cli
```

Run it:

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines/workflows/generated/genomics_full_demo_chain_cli
nextflow run main.nf -profile local -params-file params.yaml
```

For custom student data, provide a sample sheet:

```bash
cat > /tmp/my_genomics_samplesheet.csv <<'EOF'
sample_id,fastq_1,fastq_2,single_end
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,false
EOF
```

Then compile with custom inputs:

```bash
cd /data/shared/vikash/mult-omics
/home/vikash/miniconda3/envs/survom-preprocess/bin/python \
  survom-pipelines/scripts/compile_genomics_demo_chain.py \
  --samplesheet /tmp/my_genomics_samplesheet.csv \
  --genome-fasta /path/to/reference.fa \
  --outdir survom-pipelines/workflows/generated/my_genomics_chain
```

## Compile Genomics Nextflow Atomics

```bash
cd /data/shared/vikash/mult-omics
python3 survom-pipelines/scripts/smoke_genomics_workflows.py --keep /tmp/genomics_compile_smoke
```

This verifies standalone and chained generated workflows compile without running heavy external genomics tools.

## Dash Usage

1. Restart the Dash app.
2. Select `Omics type: Genomics`.
3. Upload FASTQ files or a CSV sample sheet.
4. Choose one of:
   - `Genomics 01 - Input validation`
   - `Genomics 02 - DNA reference preparation`
   - `Genomics 03 - Raw FASTQ QC`
   - `Genomics 04 - FASTQ trimming`
   - `Genomics 05-16` tool-gated downstream steps
   - `Genomics chain - demo FASTQ plus all tool checks`
   - `Genomics chain - FASTQ validation, QC, trim, provenance`
5. Run and inspect CSV/JSON outputs.

## Next Hand-Off

The next real genomics expansion should add tool-backed atomics for:

- BWA/BWA-MEM2 alignment
- samtools sort/index/flagstat/stats
- duplicate marking
- mosdepth coverage
- optional GATK BQSR
- GATK HaplotypeCaller/Mutect2
- bcftools/GATK filtering
- SV/CNV calling
- VEP/SnpEff annotation
- final MultiQC

Those should consume the manifests produced here.
