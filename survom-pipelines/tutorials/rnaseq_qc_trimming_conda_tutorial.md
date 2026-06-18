# RNA-seq QC and Trimming Tutorial

This tutorial runs the already-created modular SurvOm Nextflow workflow using the local conda environment:

```text
survom-preprocess
```

Do this before every local run:

```bash
cd /data/shared/vikash/mult-omics/survom-pipelines
source scripts/activate_survom_preprocess.sh
```

The workflow uses:

```text
FastQC v0.12.1
fastp 1.3.4
Nextflow 24.10.4 with system Java 17
```

## Test Data

Download the tiny paired-end chicken RNA-seq example:

```bash
bash testdata/download_nextflow_ggal_minimal.sh
```

Sample sheet:

```text
testdata/ggal_gut_samplesheet.csv
```

## Step 1 Only: Raw Read QC

Compile:

```bash
python3 -m builder.api compile \
  --request examples/rnaseq_qc_only_ggal.run.yaml \
  --output workflows/generated/ggal_qc_only_test
```

Run:

```bash
nextflow run workflows/generated/ggal_qc_only_test/main.nf \
  -profile local \
  -params-file workflows/generated/ggal_qc_only_test/params.yaml \
  -resume
```

Outputs:

```bash
find results/ggal_qc_only_test -maxdepth 4 -type f | sort
```

Expected main files:

```text
03_raw_qc/fastqc/ggal_gut_1_fastqc.html
03_raw_qc/fastqc/ggal_gut_2_fastqc.html
03_raw_qc/fastqc/ggal_gut_1_fastqc.zip
03_raw_qc/fastqc/ggal_gut_2_fastqc.zip
```

## Step 2 Only: Adapter and Quality Trimming

Compile:

```bash
python3 -m builder.api compile \
  --request examples/rnaseq_trim_only_ggal.run.yaml \
  --output workflows/generated/ggal_trim_only_test
```

Run:

```bash
nextflow run workflows/generated/ggal_trim_only_test/main.nf \
  -profile local \
  -params-file workflows/generated/ggal_trim_only_test/params.yaml \
  -resume
```

Outputs:

```bash
find results/ggal_trim_only_test -maxdepth 4 -type f | sort
```

Expected main files:

```text
04_trimmed/ggal_gut.fastp.html
04_trimmed/ggal_gut.fastp.json
04_trimmed/ggal_gut_R1.trimmed.fastq.gz
04_trimmed/ggal_gut_R2.trimmed.fastq.gz
04_trimmed/trimmed_fastq_manifest.tsv
```

## Workflow: Raw QC Plus Trimming

Compile:

```bash
python3 -m builder.api compile \
  --request examples/rnaseq_qc_trim_ggal.run.yaml \
  --output workflows/generated/ggal_qc_trim_test
```

Run:

```bash
nextflow run workflows/generated/ggal_qc_trim_test/main.nf \
  -profile local \
  -params-file workflows/generated/ggal_qc_trim_test/params.yaml \
  -resume
```

Outputs:

```bash
find results/ggal_qc_trim_test -maxdepth 4 -type f | sort
```

## Environment File

The conda environment is recorded here:

```text
envs/preprocess.yml
```

Recreate it:

```bash
mamba env create -f envs/preprocess.yml
```

If the environment already exists:

```bash
mamba env update -n survom-preprocess -f envs/preprocess.yml
```
