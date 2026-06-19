# Command Line Workflow Examples

Use this when Dash/API is not working, or when you want to test a workflow with real files from the terminal.

Start from the repo root:

```bash
cd /data/shared/vikash/mult-omics
```

## 1. Basic Local Run

Input folder:

```text
/data/my_test/
  sample_R1.fastq.gz
  sample_R2.fastq.gz
```

Run trimming:

```bash
./run_workflow_local.sh \
  --workflow trim_only \
  --input /data/my_test \
  --output /tmp/survom_trim_output
```

The runner creates:

```text
runs/<run_id>/
  state.json
  inputs/
  outputs/
  logs/
```

## 2. Run With Parameter Overrides

Use `--param key=value`. You can repeat it.

Example trimming with custom thresholds:

```bash
./run_workflow_local.sh \
  --workflow trim_only \
  --input /data/my_test \
  --output /tmp/survom_trim_q25 \
  --param quality_threshold=25 \
  --param minimum_read_length=30 \
  --param trim_poly_g=auto
```

Example Cutadapt fallback:

```bash
./run_workflow_local.sh \
  --workflow trim_only \
  --input /data/my_test \
  --output /tmp/survom_cutadapt_test \
  --param trimming_tool=cutadapt \
  --param quality_threshold=20 \
  --param minimum_read_length=20
```

## 3. Run With A Manifest

Instead of raw FASTQ discovery, put a manifest in the input folder.

Raw manifest:

```text
/data/my_test/raw_manifest.csv
```

Example `raw_manifest.csv`:

```csv
sample_id,fastq_1,fastq_2,single_end,strandedness
sample1,/data/my_test/sample1_R1.fastq.gz,/data/my_test/sample1_R2.fastq.gz,false,unknown
```

Run:

```bash
./run_workflow_local.sh \
  --workflow qc_trim_strandedness \
  --input /data/my_test \
  --output /tmp/survom_qc_trim_strandedness
```

## 4. Run From A Trimmed Manifest

For post-trim QC, strandedness, Salmon quantification, STAR alignment, or HISAT2 alignment, use a trimmed manifest.

Example:

```text
/data/trimmed_test/trim_manifest.tsv
```

Example `trim_manifest.tsv`:

```tsv
sample_id	fastq_1	fastq_2	single_end	strandedness
sample1	/data/trimmed_test/sample1_R1.trimmed.fastq.gz	/data/trimmed_test/sample1_R2.trimmed.fastq.gz	false	unknown
```

Run post-trim QC:

```bash
./run_workflow_local.sh \
  --workflow post_trim_qc_only \
  --input /data/trimmed_test \
  --output /tmp/survom_post_trim_qc
```

Run strandedness:

```bash
./run_workflow_local.sh \
  --workflow strandedness_only \
  --input /data/trimmed_test \
  --output /tmp/survom_strandedness
```

## 5. Run With Custom Reference Or Index Paths

Use `--param` to provide reference/index paths.

Salmon route:

```bash
./run_workflow_local.sh \
  --workflow salmon_quant_only \
  --input /data/trimmed_test \
  --output /tmp/survom_salmon_quant \
  --param selected_route=salmon \
  --param salmon_index=/data/refs/my_salmon_index
```

STAR route:

```bash
./run_workflow_local.sh \
  --workflow star_align_only \
  --input /data/trimmed_test \
  --output /tmp/survom_star_align \
  --param selected_route=star \
  --param star_index=/data/refs/my_star_index \
  --param gtf=/data/refs/genes.gtf
```

HISAT2 route:

```bash
./run_workflow_local.sh \
  --workflow hisat2_align_only \
  --input /data/trimmed_test \
  --output /tmp/survom_hisat2_align \
  --param selected_route=hisat2 \
  --param hisat2_index=/data/refs/hisat2/genome
```

If an index path is missing or invalid, the compiler should fail early with a clear message.

## 6. Run Count Matrix Routes

Salmon transcript count route:

```bash
./run_workflow_local.sh \
  --workflow salmon_count_matrix \
  --input /data/trimmed_test \
  --output /tmp/survom_salmon_count_matrix \
  --param salmon_index=/data/refs/my_salmon_index
```

STAR alignment and gene counting route:

```bash
./run_workflow_local.sh \
  --workflow star_count_matrix \
  --input /data/trimmed_test \
  --output /tmp/survom_star_count_matrix \
  --param star_index=/data/refs/my_star_index \
  --param gtf=/data/refs/genes.gtf
```

## 7. Resume

The command prints a run ID. Resume with:

```bash
./run_workflow_local.sh --resume <run_id>
```

Example:

```bash
./run_workflow_local.sh --resume trim_only_20260619_120000
```

Resume reads:

```text
runs/<run_id>/state.json
```

and relaunches Nextflow with `-resume`.

## 8. Use A Different Nextflow Profile

Default:

```text
local,docker
```

Override:

```bash
./run_workflow_local.sh \
  --workflow trim_only \
  --input /data/my_test \
  --output /tmp/survom_trim_conda \
  --profile local,conda
```

Only use a profile if it is configured in the generated Nextflow config and the environment supports it.

## 9. Check Logs

For a local run:

```bash
cat runs/<run_id>/state.json
tail -n 100 runs/<run_id>/logs/nextflow.stdout.log
tail -n 100 runs/<run_id>/logs/nextflow.stderr.log
tail -n 100 runs/<run_id>/logs/.nextflow.log
```

For all-test runs:

```bash
ls -ltr .test-runs
tail -n 100 .test-runs/$(ls -1 .test-runs | tail -1)/workflow_real_Nextflow_smoke.log
```

## 10. Good Real-Test Command

Use this as a simple real command-line test:

```bash
./run_workflow_local.sh \
  --workflow qc_trim_strandedness \
  --input /data/shared/vikash/mult-omics/survom-pipelines/testdata \
  --output /tmp/survom_real_cli_test
```

If your input folder does not contain FASTQ files or a valid manifest, the runner will stop with a clear input error.
