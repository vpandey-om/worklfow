# Example Assets

`samplesheet.csv` is a header-only placeholder used by the example run requests.

For paired-end data, use:

```csv
sample_id,library_layout,fastq_1,fastq_2,platform
sample_a,paired,/path/sample_a_R1.fastq.gz,/path/sample_a_R2.fastq.gz,illumina
```

For single-end data, leave `fastq_2` empty:

```csv
sample_id,library_layout,fastq_1,fastq_2,platform
sample_b,single,/path/sample_b.fastq.gz,,illumina
```
