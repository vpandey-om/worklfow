#!/usr/bin/env bash
set -euo pipefail

mkdir -p testdata/data/ggal

BASE="https://raw.githubusercontent.com/nextflow-io/rnaseq-nf/master/data/ggal"

curl -L "$BASE/ggal_gut_1.fq" -o testdata/data/ggal/ggal_gut_1.fq
curl -L "$BASE/ggal_gut_2.fq" -o testdata/data/ggal/ggal_gut_2.fq

cat > testdata/ggal_gut_samplesheet.csv <<'CSV'
sample_id,fastq_1,fastq_2,library_layout,condition,organism
ggal_gut,testdata/data/ggal/ggal_gut_1.fq,testdata/data/ggal/ggal_gut_2.fq,paired,gut,Gallus_gallus
CSV

echo "Downloaded small paired-end chicken RNA-seq test data."
echo "Sample sheet: testdata/ggal_gut_samplesheet.csv"
