#!/usr/bin/env bash
set -euo pipefail

# Usage:
# bash make_mini_gallus_ref.sh
#
# Optional:
# bash make_mini_gallus_ref.sh 1000000
#
# Argument is minimum contig length. Default = 1 Mb.

MINLEN="${1:-1000000}"

FASTA="Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz"
GTF="Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz"
OUTDIR="mini_ref"

mkdir -p "$OUTDIR"

echo "Input FASTA: $FASTA"
echo "Input GTF:   $GTF"
echo "Output dir:  $OUTDIR"
echo "Min length:  $MINLEN"

if [[ ! -f "$FASTA" ]]; then
  echo "ERROR: FASTA not found: $FASTA"
  exit 1
fi

if [[ ! -f "$GTF" ]]; then
  echo "ERROR: GTF not found: $GTF"
  exit 1
fi

echo "Finding contigs with annotation..."
zcat "$GTF" \
  | awk '$0 !~ /^#/ {print $1}' \
  | sort -u \
  > "$OUTDIR/gtf_contigs.txt"

echo "Calculating FASTA contig lengths..."
zcat "$FASTA" \
  | awk '
      /^>/ {
        if (name != "") print name "\t" len
        name=$1
        sub(/^>/, "", name)
        len=0
        next
      }
      {
        len += length($0)
      }
      END {
        if (name != "") print name "\t" len
      }
    ' \
  > "$OUTDIR/fasta_lengths.tsv"

echo "Choosing smallest annotated contig >= ${MINLEN} bp..."
awk -v m="$MINLEN" '
  NR==FNR {
    annot[$1]=1
    next
  }
  ($1 in annot) && $2 >= m {
    print $1 "\t" $2
  }
' "$OUTDIR/gtf_contigs.txt" "$OUTDIR/fasta_lengths.tsv" \
  | sort -k2,2n \
  > "$OUTDIR/annotated_contigs_by_size.tsv"

if [[ ! -s "$OUTDIR/annotated_contigs_by_size.tsv" ]]; then
  echo "ERROR: No annotated contig found with length >= $MINLEN"
  echo "Try smaller value, e.g.: bash make_mini_gallus_ref.sh 100000"
  exit 1
fi

CONTIG=$(head -1 "$OUTDIR/annotated_contigs_by_size.tsv" | cut -f1)
LEN=$(head -1 "$OUTDIR/annotated_contigs_by_size.tsv" | cut -f2)

echo "Selected contig: $CONTIG"
echo "Length:          $LEN bp"

echo "Writing mini FASTA..."
zcat "$FASTA" \
  | awk -v c="$CONTIG" '
      /^>/ {
        name=$1
        sub(/^>/, "", name)
        keep=(name==c)
      }
      keep
    ' \
  | gzip -c \
  > "$OUTDIR/Gallus_gallus.mini.fa.gz"

echo "Writing mini GTF..."
zcat "$GTF" \
  | awk -v c="$CONTIG" '
      /^#/ || $1==c
    ' \
  | gzip -c \
  > "$OUTDIR/Gallus_gallus.mini.gtf.gz"

echo
echo "Done."
echo "Mini reference files:"
ls -lh "$OUTDIR/Gallus_gallus.mini.fa.gz" "$OUTDIR/Gallus_gallus.mini.gtf.gz"

echo
echo "Check FASTA contig:"
zgrep '^>' "$OUTDIR/Gallus_gallus.mini.fa.gz" | head -1

echo
echo "Check GTF contig:"
zcat "$OUTDIR/Gallus_gallus.mini.gtf.gz" \
  | awk '$0 !~ /^#/ {print $1; exit}'

echo
echo "Use these in nf-core/rnaseq:"
echo "--fasta $(pwd)/$OUTDIR/Gallus_gallus.mini.fa.gz \\"
echo "--gtf   $(pwd)/$OUTDIR/Gallus_gallus.mini.gtf.gz"