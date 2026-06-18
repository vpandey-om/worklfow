# Senior Review Notes For RNA-seq Tutorial Materials

## Bottom Line

Do not share `design/rna-seq/first_try.md` as the final student-facing or reviewer-facing document.

It is useful as a source draft, but it reads like a tool encyclopedia. A senior reviewer may question it because it exposes too many tools, mixes production workflow design with low-level command flags, and includes several claims that should be verified against official tool documentation before publication.

Use these two curated files instead:

- `docs/rnaseq_atomic_step_registry_curated.md`
- `docs/rnaseq_tool_selection_policy.md`

## Main Problems In The Draft

1. It starts directly with raw read QC and does not clearly define project setup, sample-sheet validation, run metadata, or reference validation as first-class production steps.

2. It lists too many tools per step. That is useful for teaching tool history, but weak for a production platform. The platform should expose one default and one fallback per step.

3. Reference compatibility is not emphasized enough. Genome FASTA, GTF, STAR index, Salmon index, transcriptome FASTA, and tx2gene mapping must be checked together. This is one of the easiest places to produce silently wrong RNA-seq results.

4. Some tool claims are too specific for an architecture document, for example exact speed claims, default flags, and precise option behavior. These should be marked `VERIFY_IN_OFFICIAL_DOCS` unless checked.

5. Trimming includes methylation/RRBS-specific guidance. That belongs in a methylation workflow note, not in the main RNA-seq MVP path.

6. The document mixes alignment-based RNA-seq, pseudoalignment RNA-seq, reference building, QC, DE, enrichment, and visualization without a strong route model. Students may copy the wrong combination of steps.

7. Visualization is treated as one broad plotting step. For a dynamic workflow UI, volcano, MA plot, heatmap, PCA, and enrichment plots should be separate visible tasks or child nodes.

8. MultiQC inputs must be route-aware. A Salmon-only run should not require STAR or featureCounts logs.

9. Report generation should not depend on provenance output if provenance is produced after the report. That creates a circular workflow dependency.

10. Differential expression guidance should repeatedly state that DESeq2/edgeR require raw integer counts, not TPM, FPKM, RPKM, or normalized values.

## Corrections Already Applied

The curated registry now includes:

- explicit reference-bundle validation,
- route-aware Salmon and STAR dependencies,
- `analysis_fastq_manifest` instead of ambiguous raw/trimmed FASTQ input,
- one-of dependency for count matrix QC from either tximport or featureCounts,
- route-aware MultiQC optional inputs,
- report rendering without a circular provenance dependency,
- clearer visualization-bundle wording for production UI design.

## What Is Safe To Give Students

Give students:

- the curated atomic registry,
- the tool selection policy,
- a small route diagram,
- one clean hands-on exercise using Salmon -> tximport -> DESeq2,
- one optional alignment route exercise using STAR -> featureCounts -> DESeq2.

Keep `first_try.md` as an internal reference draft only.

