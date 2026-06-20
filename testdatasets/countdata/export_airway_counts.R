#!/usr/bin/env Rscript

dir.create("test_data/airway", recursive = TRUE, showWarnings = FALSE)

library(airway)
library(SummarizedExperiment)

data("airway")

counts <- assay(airway)

write.csv(
  data.frame(gene_id = rownames(counts), counts, check.names = FALSE),
  "test_data/airway/airway_raw_counts.csv",
  row.names = FALSE,
  quote = FALSE
)

meta <- as.data.frame(colData(airway))
meta$sample_id <- rownames(meta)
meta$condition <- as.character(meta$dex)
meta$pair_id <- as.character(meta$cell)

meta <- meta[match(colnames(counts), meta$sample_id),
             c("sample_id", "condition", "pair_id", "cell", "dex"),
             drop = FALSE]

write.csv(
  meta,
  "test_data/airway/airway_sample_metadata.csv",
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(
    contrast_id = "treated_vs_untreated",
    numerator = "trt",
    denominator = "untrt",
    design_column = "condition"
  ),
  "test_data/airway/airway_contrasts.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("Done: test_data/airway/\n")