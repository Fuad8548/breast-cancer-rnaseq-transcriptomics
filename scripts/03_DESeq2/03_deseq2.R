# Stage 03: Fit DESeq2 model and extract T3 vs T1 results.

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(tidyverse)
})


# Load the normalized/prepared DESeq2 object created in Stage 02.
dds <- readRDS("r_objects/dds_normalized.rds")

dir.create(
  "results/differential_expression",
  recursive = TRUE,
  showWarnings = FALSE
)

# Fit the negative-binomial model.
dds <- DESeq(dds)

# Inspect available model coefficients.
writeLines(
  resultsNames(dds),
  "results/differential_expression/results_names.txt"
)

# Explicitly extract T3 relative to T1.
res <- results(
  dds,
  contrast = c("timepoint", "T3", "T1")
)

# Save complete result table.
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
write.csv(
  res_df,
  "results/differential_expression/DESeq2_T3_vs_T1_all_genes.csv",
  row.names = FALSE
)

# Save sorted views for inspection.
res_by_padj <- res[order(res$padj), ]
write.csv(
  cbind(gene_id = rownames(res_by_padj), as.data.frame(res_by_padj)),
  "results/differential_expression/DESeq2_T3_vs_T1_sorted_by_padj.csv",
  row.names = FALSE
)

res_by_abs_lfc <- res[order(abs(res$log2FoldChange), decreasing = TRUE), ]
write.csv(
  cbind(gene_id = rownames(res_by_abs_lfc), as.data.frame(res_by_abs_lfc)),
  "results/differential_expression/DESeq2_T3_vs_T1_sorted_by_abs_log2FC.csv",
  row.names = FALSE
)

# Save a compact summary.
summary_lines <- capture.output(summary(res))
writeLines(
  summary_lines,
  "results/differential_expression/DESeq2_T3_vs_T1_summary.txt"
)

# Save the fitted object for downstream analyses.
saveRDS(dds, "r_objects/dds_fitted.rds")

message("Stage 03 completed successfully.")
message("T3 vs T1 results saved.")
