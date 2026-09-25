# Project 1: Diagnose DESeq2 T3 vs T1 results
# Run from the project root after creating/fitting dds.

suppressPackageStartupMessages(library(DESeq2))

counts_filtered <- read.delim(
  "data/processed/counts_T1_T3_filtered.tsv",
  row.names = 1,
  check.names = FALSE
)
metadata <- read.csv(
  "data/processed/metadata_T1_T3.csv",
  row.names = 1,
  check.names = FALSE
)
metadata$patient_id <- factor(metadata$patient_id)
metadata$timepoint <- factor(metadata$timepoint, levels = c("T1", "T3"))

stopifnot(all(colnames(counts_filtered) == rownames(metadata)))

int_counts_fil = round(counts_filtered)

dds <- DESeqDataSetFromMatrix(
  countData = int_counts_fil,
  colData = metadata,
  design = ~ patient_id + timepoint
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("timepoint", "T3", "T1"))

# Create output directory
dir.create("results/differential_expression", recursive = TRUE, showWarnings = FALSE)

# Basic diagnostics
n_p05 <- sum(res$pvalue < 0.05, na.rm = TRUE)
n_p01 <- sum(res$pvalue < 0.01, na.rm = TRUE)
n_p001 <- sum(res$pvalue < 0.001, na.rm = TRUE)
max_abs_lfc <- max(abs(res$log2FoldChange), na.rm = TRUE)

cat("Genes with nominal p < 0.05:", n_p05, "\n")
cat("Genes with nominal p < 0.01:", n_p01, "\n")
cat("Genes with nominal p < 0.001:", n_p001, "\n")
cat("Maximum absolute log2FC:", max_abs_lfc, "\n")

# Top genes by absolute effect size
res_lfc <- res[order(abs(res$log2FoldChange), decreasing = TRUE), ]
write.csv(
  as.data.frame(res_lfc),
  "results/differential_expression/top_by_abs_log2FC.csv",
  row.names = TRUE
)

# MA plot
pdf("results/differential_expression/MA_plot_T3_vs_T1.pdf", width = 8, height = 6)
plotMA(res, alpha = 0.05, main = "MA Plot: T3 vs T1")
dev.off()

# Save complete result
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
write.csv(
  res_df,
  "results/differential_expression/DESeq2_T3_vs_T1_all_genes.csv",
  row.names = FALSE
)

# Save fitted object for later stages
saveRDS(dds, "data/processed/dds_T1_T3_fitted.rds")
