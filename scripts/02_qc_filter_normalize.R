# Project 1: Breast-cancer RNA-seq transcriptomics
# Stage 02: QC, low-expression filtering, DESeq2 object, normalization, VST, PCA.

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))

# Paths
count_file <- "data/processed/counts_T1_T3.tsv"
metadata_file <- "data/processed/metadata_T1_T3.csv"

qc_dir <- "results/qc"
diff_dir <- "results/differential_expression"

for (d in c(qc_dir, diff_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Load data
counts <- read.delim(
  count_file,
  row.names = 1,
  check.names = FALSE
)

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE
)

# Experimental variables
metadata$patient_id <- factor(metadata$patient_id)
metadata$timepoint <- factor(metadata$timepoint, levels = c("T1", "T3"))

# Validate structure
stopifnot(ncol(counts) == nrow(metadata))
stopifnot(all(colnames(counts) == rownames(metadata)))
stopifnot(all(counts >= 0))
stopifnot(all(counts == round(as.matrix(counts))))
stopifnot(identical(as.integer(table(metadata$timepoint)), c(24L, 24L)))

# -------------------------
# Basic sample-level QC
# -------------------------

library_sizes <- colSums(counts)
write.table(
  data.frame(sample_id = names(library_sizes), library_size = as.numeric(library_sizes)),
  file = file.path(qc_dir, "library_sizes.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

pdf(file.path(qc_dir, "library_sizes.pdf"), width = 12, height = 6)
barplot(
  library_sizes / 1e6,
  las = 2,
  cex.names = 0.6,
  ylab = "Total counts (millions)",
  main = "RNA-seq library size"
)
dev.off()

detected_genes <- colSums(counts > 0)
write.table(
  data.frame(sample_id = names(detected_genes), detected_genes = as.numeric(detected_genes)),
  file = file.path(qc_dir, "detected_genes.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

pdf(file.path(qc_dir, "detected_genes.pdf"), width = 12, height = 6)
barplot(
  detected_genes / 1000,
  las = 2,
  cex.names = 0.6,
  ylab = "Detected genes (thousands)",
  main = "Genes detected per sample"
)
dev.off()

# -------------------------
# Low-expression filtering
# Keep genes with >=10 counts in >=24 of 48 samples.
# -------------------------

keep <- rowSums(counts >= 10) >= 24
counts_filtered <- counts[keep, , drop = FALSE]

write.table(
  counts_filtered,
  file = "data/processed/counts_T1_T3_filtered.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.csv(
  data.frame(
    total_genes = nrow(counts),
    retained_genes = sum(keep),
    filtered_genes = sum(!keep),
    percent_retained = 100 * sum(keep) / nrow(counts)
  ),
  file = file.path(qc_dir, "filtering_summary.csv"),
  row.names = FALSE
)

# -------------------------
# DESeq2 object
# -------------------------
int_counts_fil = round(counts_filtered)

dds <- DESeqDataSetFromMatrix(
  countData = int_counts_fil,
  colData = metadata,
  design = ~ patient_id + timepoint
)

# -------------------------
# Normalization
# -------------------------

dds <- estimateSizeFactors(dds)
size_factor_table <- data.frame(
  sample_id = names(sizeFactors(dds)),
  size_factor = as.numeric(sizeFactors(dds))
)
write.csv(
  size_factor_table,
  file = file.path(qc_dir, "normalization_size_factors.csv"),
  row.names = FALSE
)

pdf(file.path(qc_dir, "size_factors.pdf"), width = 12, height = 6)
barplot(
  sizeFactors(dds),
  names.arg = names(sizeFactors(dds)),
  las = 2,
  cex.names = 0.6,
  ylab = "DESeq2 size factor",
  main = "DESeq2 normalization size factors"
)
dev.off()

# Normalized counts for inspection/visualization.
normalized_counts <- counts(dds, normalized = TRUE)
write.table(
  normalized_counts,
  file = "data/processed/counts_T1_T3_normalized.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Variance-stabilizing transformation for exploratory analyses.
vsd <- vst(dds, blind = FALSE)

# -------------------------
# PCA
# -------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = "timepoint",
  returnData = TRUE
)
pca_data$patient_id <- metadata[rownames(pca_data), "patient_id"]
percent_var <- attr(pca_data, "percentVar")

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of Breast Cancer RNA-seq Samples",
    x = paste0("PC1: ", round(100 * percent_var[1]), "% variance"),
    y = paste0("PC2: ", round(100 * percent_var[2]), "% variance"),
    color = "Timepoint"
  ) +
  theme_minimal()

ggsave(
  file.path(qc_dir, "PCA_T1_T3.pdf"),
  p, width = 8, height = 6
)

p_pair <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = group, group = patient_id)
) +
  geom_line(color = "grey70", linewidth = 0.6) +
  geom_point(size = 3) +
  labs(
    title = "Paired PCA: T1 to T3",
    x = paste0("PC1: ", round(100 * percent_var[1]), "% variance"),
    y = paste0("PC2: ", round(100 * percent_var[2]), "% variance"),
    color = "Timepoint"
  ) +
  theme_minimal()

ggsave(
  file.path(qc_dir, "PCA_paired_T1_T3.pdf"),
  p_pair, width = 8, height = 6
)

# Save the DESeq2 object for continuation without rerunning this stage.
saveRDS(dds, "data/processed/dds_normalized.rds")

message("Stage 02 completed successfully.")
message("Filtered genes: ", nrow(counts_filtered))
message("Samples: ", ncol(counts_filtered))
message("DESeq2 object saved to data/processed/dds_normalized.rds")
