# This script takes the fitted DESeq2 result, maps Ensembl IDs to human gene identifiers, prepares a signed Wald-statistic ranking for GSEA, and saves the mapping/ranking tables.

# ============================================================
# Project 1 — Paired T1 vs T3 Transcriptomics
# Script 05 — Annotation and Gene Ranking
# ============================================================

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------
# 2. Load fitted DESeq2 object
# ------------------------------------------------------------

dds_file <- "r_objects/dds_T1_T3_fitted.rds"

if (!file.exists(dds_file)) {
  stop(
    "Could not find fitted DESeq2 object: ",
    dds_file,
    "\nRun 03_deseq2.R and 04_diagnose_DE.R first."
  )
}

dds <- readRDS(dds_file)

# ------------------------------------------------------------
# 3. Extract DESeq2 results
# ------------------------------------------------------------

res <- results(
  dds,
  contrast = c("timepoint", "T3", "T1")
)

res_df <- as.data.frame(res)

# Preserve the original Ensembl identifier.
res_df$ENSEMBL <- rownames(res_df)

# ------------------------------------------------------------
# 4. Remove Ensembl version suffixes
#
# Example:
# ENSG00000123456.8
# becomes
# ENSG00000123456
# ------------------------------------------------------------

res_df$ENSEMBL <- sub(
  "\\..*$",
  "",
  res_df$ENSEMBL
)

# ------------------------------------------------------------
# 5. Check basic dimensions
# ------------------------------------------------------------

cat("\n============================================\n")
cat("Annotation and ranking summary\n")
cat("============================================\n")

cat("Number of DESeq2 genes:", nrow(res_df), "\n")
cat(
  "Unique Ensembl IDs:",
  length(unique(res_df$ENSEMBL)),
  "\n"
)

# ------------------------------------------------------------
# 6. Annotate Ensembl IDs
#
# SYMBOL    = HGNC gene symbol
# ENTREZID  = NCBI Entrez Gene ID
# GENENAME  = full gene name
# ------------------------------------------------------------

annotation_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(res_df$ENSEMBL),
  keytype = "ENSEMBL",
  columns = c(
    "ENSEMBL",
    "SYMBOL",
    "ENTREZID",
    "GENENAME"
  )
)

# Remove exact duplicate annotation rows.
annotation_raw <- unique(annotation_raw)

# ------------------------------------------------------------
# 7. Handle multiple annotation records
#
# A single Ensembl ID can sometimes map to multiple records.
# For the main annotation table we retain one record per Ensembl ID.
#
# Priority:
#   1. record with Entrez ID
#   2. record with gene symbol
#   3. otherwise first available record
# ------------------------------------------------------------

annotation_map <- annotation_raw[
  !is.na(annotation_raw$ENTREZID) |
    !is.na(annotation_raw$SYMBOL),
]

annotation_map <- annotation_map[
  order(
    annotation_map$ENSEMBL,
    is.na(annotation_map$ENTREZID),
    is.na(annotation_map$SYMBOL)
  ),
]

annotation_map <- annotation_map[
  !duplicated(annotation_map$ENSEMBL),
]

# ------------------------------------------------------------
# 8. Merge annotation with DESeq2 results
# ------------------------------------------------------------

annotated_df <- merge(
  res_df,
  annotation_map,
  by = "ENSEMBL",
  all.x = TRUE,
  sort = FALSE
)

# Put ENSEMBL first for readability.
preferred_order <- c(
  "ENSEMBL",
  "SYMBOL",
  "ENTREZID",
  "GENENAME",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

preferred_order <- preferred_order[
  preferred_order %in% colnames(annotated_df)
]

remaining_columns <- setdiff(
  colnames(annotated_df),
  preferred_order
)

annotated_df <- annotated_df[
  ,
  c(preferred_order, remaining_columns)
]

# ------------------------------------------------------------
# 9. Save annotation map and annotated DESeq2 results
# ------------------------------------------------------------

dir.create(
  "results/differential_expression",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "r_objects",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  annotation_map,
  "results/differential_expression/ensembl_annotation_map.csv",
  row.names = FALSE
)

write.csv(
  annotated_df,
  "results/differential_expression/DESeq2_T3_vs_T1_annotated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Annotation mapping summary
# ------------------------------------------------------------

mapping_summary <- data.frame(
  metric = c(
    "DESeq2 genes",
    "Unique Ensembl IDs",
    "Ensembl IDs with SYMBOL",
    "Ensembl IDs with ENTREZID",
    "Ensembl IDs with GENENAME"
  ),
  value = c(
    nrow(res_df),
    length(unique(res_df$ENSEMBL)),
    sum(!is.na(annotated_df$SYMBOL)),
    sum(!is.na(annotated_df$ENTREZID)),
    sum(!is.na(annotated_df$GENENAME))
  )
)

write.csv(
  mapping_summary,
  "results/differential_expression/annotation_mapping_summary.csv",
  row.names = FALSE
)

print(mapping_summary)

# ------------------------------------------------------------
# 11. Prepare ranked gene list for GSEA
#
# Ranking variable:
#
#   signed DESeq2 Wald statistic
#
# Positive values = T3-associated
# Negative values = T1-associated
#
# We deliberately do NOT filter genes using p < 0.05 here.
# GSEA should use the complete ranked signal.
# ------------------------------------------------------------

ranking_df <- annotated_df[
  !is.na(annotated_df$ENTREZID) &
    !is.na(annotated_df$stat),
]

ranking_df$ENTREZID <- as.character(
  ranking_df$ENTREZID
)

# ------------------------------------------------------------
# 12. Remove duplicated Entrez IDs
#
# An Entrez gene can occasionally correspond to multiple
# Ensembl entries.
#
# Keep the record with the largest absolute Wald statistic.
# ------------------------------------------------------------

ranking_df <- ranking_df[
  order(
    abs(ranking_df$stat),
    decreasing = TRUE
  ),
]

ranking_unique <- ranking_df[
  !duplicated(ranking_df$ENTREZID),
]

# ------------------------------------------------------------
# 13. Create named numeric ranking vector
# ------------------------------------------------------------

gene_rank <- ranking_unique$stat

names(gene_rank) <- ranking_unique$ENTREZID

# Remove non-finite values.
gene_rank <- gene_rank[
  is.finite(gene_rank)
]

# Sort from strongest T3-associated signal
# to strongest T1-associated signal.
gene_rank <- sort(
  gene_rank,
  decreasing = TRUE
)

# ------------------------------------------------------------
# 14. Save ranking table
# ------------------------------------------------------------

ranked_gene_table <- data.frame(
  ENTREZID = names(gene_rank),
  Wald_statistic = as.numeric(gene_rank)
)

# Add SYMBOL where possible.
symbol_lookup <- ranking_unique[
  match(
    ranked_gene_table$ENTREZID,
    ranking_unique$ENTREZID
  ),
  "SYMBOL"
]

ranked_gene_table$SYMBOL <- symbol_lookup

ranked_gene_table <- ranked_gene_table[
  ,
  c("ENTREZID", "SYMBOL", "Wald_statistic")
]

write.csv(
  ranked_gene_table,
  "results/differential_expression/ranked_gene_list_for_GSEA.csv",
  row.names = FALSE
)

saveRDS(
  gene_rank,
  "r_objects/gene_list_T3_vs_T1_Entrez.rds"
)

# ------------------------------------------------------------
# 15. Final summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("Final annotation statistics\n")
cat("============================================\n")

cat(
  "Genes with SYMBOL:",
  sum(!is.na(annotated_df$SYMBOL)),
  "\n"
)

cat(
  "Genes with ENTREZID:",
  sum(!is.na(annotated_df$ENTREZID)),
  "\n"
)

cat(
  "Unique Entrez genes in ranking:",
  length(gene_rank),
  "\n"
)

cat(
  "Highest Wald statistic:",
  max(gene_rank),
  "\n"
)

cat(
  "Lowest Wald statistic:",
  min(gene_rank),
  "\n"
)

cat("\nFiles written:\n")
cat(
  "  results/differential_expression/ensembl_annotation_map.csv\n"
)
cat(
  "  results/differential_expression/DESeq2_T3_vs_T1_annotated.csv\n"
)
cat(
  "  results/differential_expression/annotation_mapping_summary.csv\n"
)
cat(
  "  results/differential_expression/ranked_gene_list_for_GSEA.csv\n"
)
cat(
  "  r_objects/gene_list_T3_vs_T1_Entrez.rds\n"
)

cat("\nScript 05 completed successfully.\n")
