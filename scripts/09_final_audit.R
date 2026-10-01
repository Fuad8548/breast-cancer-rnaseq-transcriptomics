# ============================================================
# Project 1 — Final Reproducibility Audit
# Script 09 — Consistency and Output Audit
# ============================================================

cat("\n============================================\n")
cat("PROJECT 1 FINAL AUDIT\n")
cat("============================================\n\n")

# ------------------------------------------------------------
# 1. Required files
# ------------------------------------------------------------

required_files <- c(
    # Processed data
    "data/processed/dds_T1_T3_fitted.rds",
    "data/processed/gene_list_T3_vs_T1_Entrez.rds",

    # Differential expression
    "results/differential_expression/DESeq2_T3_vs_T1_all_genes.csv",
    "results/differential_expression/ensembl_annotation_map.csv",
    "results/differential_expression/annotation_mapping_summary.csv",
    "results/differential_expression/ranked_gene_list_for_GSEA.csv",

    # Main scripts
    "scripts/01_extract_T1_T3.sh",
    "scripts/02_qc_filter_normalize.R",
    "scripts/03_deseq2.R",
    "scripts/04_diagnose_DE.R",
    "scripts/05_annotation_and_ranking.R",
    "scripts/06_gsea_and_pathways.R",
    "scripts/07_final_candidate_validation_and_figures.R",
    "scripts/08_final_figures.R"
)

missing_files <- required_files[
    !file.exists(required_files)
]

cat("Required-file check:\n")

if (length(missing_files) == 0) {
    cat("  PASS — all required files exist.\n\n")
} else {
    cat("  WARNING — missing files:\n")
    print(missing_files)
    cat("\n")
}

# ------------------------------------------------------------
# 2. Load DESeq2 object
# ------------------------------------------------------------

dds <- readRDS(
    "data/processed/dds_T1_T3_fitted.rds"
)

cat("DESeq2 object check:\n")

cat(
    "  Genes:",
    nrow(dds),
    "\n"
)

cat(
    "  Samples:",
    ncol(dds),
    "\n"
)

# ------------------------------------------------------------
# 3. Sample-design check
# ------------------------------------------------------------

sample_info <- as.data.frame(
    colData(dds)
)

cat("\nSample design:\n")

print(
    table(
        sample_info$timepoint
    )
)

cat(
    "Patients:",
    length(unique(sample_info$patient_id)),
    "\n"
)

# ------------------------------------------------------------
# 4. Verify paired structure
# ------------------------------------------------------------

pair_table <- table(
    sample_info$patient_id,
    sample_info$timepoint
)

complete_pairs <- sum(
    rowSums(pair_table > 0) == 2 &
        rowSums(pair_table) == 2
)

cat("\nPaired-sample check:\n")

cat(
    "  Complete T1/T3 patient pairs:",
    complete_pairs,
    "\n"
)

if (complete_pairs == 24) {
    cat("  PASS — 24 complete pairs detected.\n")
} else {
    cat("  WARNING — expected 24 complete pairs.\n")
}

# ------------------------------------------------------------
# 5. Re-run primary DESeq2 contrast
# ------------------------------------------------------------

library(DESeq2)

res <- results(
    dds,
    contrast = c("timepoint", "T3", "T1")
)

cat("\nPrimary DESeq2 result:\n")

cat(
    "  Total tested genes:",
    nrow(res),
    "\n"
)

cat(
    "  padj < 0.10:",
    sum(res$padj < 0.10, na.rm = TRUE),
    "\n"
)

cat(
    "  padj < 0.05:",
    sum(res$padj < 0.05, na.rm = TRUE),
    "\n"
)

cat(
    "  Nominal p < 0.05:",
    sum(res$pvalue < 0.05, na.rm = TRUE),
    "\n"
)

# ------------------------------------------------------------
# 6. Ranking check
# ------------------------------------------------------------

gene_rank <- readRDS(
    "data/processed/gene_list_T3_vs_T1_Entrez.rds"
)

cat("\nGSEA ranking check:\n")

cat(
    "  Ranked genes:",
    length(gene_rank),
    "\n"
)

cat(
    "  Positive-ranked genes:",
    sum(gene_rank > 0),
    "\n"
)

cat(
    "  Negative-ranked genes:",
    sum(gene_rank < 0),
    "\n"
)

if (
    all(is.finite(gene_rank)) &&
        !is.null(names(gene_rank)) &&
        length(unique(names(gene_rank))) == length(gene_rank)
) {
    cat("  PASS — ranking is finite and Entrez IDs are unique.\n")
} else {
    cat("  WARNING — ranking requires inspection.\n")
}

# ------------------------------------------------------------
# 7. Candidate genes
# ------------------------------------------------------------

candidate_genes <- c(
    "KNL1",
    "CDK1",
    "CCNB1",
    "RAD51",
    "GYS1",
    "SCD",
    "INSR",
    "TNC"
)

annotated <- read.csv(
    "results/differential_expression/DESeq2_T3_vs_T1_all_genes.csv",
    stringsAsFactors = FALSE
)

candidate_check <- annotated[
    annotated$SYMBOL %in% candidate_genes,
    c(
        "SYMBOL",
        "log2FoldChange",
        "pvalue",
        "padj",
        "stat"
    )
]

cat("\nCandidate check:\n")

print(candidate_check)

# ------------------------------------------------------------
# 8. Final statement
# ------------------------------------------------------------

cat("\n============================================\n")
cat("INTERPRETATION CONSISTENCY CHECK\n")
cat("============================================\n")

cat(
    "\nThe project should consistently state:\n\n"
)

cat(
    "1. No individual genes reached genome-wide FDR significance.\n"
)

cat(
    "2. Pathway-level GSEA revealed coordinated biological enrichment.\n"
)

cat(
    "3. T3-associated biology centers strongly on proliferation,\n"
)
cat(
    "   chromosome segregation, DNA replication and genome maintenance.\n"
)

cat(
    "4. T1-associated biology includes metabolic and tissue-interaction\n"
)
cat(
    "   programs.\n"
)

cat(
    "5. Candidate genes are exploratory and are not validated biomarkers.\n"
)

cat("\n============================================\n")
cat("AUDIT COMPLETE\n")
cat("============================================\n")
