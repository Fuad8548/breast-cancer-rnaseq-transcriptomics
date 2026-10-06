# Project 1: GO + KEGG GSEA, pathway convergence, and candidate modules
# Run from the project root.
#
# Biological questions:
# 1) Which GO Biological Processes are enriched toward T3 or T1?
# 2) Which KEGG pathways show the same directional enrichment?
# 3) Which core/leading-edge genes are shared by GO and KEGG?
#
# Important:
# - The ranked list comes from the DESeq2 Wald statistic.
# - Positive values = T3-associated; negative values = T1-associated.
# - No individual gene passed DESeq2 FDR < 0.05 in the genome-wide analysis.
# - Candidate genes are exploratory, not validated biomarkers.

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results/candidates", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# 1. Load the saved ranked gene list
# -------------------------------------------------------------------

gene_list <- readRDS(
  "r_objects/gene_list_T3_vs_T1_Entrez.rds"
)

stopifnot(length(gene_list) > 1000)
stopifnot(!is.null(names(gene_list)))

cat("Ranked genes:", length(gene_list), "\n")
cat("Top ranks:\n")
print(head(gene_list))
cat("Bottom ranks:\n")
print(tail(gene_list))

# -------------------------------------------------------------------
# 2. GO Biological Process GSEA
# -------------------------------------------------------------------

gsea_go_bp <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE
)

go_df <- as.data.frame(gsea_go_bp)

saveRDS(
  gsea_go_bp,
  "results/enrichment/gsea_GO_BP_T3_vs_T1.rds"
)

write.csv(
  go_df,
  "results/enrichment/GO_BP_GSEA_T3_vs_T1_all_terms.csv",
  row.names = FALSE
)

cat(
  "GO BP terms:", nrow(go_df),
  "\nGO BP FDR < 0.05:",
  sum(go_df$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

# -------------------------------------------------------------------
# 3. GO redundancy reduction
# -------------------------------------------------------------------

gsea_go_simplified <- simplify(
  gsea_go_bp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

go_simple_df <- as.data.frame(gsea_go_simplified)

saveRDS(
  gsea_go_simplified,
  "results/enrichment/gsea_GO_BP_T3_vs_T1_simplified.rds"
)

write.csv(
  go_simple_df,
  "results/enrichment/GO_BP_GSEA_T3_vs_T1_simplified.csv",
  row.names = FALSE
)

cat(
  "Simplified GO terms:", nrow(go_simple_df),
  "\n"
)

# -------------------------------------------------------------------
# 4. Directional GO summaries
# -------------------------------------------------------------------

go_t3 <- go_simple_df[
  go_simple_df$NES > 0 & go_simple_df$p.adjust < 0.05,
]

go_t1 <- go_simple_df[
  go_simple_df$NES < 0 & go_simple_df$p.adjust < 0.05,
]

go_t3 <- go_t3[order(go_t3$p.adjust), ]
go_t1 <- go_t1[order(go_t1$p.adjust), ]

write.csv(
  head(go_t3, 30),
  "results/enrichment/GO_BP_GSEA_T3_top30.csv",
  row.names = FALSE
)

write.csv(
  head(go_t1, 30),
  "results/enrichment/GO_BP_GSEA_T1_top30.csv",
  row.names = FALSE
)

# -------------------------------------------------------------------
# 5. KEGG GSEA
# -------------------------------------------------------------------

gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE
)

kegg_df <- as.data.frame(gsea_kegg)

saveRDS(
  gsea_kegg,
  "results/enrichment/gsea_KEGG_T3_vs_T1.rds"
)

write.csv(
  kegg_df,
  "results/enrichment/KEGG_GSEA_T3_vs_T1_all_pathways.csv",
  row.names = FALSE
)

kegg_t3 <- kegg_df[
  kegg_df$NES > 0 & kegg_df$p.adjust < 0.05,
]

kegg_t1 <- kegg_df[
  kegg_df$NES < 0 & kegg_df$p.adjust < 0.05,
]

kegg_t3 <- kegg_t3[order(kegg_t3$p.adjust), ]
kegg_t1 <- kegg_t1[order(kegg_t1$p.adjust), ]

write.csv(
  head(kegg_t3, 30),
  "results/enrichment/KEGG_GSEA_T3_top30.csv",
  row.names = FALSE
)

write.csv(
  head(kegg_t1, 30),
  "results/enrichment/KEGG_GSEA_T1_top30.csv",
  row.names = FALSE
)

# -------------------------------------------------------------------
# 6. Helper: extract unique core-enrichment genes
# -------------------------------------------------------------------

extract_core_genes <- function(df) {
  if (!"core_enrichment" %in% names(df)) {
    stop("core_enrichment column not found.")
  }
  ids <- unlist(strsplit(df$core_enrichment, "/"))
  unique(as.character(ids))
}

# Top 10 pathways/processes for convergence
go_t3_top10 <- head(go_t3, 10)
go_t1_top10 <- head(go_t1, 10)
kegg_t3_top10 <- head(kegg_t3, 10)
kegg_t1_top10 <- head(kegg_t1, 10)

go_t3_genes <- extract_core_genes(go_t3_top10)
go_t1_genes <- extract_core_genes(go_t1_top10)
kegg_t3_genes <- extract_core_genes(kegg_t3_top10)
kegg_t1_genes <- extract_core_genes(kegg_t1_top10)

shared_t3 <- intersect(go_t3_genes, kegg_t3_genes)
shared_t1 <- intersect(go_t1_genes, kegg_t1_genes)

cat("Shared T3 GO/KEGG core genes:", length(shared_t3), "\n")
cat("Shared T1 GO/KEGG core genes:", length(shared_t1), "\n")

saveRDS(shared_t3, "results/candidates/shared_T3_GO_KEGG_genes.rds")
saveRDS(shared_t1, "results/candidates/shared_T1_GO_KEGG_genes.rds")

# -------------------------------------------------------------------
# 7. Annotate shared genes
# -------------------------------------------------------------------

annotate_entrez <- function(ids) {
  if (length(ids) == 0) {
    return(data.frame(
      ENTREZID = character(),
      SYMBOL = character(),
      GENENAME = character()
    ))
  }

  x <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ids,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
  )

  x <- x[!duplicated(x$ENTREZID), ]
  x
}

shared_t3_map <- annotate_entrez(shared_t3)
shared_t1_map <- annotate_entrez(shared_t1)

write.csv(
  shared_t3_map,
  "results/candidates/T3_GO_KEGG_shared_genes.csv",
  row.names = FALSE
)

write.csv(
  shared_t1_map,
  "results/candidates/T1_GO_KEGG_shared_genes.csv",
  row.names = FALSE
)

# -------------------------------------------------------------------
# 8. Rebuild DESeq2 results and Ensembl -> Entrez mapping
# -------------------------------------------------------------------

dds <- readRDS(
  "r_objects/dds_T1_T3_fitted.rds"
)

res <- results(
  dds,
  contrast = c("timepoint", "T3", "T1")
)

res_df <- as.data.frame(res)
res_df$ENSEMBL <- sub("\\..*$", "", rownames(res_df))

de_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(res_df$ENSEMBL),
  columns = "ENTREZID",
  keytype = "ENSEMBL"
)

de_map <- de_map[
  !is.na(de_map$ENTREZID),
]

de_map <- de_map[
  !duplicated(de_map$ENSEMBL),
]

res_df <- merge(
  res_df,
  de_map,
  by = "ENSEMBL",
  all.x = TRUE
)

# Keep one DESeq2 row per Entrez ID, choosing the strongest absolute Wald statistic.
de_entrez <- res_df[
  !is.na(res_df$ENTREZID),
]

de_entrez <- de_entrez[
  order(
    de_entrez$ENTREZID,
    -abs(de_entrez$stat)
  ),
]

de_entrez <- de_entrez[
  !duplicated(de_entrez$ENTREZID),
]

# -------------------------------------------------------------------
# 9. Build GO+KEGG convergence candidate tables
# -------------------------------------------------------------------

t3_convergence <- merge(
  shared_t3_map,
  de_entrez[, c(
    "ENTREZID",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj",
    "stat"
  )],
  by = "ENTREZID",
  all.x = TRUE
)

t3_convergence <- t3_convergence[
  t3_convergence$log2FoldChange > 0,
]

t3_convergence <- t3_convergence[
  order(-abs(t3_convergence$stat)),
]

t1_convergence <- merge(
  shared_t1_map,
  de_entrez[, c(
    "ENTREZID",
    "baseMean",
    "log2FoldChange",
    "pvalue",
    "padj",
    "stat"
  )],
  by = "ENTREZID",
  all.x = TRUE
)

t1_convergence <- t1_convergence[
  t1_convergence$log2FoldChange < 0,
]

t1_convergence <- t1_convergence[
  order(-abs(t1_convergence$stat)),
]

write.csv(
  t3_convergence,
  "results/candidates/T3_GO_KEGG_convergence_candidates.csv",
  row.names = FALSE
)

write.csv(
  t1_convergence,
  "results/candidates/T1_GO_KEGG_convergence_candidates.csv",
  row.names = FALSE
)

# -------------------------------------------------------------------
# 10. Summary table
# -------------------------------------------------------------------

summary_table <- data.frame(
  analysis = c(
    "GO BP GSEA",
    "GO BP simplified",
    "KEGG GSEA",
    "Shared T3 GO/KEGG core genes",
    "Shared T1 GO/KEGG core genes"
  ),
  value = c(
    nrow(go_df),
    nrow(go_simple_df),
    nrow(kegg_df),
    length(shared_t3),
    length(shared_t1)
  )
)

write.csv(
  summary_table,
  "results/enrichment/pathway_analysis_summary.csv",
  row.names = FALSE
)

cat("\nPathway-analysis stage complete.\n")
