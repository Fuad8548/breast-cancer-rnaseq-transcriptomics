# ============================================================
# Project 1 — Paired T1 vs T3 Transcriptomics
# Script 11 — Leading-edge Gene Heatmap
#
# Purpose:
#   Identify genes repeatedly contributing to the dominant
#   T3-associated GO Biological Process enrichment and visualize
#   their expression across the matched T1/T3 samples.
#
# Important:
#   These genes are pathway-leading-edge genes, not necessarily
#   genome-wide significant DEGs.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
    library(DESeq2)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(pheatmap)
})


# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------

dds_file <- "r_objects/dds_T1_T3_fitted.rds"

rank_file <- "r_objects/gene_list_T3_vs_T1_Entrez.rds"

if (!file.exists(dds_file)) {
    stop(
        "Missing DESeq2 object: ",
        dds_file
    )
}

if (!file.exists(rank_file)) {
    stop(
        "Missing ranked gene list: ",
        rank_file
    )
}


# ------------------------------------------------------------
# 3. Output directories
# ------------------------------------------------------------

dir.create(
    "results/enrichment",
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    "results/figures",
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Load DESeq2 object and ranked gene list
# ------------------------------------------------------------

dds <- readRDS(
    dds_file
)

gene_rank <- readRDS(
    rank_file
)


# ------------------------------------------------------------
# 5. Run GO Biological Process GSEA
#
# This reproduces the ranked pathway analysis using the saved
# ranking. The purpose here is specifically to retrieve the
# core-enrichment / leading-edge genes.
# ------------------------------------------------------------

cat("\nRunning GO Biological Process GSEA...\n")

gsea_go_bp <- gseGO(
    geneList = gene_rank,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE
)


# ------------------------------------------------------------
# 6. Remove redundant GO terms
# ------------------------------------------------------------

gsea_go_bp_simplified <- simplify(
    gsea_go_bp,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
)


# ------------------------------------------------------------
# 7. Select dominant T3-associated terms
#
# We focus on positively enriched terms because the purpose of
# this figure is to explain the dominant T3 biological program.
# ------------------------------------------------------------

gsea_table <- as.data.frame(
    gsea_go_bp_simplified
)

t3_terms <- gsea_table[
    gsea_table$NES > 0 &
        gsea_table$p.adjust < 0.05,
]

t3_terms <- t3_terms[
    order(
        t3_terms$p.adjust,
        -t3_terms$NES
    ),
]

# Use the strongest 10 simplified terms.
n_terms <- min(
    10,
    nrow(t3_terms)
)

if (n_terms < 1) {
    stop(
        "No significant positive GO Biological Process terms were found."
    )
}

selected_terms <- t3_terms[
    seq_len(n_terms),
]


# ------------------------------------------------------------
# 8. Save selected terms
# ------------------------------------------------------------

write.csv(
    selected_terms,
    "results/enrichment/leading_edge_selected_T3_GO_terms.csv",
    row.names = FALSE
)


# ------------------------------------------------------------
# 9. Extract leading-edge genes
#
# `core_enrichment` contains the genes that contribute most
# strongly to each enriched gene set.
#
# IDs are stored as slash-separated Entrez IDs.
# ------------------------------------------------------------

leading_edge_list <- lapply(
    seq_len(nrow(selected_terms)),
    function(i) {
        ids <- unlist(
            strsplit(
                selected_terms$core_enrichment[i],
                split = "/"
            )
        )

        unique(ids)
    }
)

names(leading_edge_list) <-
    selected_terms$ID


# ------------------------------------------------------------
# 10. Count recurrence
#
# A gene appearing in many of the top enriched terms is more
# representative of the shared biological program than a gene
# appearing in only one term.
# ------------------------------------------------------------

all_leading_edges <- unique(
    unlist(
        leading_edge_list
    )
)

recurrence_table <- data.frame(
    ENTREZID = all_leading_edges,
    n_terms = sapply(
        all_leading_edges,
        function(gene) {
            sum(
                sapply(
                    leading_edge_list,
                    function(set) gene %in% set
                )
            )
        }
    ),
    stringsAsFactors = FALSE
)

recurrence_table <- recurrence_table[
    order(
        -recurrence_table$n_terms,
        recurrence_table$ENTREZID
    ),
]


# ------------------------------------------------------------
# 11. Annotate Entrez IDs
# ------------------------------------------------------------

annotation <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = recurrence_table$ENTREZID,
    keytype = "ENTREZID",
    columns = c(
        "ENTREZID",
        "SYMBOL",
        "ENSEMBL",
        "GENENAME"
    )
)

annotation <- annotation[
    !duplicated(
        annotation$ENTREZID
    ),
]


# ------------------------------------------------------------
# 12. Merge recurrence + annotation
# ------------------------------------------------------------

recurrence_table <- merge(
    recurrence_table,
    annotation,
    by = "ENTREZID",
    all.x = TRUE,
    sort = FALSE
)

recurrence_table <- recurrence_table[
    order(
        -recurrence_table$n_terms
    ),
]


# ------------------------------------------------------------
# 13. Save complete leading-edge table
# ------------------------------------------------------------

write.csv(
    recurrence_table,
    "results/enrichment/T3_leading_edge_gene_recurrence.csv",
    row.names = FALSE
)


# ------------------------------------------------------------
# 14. Select top recurrent genes
#
# Limit the heatmap to the 20 most recurrent genes for
# interpretability.
# ------------------------------------------------------------

candidate_leading_edge <- recurrence_table[
    !is.na(recurrence_table$ENSEMBL) &
        !is.na(recurrence_table$SYMBOL),
]

n_genes <- min(
    20,
    nrow(candidate_leading_edge)
)

if (n_genes < 5) {
    stop(
        "Too few annotated leading-edge genes available for the heatmap."
    )
}

candidate_leading_edge <-
    candidate_leading_edge[
        seq_len(n_genes),
    ]


# ------------------------------------------------------------
# 15. Load variance-stabilized expression
# ------------------------------------------------------------

vsd <- vst(
    dds,
    blind = FALSE
)

expr <- assay(
    vsd
)

rownames(expr) <- sub(
    "\\..*$",
    "",
    rownames(expr)
)


# ------------------------------------------------------------
# 16. Match selected leading-edge genes to expression matrix
# ------------------------------------------------------------

keep <- candidate_leading_edge$ENSEMBL %in%
    rownames(expr)

candidate_leading_edge <-
    candidate_leading_edge[keep, ]

if (nrow(candidate_leading_edge) < 5) {
    stop(
        "Too few leading-edge genes matched the VST expression matrix."
    )
}

leading_expr <- expr[
    candidate_leading_edge$ENSEMBL, ,
    drop = FALSE
]


# ------------------------------------------------------------
# 17. Replace Ensembl row labels with gene symbols
# ------------------------------------------------------------

rownames(leading_expr) <-
    candidate_leading_edge$SYMBOL


# ------------------------------------------------------------
# 18. Z-score each gene
#
# Row scaling makes the heatmap show relative expression
# patterns rather than absolute expression magnitude.
# ------------------------------------------------------------

leading_expr_scaled <- t(
    scale(
        t(leading_expr)
    )
)


# ------------------------------------------------------------
# 19. Sample annotation
# ------------------------------------------------------------

sample_info <- as.data.frame(
    colData(dds)
)

sample_info$patient_id <-
    as.character(
        sample_info$patient_id
    )

sample_info$timepoint <-
    as.character(
        sample_info$timepoint
    )

sample_info$sample_id <-
    rownames(sample_info)


# ------------------------------------------------------------
# 20. Order samples by patient and then timepoint
#
# This creates adjacent T1/T3 pairs.
# ------------------------------------------------------------

sample_info$timepoint_order <- ifelse(
    sample_info$timepoint == "T1",
    1,
    2
)

sample_info <- sample_info[
    order(
        sample_info$patient_id,
        sample_info$timepoint_order
    ),
]

sample_order <- sample_info$sample_id

leading_expr_scaled <-
    leading_expr_scaled[
        ,
        sample_order,
        drop = FALSE
    ]


# ------------------------------------------------------------
# 21. Create sample annotation for heatmap
# ------------------------------------------------------------

annotation_col <- sample_info[
    ,
    c(
        "patient_id",
        "timepoint"
    ),
    drop = FALSE
]

rownames(annotation_col) <-
    annotation_col$sample_id

annotation_col$sample_id <- NULL


# ------------------------------------------------------------
# 22. Create heatmap
# ------------------------------------------------------------
library(reshape2) # for melt()

# 1. GET CLUSTERED ROW ORDER (To mimic pheatmap's row clustering)
row_dist <- dist(leading_expr_scaled)
row_hc <- hclust(row_dist)
row_order <- rownames(leading_expr_scaled)[row_hc$order]

# 2. PREPARE THE EXPRESSION DATA FOR GGPLOT
expr_df <- as.data.frame(leading_expr_scaled)
expr_df$Gene <- rownames(expr_df)
expr_melted <- melt(expr_df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")

# Set the row factor levels so the genes are clustered properly
expr_melted$Gene <- factor(expr_melted$Gene, levels = row_order)

# 3. PREPARE THE ANNOTATION DATA FOR GGPLOT
# Make sure sample names are a column, not just row names
anno_df <- as.data.frame(annotation_col)
anno_df$Sample <- rownames(anno_df)
# Select your main annotation column name here (replace 'Your_Group_Column' with yours)
group_col <- colnames(annotation_col)[1]

# Merge expression data with annotations
plot_data <- merge(expr_melted, anno_df, by = "Sample", all.x = TRUE)

# 4. BUILD THE HEATMAP
p <- ggplot(plot_data, aes(x = Sample, y = Gene)) +
    # The expression heatmap
    geom_tile(aes(fill = Expression)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +

    # The top annotation bar (drawn using a secondary tile layer grouped at the top)
    # Note: To separate them cleanly, you can also use cowplot/patchwork,
    # but faceting by your annotation group is often much cleaner:
    facet_grid(. ~ get(group_col), scales = "free_x", space = "free_x") +

    # Styling to match your original pheatmap aesthetics
    theme_minimal() +
    theme(
        axis.text.x = element_blank(), # show_colnames = FALSE
        axis.text.y = element_text(size = 9), # fontsize_row = 9
        axis.ticks = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        strip.background = element_rect(fill = "gray90", color = NA), # Annotation labels
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
        title = "Leading-edge genes of dominant T3-associated GO programs",
        x = NULL,
        y = NULL
    )

# 5. GGSAVE IT WITHOUT FEAR
ggsave(
    filename = "Figure8_T3_leading_edge_heatmap.pdf",
    plot = p,
    width = 12,
    height = 9,
    units = "in"
)


# ------------------------------------------------------------
# 23. Save selected heatmap genes
# ------------------------------------------------------------

write.csv(
    candidate_leading_edge,
    "results/enrichment/T3_leading_edge_heatmap_genes.csv",
    row.names = FALSE
)


# ------------------------------------------------------------
# 24. Print summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("LEADING-EDGE ANALYSIS SUMMARY\n")
cat("============================================\n\n")

cat(
    "Selected GO terms:",
    nrow(selected_terms),
    "\n"
)

cat(
    "Leading-edge genes identified:",
    nrow(recurrence_table),
    "\n"
)

cat(
    "Genes displayed in heatmap:",
    nrow(candidate_leading_edge),
    "\n\n"
)

cat("Selected biological terms:\n")

print(
    selected_terms[
        ,
        c(
            "ID",
            "Description",
            "NES",
            "p.adjust"
        )
    ]
)

cat("\nTop recurrent leading-edge genes:\n")

print(
    candidate_leading_edge[
        ,
        c(
            "SYMBOL",
            "n_terms",
            "GENENAME"
        )
    ]
)

cat("\nFiles written:\n")

cat(
    "  results/enrichment/leading_edge_selected_T3_GO_terms.csv\n"
)

cat(
    "  results/enrichment/T3_leading_edge_gene_recurrence.csv\n"
)

cat(
    "  results/enrichment/T3_leading_edge_heatmap_genes.csv\n"
)

cat(
    "  results/figures/Figure8_T3_leading_edge_heatmap.pdf\n"
)

cat("\n============================================\n")
cat("SCRIPT 11 COMPLETED\n")
cat("============================================\n")
