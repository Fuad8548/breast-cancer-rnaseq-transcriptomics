# ================================================================
# Project 1: Final figure generation
# ================================================================
#
# This script generates the final portfolio-ready figures from
# previously saved analysis objects/results.
#
# No new genome-wide statistical analysis is introduced here.
#
# Figures:
# 1. PCA
# 2. Volcano plot
# 3. GO GSEA summary
# 4. KEGG GSEA summary
# 5. Representative candidate paired-expression plot
# 6. Candidate expression heatmap
#
# Run from the project root.
# ================================================================


suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
})


dir.create(
    "results/figures",
    recursive = TRUE,
    showWarnings = FALSE
)


# ================================================================
# 1. Load DESeq2 object
# ================================================================

dds <- readRDS(
    "data/processed/dds_T1_T3_fitted.rds"
)


# ================================================================
# 2. VST transformation
# ================================================================

vsd <- vst(
    dds,
    blind = FALSE
)

vst_mat <- assay(vsd)

meta <- as.data.frame(
    colData(dds)
)


# ================================================================
# FIGURE 1
# PCA
# ================================================================

pca <- prcomp(
    t(vst_mat),
    center = TRUE,
    scale. = FALSE
)

percent_var <- 100 *
    (pca$sdev^2 /
        sum(pca$sdev^2))

pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
)

pca_df$timepoint <-
    meta[
        pca_df$sample_id,
        "timepoint"
    ]

pca_df$patient_id <-
    meta[
        pca_df$sample_id,
        "patient_id"
    ]


p_pca <- ggplot(
    pca_df,
    aes(
        x = PC1,
        y = PC2,
        label = patient_id,
        shape = timepoint
    )
) +
    geom_point(
        size = 3
    ) +
    labs(
        title = "PCA of breast cancer RNA-seq samples",
        x = paste0(
            "PC1: ",
            round(percent_var[1], 1),
            "% variance"
        ),
        y = paste0(
            "PC2: ",
            round(percent_var[2], 1),
            "% variance"
        ),
        shape = "Timepoint"
    ) +
    theme_bw()


ggsave(
    "results/figures/Figure1_PCA.pdf",
    p_pca,
    width = 8,
    height = 6
)


# ================================================================
# FIGURE 2
# Volcano plot
# ================================================================

res <- results(
    dds,
    contrast = c(
        "timepoint",
        "T3",
        "T1"
    )
)

volcano_df <- as.data.frame(
    res
)

volcano_df$gene_id <-
    rownames(volcano_df)

volcano_df$neglog10p <-
    -log10(
        pmax(
            volcano_df$pvalue,
            .Machine$double.xmin
        )
    )

volcano_df$direction <-
    ifelse(
        volcano_df$log2FoldChange > 0,
        "T3",
        "T1"
    )


candidate_symbols <- c(
    "KNL1",
    "CDK1",
    "CCNB1",
    "RAD51",
    "GYS1",
    "SCD",
    "INSR",
    "TNC"
)

gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(volcano_df$gene_id),
    columns = "SYMBOL",
    keytype = "ENSEMBL"
)

gene_map <- gene_map[
    !duplicated(gene_map$ENSEMBL),
]

volcano_df$ENSEMBL <-
    sub(
        "\\..*$",
        "",
        volcano_df$gene_id
    )

volcano_df <- merge(
    volcano_df,
    gene_map,
    by = "ENSEMBL",
    all.x = TRUE
)

volcano_df$label <-
    ifelse(
        volcano_df$SYMBOL %in%
            candidate_symbols,
        volcano_df$SYMBOL,
        NA
    )


p_volcano <- ggplot(
    volcano_df,
    aes(
        x = log2FoldChange,
        y = neglog10p
    )
) +
    geom_point(
        alpha = 0.55,
        size = 1.3
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    geom_text(
        data = subset(
            volcano_df,
            !is.na(label)
        ),
        aes(
            label = label
        ),
        vjust = -0.7,
        size = 3
    ) +
    labs(
        title = "T3 vs T1 differential expression",
        x = "log2 fold change",
        y = "-log10(p-value)"
    ) +
    theme_bw()


ggsave(
    "results/figures/Figure2_Volcano.pdf",
    p_volcano,
    width = 8,
    height = 6
)


# ================================================================
# FIGURE 3
# GO GSEA
# ================================================================

go_simple_df <- as.data.frame(
    readRDS(
        "results/enrichment/gsea_GO_BP_T3_vs_T1_simplified.rds"
    )
)

go_t3 <- go_simple_df[
    go_simple_df$NES > 0 &
        go_simple_df$p.adjust < 0.05,
]

go_t1 <- go_simple_df[
    go_simple_df$NES < 0 &
        go_simple_df$p.adjust < 0.05,
]

go_t3 <- go_t3[
    order(go_t3$p.adjust),
]

go_t1 <- go_t1[
    order(go_t1$p.adjust),
]


go_plot_df <- rbind(
    head(go_t3, 8),
    head(go_t1, 8)
)

go_plot_df$direction <-
    ifelse(
        go_plot_df$NES > 0,
        "T3",
        "T1"
    )


p_go <- ggplot(
    go_plot_df,
    aes(
        x = NES,
        y = reorder(
            Description,
            NES
        ),
        size = -log10(p.adjust)
    )
) +
    geom_point() +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    labs(
        title = "GO Biological Process GSEA",
        x = "Normalized enrichment score",
        y = NULL,
        size = "-log10(FDR)"
    ) +
    theme_bw()


ggsave(
    "results/figures/Figure3_GO_GSEA.pdf",
    p_go,
    width = 10,
    height = 8
)


# ================================================================
# FIGURE 4
# KEGG GSEA
# ================================================================

kegg_df <- read.csv(
    "results/enrichment/KEGG_GSEA_T3_vs_T1_all_pathways.csv",
    check.names = FALSE
)

kegg_t3 <- kegg_df[
    kegg_df$NES > 0 &
        kegg_df$p.adjust < 0.05,
]

kegg_t1 <- kegg_df[
    kegg_df$NES < 0 &
        kegg_df$p.adjust < 0.05,
]

kegg_t3 <- kegg_t3[
    order(kegg_t3$p.adjust),
]

kegg_t1 <- kegg_t1[
    order(kegg_t1$p.adjust),
]

kegg_plot_df <- rbind(
    head(kegg_t3, 8),
    head(kegg_t1, 8)
)


p_kegg <- ggplot(
    kegg_plot_df,
    aes(
        x = NES,
        y = reorder(
            Description,
            NES
        ),
        size = -log10(p.adjust)
    )
) +
    geom_point() +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    labs(
        title = "KEGG GSEA",
        x = "Normalized enrichment score",
        y = NULL,
        size = "-log10(FDR)"
    ) +
    theme_bw()


ggsave(
    "results/figures/Figure4_KEGG_GSEA.pdf",
    p_kegg,
    width = 10,
    height = 8
)


# ================================================================
# FIGURE 5
# Candidate paired expression
# ================================================================

candidate_ens <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = candidate_symbols,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
)

candidate_df_list <- lapply(
    candidate_symbols,
    function(sym) {
        ens <- unname(
            candidate_ens[sym]
        )

        data.frame(
            SYMBOL = sym,
            expression =
                as.numeric(
                    vst_mat[
                        ens,
                        colnames(vst_mat)
                    ]
                ),
            patient_id =
                as.character(
                    meta[
                        colnames(vst_mat),
                        "patient_id"
                    ]
                ),
            timepoint =
                as.character(
                    meta[
                        colnames(vst_mat),
                        "timepoint"
                    ]
                )
        )
    }
)

candidate_df <- do.call(
    rbind,
    candidate_df_list
)

candidate_df$timepoint <- factor(
    candidate_df$timepoint,
    levels = c(
        "T1",
        "T3"
    )
)


p_candidates <- ggplot(
    candidate_df,
    aes(
        x = timepoint,
        y = expression,
        group = patient_id
    )
) +
    geom_line(
        alpha = 0.35
    ) +
    geom_point(
        size = 1.8
    ) +
    facet_wrap(
        ~SYMBOL,
        scales = "free_y",
        ncol = 4
    ) +
    labs(
        title =
            "Patient-level validation of representative candidates",
        subtitle =
            "24 matched patients; VST-transformed expression",
        x = "Timepoint",
        y = "VST expression"
    ) +
    theme_bw()


ggsave(
    "results/figures/Figure5_Candidate_Paired_Expression.pdf",
    p_candidates,
    width = 12,
    height = 7
)


# ================================================================
# FIGURE 6
# Candidate heatmap
# ================================================================

candidate_matrix <- sapply(
    candidate_symbols,
    function(sym) {
        ens <- unname(
            candidate_ens[sym]
        )

        as.numeric(
            vst_mat[
                ens,
                colnames(vst_mat)
            ]
        )
    }
)

candidate_matrix <- t(
    candidate_matrix
)

rownames(candidate_matrix) <-
    candidate_symbols

colnames(candidate_matrix) <-
    colnames(vst_mat)


sample_order <- order(
    factor(
        meta[colnames(vst_mat), "timepoint"],
        levels = c("T1", "T3")
    ),
    meta[colnames(vst_mat), "patient_id"]
)

candidate_matrix_scaled <- t(
    scale(
        t(candidate_matrix)
    )
)

pdf(
    "results/figures/Figure6_Candidate_Heatmap.pdf",
    width = 12,
    height = 5
)

heatmap(
    candidate_matrix_scaled[
        ,
        sample_order
    ],
    Rowv = NA,
    Colv = NA,
    scale = "none",
    margins = c(
        10,
        8
    ),
    main =
        "Representative candidate genes"
)

dev.off()


cat(
    "\nFinal figure generation complete.\n"
)
