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
    "r_objects/dds_T1_T3_fitted.rds"
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
