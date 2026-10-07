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
