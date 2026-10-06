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
