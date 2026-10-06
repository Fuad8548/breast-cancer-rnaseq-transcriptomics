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
