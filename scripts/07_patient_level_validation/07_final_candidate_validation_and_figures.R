# Project 1: Final candidate validation and figures
# Run from the project root after 06_gsea_and_pathways.R.
#
# Biological question:
# Do representative pathway-convergent genes show their expected
# T1 -> T3 direction consistently across the 24 matched patients?
#
# These are exploratory candidate-level analyses. They do not replace
# the genome-wide DESeq2 inference.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

dir.create("results/candidates", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# 1. Load fitted DESeq2 object and VST expression
# -------------------------------------------------------------------

dds <- readRDS(
  "data/processed/dds_T1_T3_fitted.rds"
)

vsd <- vst(
  dds,
  blind = FALSE
)

vst_mat <- assay(vsd)

meta <- as.data.frame(colData(dds))

# -------------------------------------------------------------------
# 2. Representative genes selected from GO+KEGG convergence
# -------------------------------------------------------------------

candidate_info <- data.frame(
  SYMBOL = c(
    "KNL1", "CDK1", "CCNB1", "RAD51",
    "GYS1", "SCD", "INSR", "TNC"
  ),
  expected_direction = c(
    "T3", "T3", "T3", "T3",
    "T1", "T1", "T1", "T1"
  ),
  stringsAsFactors = FALSE
)

candidate_symbols <- candidate_info$SYMBOL

# Map SYMBOL -> Ensembl directly
candidate_ens <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = candidate_symbols,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)

candidate_check <- data.frame(
  SYMBOL = candidate_symbols,
  ENSEMBL = unname(candidate_ens),
  present_in_vst = unname(candidate_ens) %in% rownames(vst_mat)
)

write.csv(
  candidate_check,
  "results/candidates/candidate_gene_availability.csv",
  row.names = FALSE
)

if (!all(candidate_check$present_in_vst)) {
  missing <- candidate_check$SYMBOL[
    !candidate_check$present_in_vst
  ]
  stop(
    paste(
      "Candidate genes missing from VST matrix:",
      paste(missing, collapse = ", ")
    )
  )
}

# -------------------------------------------------------------------
# 3. Build long-format paired expression table
# -------------------------------------------------------------------

plot_df_list <- lapply(
  candidate_symbols,
  function(sym) {

    ens <- unname(candidate_ens[sym])

    data.frame(
      sample_id = colnames(vst_mat),
      SYMBOL = sym,
      expression = as.numeric(
        vst_mat[ens, colnames(vst_mat)]
      ),
      patient_id = as.character(
        meta[colnames(vst_mat), "patient_id"]
      ),
      timepoint = as.character(
        meta[colnames(vst_mat), "timepoint"]
      ),
      stringsAsFactors = FALSE
    )
  }
)

plot_df <- do.call(rbind, plot_df_list)

plot_df$timepoint <- factor(
  plot_df$timepoint,
  levels = c("T1", "T3")
)

# -------------------------------------------------------------------
# 4. Candidate-level paired statistics
# -------------------------------------------------------------------

validation_results <- lapply(
  candidate_symbols,
  function(sym) {

    sub <- plot_df[
      plot_df$SYMBOL == sym,
      c("patient_id", "timepoint", "expression")
    ]

    t1 <- sub[
      sub$timepoint == "T1",
      c("patient_id", "expression")
    ]

    t3 <- sub[
      sub$timepoint == "T3",
      c("patient_id", "expression")
    ]

    names(t1)[2] <- "T1"
    names(t3)[2] <- "T3"

    paired <- merge(
      t1,
      t3,
      by = "patient_id"
    )

    paired <- paired[
      complete.cases(paired),
    ]

    delta <- paired$T3 - paired$T1

    expected <- candidate_info[
      candidate_info$SYMBOL == sym,
      "expected_direction"
    ]

    concordant <- if (expected == "T3") {
      delta > 0
    } else {
      delta < 0
    }

    wt <- wilcox.test(
      paired$T3,
      paired$T1,
      paired = TRUE,
      exact = FALSE
    )

    data.frame(
      SYMBOL = sym,
      expected_direction = expected,
      n_pairs = length(delta),
      median_delta = median(delta),
      mean_delta = mean(delta),
      n_concordant = sum(concordant),
      concordance_percent = 100 * mean(concordant),
      paired_wilcox_p = wt$p.value
    )
  }
)

candidate_validation <- do.call(
  rbind,
  validation_results
)

# Candidate-level BH adjustment across the 8 preselected tests.
candidate_validation$paired_wilcox_p_adj <- p.adjust(
  candidate_validation$paired_wilcox_p,
  method = "BH"
)

write.csv(
  candidate_validation,
  "results/candidates/paired_candidate_validation.csv",
  row.names = FALSE
)

# -------------------------------------------------------------------
# 5. Paired expression plot
# -------------------------------------------------------------------

p <- ggplot(
  plot_df,
  aes(
    x = timepoint,
    y = expression,
    group = patient_id
  )
) +
  geom_line(alpha = 0.35) +
  geom_point(size = 1.8) +
  facet_wrap(
    ~ SYMBOL,
    scales = "free_y",
    ncol = 4
  ) +
  labs(
    title = "Paired expression of representative candidate genes",
    subtitle = "24 matched patients; VST-transformed expression",
    x = "Timepoint",
    y = "VST expression"
  ) +
  theme_bw()

ggsave(
  "results/figures/paired_candidate_expression.pdf",
  plot = p,
  width = 12,
  height = 7
)

# -------------------------------------------------------------------
# 6. Candidate heatmap data
# -------------------------------------------------------------------

candidate_matrix <- sapply(
  candidate_symbols,
  function(sym) {
    ens <- unname(candidate_ens[sym])
    as.numeric(vst_mat[ens, colnames(vst_mat)])
  }
)

candidate_matrix <- t(candidate_matrix)
rownames(candidate_matrix) <- candidate_symbols
colnames(candidate_matrix) <- colnames(vst_mat)

write.csv(
  candidate_matrix,
  "results/candidates/representative_candidate_VST_matrix.csv",
  row.names = TRUE
)

# -------------------------------------------------------------------
# 7. Compact final candidate summary
# -------------------------------------------------------------------

candidate_summary <- candidate_validation

candidate_summary$interpretation_role <- c(
  "T3 kinetochore/chromosome-segregation representative",
  "T3 cell-cycle representative",
  "T3 mitotic cell-cycle representative",
  "T3 DNA-repair/genome-maintenance representative",
  "T1 carbohydrate/glycogen metabolism representative",
  "T1 lipid metabolism representative",
  "T1 metabolic signaling representative",
  "T1 extracellular/tissue-organization representative"
)

write.csv(
  candidate_summary,
  "results/candidates/final_representative_candidates.csv",
  row.names = FALSE
)

cat("\nFinal candidate-validation stage complete.\n")
print(candidate_summary)
