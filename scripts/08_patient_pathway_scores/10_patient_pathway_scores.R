# ============================================================
# Project 1 — Paired T1 vs T3 Transcriptomics
# Script 10 — Patient-level Pathway Scores
#
# Purpose:
#   Quantify selected biological programs in each individual
#   sample using ssGSEA, then compare T3 vs T1 within patients.
#
# Important:
#   These pathway signatures are exploratory and were selected
#   from the biological themes observed in the GSEA analysis.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
    library(DESeq2)
    library(GSVA)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(ggplot2)
    library(tidyr)
    library(dplyr)
})


# ------------------------------------------------------------
# 2. Input / output paths
# ------------------------------------------------------------

dds_file <- "r_objects/dds_T1_T3_fitted.rds"

if (!file.exists(dds_file)) {
    stop(
        "Missing file: ", dds_file,
        "\nRun the upstream DESeq2 analysis first."
    )
}

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

dir.create(
    "r_objects",
    recursive = TRUE,
    showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Load fitted DESeq2 object
# ------------------------------------------------------------

dds <- readRDS(dds_file)


# ------------------------------------------------------------
# 4. Create variance-stabilized expression matrix
#
# ssGSEA is applied to continuous expression values.
# VST is appropriate here because the original DESeq2 object
# contains RNA-seq count information.
# ------------------------------------------------------------

vsd <- vst(
    dds,
    blind = FALSE
)

expr <- assay(vsd)

cat("\n============================================\n")
cat("Expression matrix\n")
cat("============================================\n")

cat("Genes:", nrow(expr), "\n")
cat("Samples:", ncol(expr), "\n")


# ------------------------------------------------------------
# 5. Remove Ensembl version suffixes
#
# This makes the identifiers compatible with annotation data.
# ------------------------------------------------------------

rownames(expr) <- sub(
    "\\..*$",
    "",
    rownames(expr)
)


# ------------------------------------------------------------
# 6. Define biological signatures
#
# These signatures summarize the major themes found in the
# previous GO analysis.
#
# T3-associated:
#   A. chromosome segregation / kinetochore biology
#   B. DNA replication / homologous-recombination repair
#
# T1-associated:
#   C. carbohydrate / small-molecule metabolism
#   D. extracellular matrix / cell-junction biology
#
# Several related GO terms are combined into broader signatures.
# ------------------------------------------------------------

signature_terms <- list(
    T3_Chromosome_Segregation = c(
        "GO:0007059", # chromosome segregation
        "GO:0008608", # spindle attachment to kinetochore
        "GO:0051304" # chromosome separation
    ),
    T3_DNA_Replication_Repair = c(
        "GO:0006261", # DNA-templated DNA replication
        "GO:0000724" # homologous recombination repair
    ),
    T1_Metabolism = c(
        "GO:0005975", # carbohydrate metabolic process
        "GO:0019318", # hexose metabolic process
        "GO:0044282" # small molecule catabolic process
    ),
    T1_ECM_Junction = c(
        "GO:0030198", # extracellular matrix organization
        "GO:0034329", # cell junction assembly
        "GO:0002009" # epithelial morphogenesis
    )
)


# ------------------------------------------------------------
# 7. Retrieve Ensembl genes belonging to each GO term
# ------------------------------------------------------------

cat("\nBuilding gene sets...\n")

build_go_gene_set <- function(go_ids, set_name) {
    annotation <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = go_ids,
        keytype = "GOALL",
        columns = c(
            "GOALL",
            "ENSEMBL",
            "SYMBOL"
        )
    )

    annotation <- annotation[
        !is.na(annotation$ENSEMBL),
    ]

    genes <- unique(
        annotation$ENSEMBL
    )

    genes <- intersect(
        genes,
        rownames(expr)
    )

    cat(
        set_name,
        ":",
        length(genes),
        "genes represented in expression matrix\n"
    )

    if (length(genes) < 5) {
        warning(
            "Very small gene set for ",
            set_name,
            ": ",
            length(genes),
            " genes."
        )
    }

    genes
}


gene_sets <- lapply(
    names(signature_terms),
    function(x) {
        build_go_gene_set(
            signature_terms[[x]],
            x
        )
    }
)

names(gene_sets) <- names(
    signature_terms
)

saveRDS(
    gene_sets,
    file = "r_objects/frozen_pathway_signatures.rds"
)


# ------------------------------------------------------------
# 8. Check gene-set overlap
# ------------------------------------------------------------

cat("\nGene-set sizes:\n")

gene_set_sizes <- data.frame(
    signature = names(gene_sets),
    n_genes = sapply(
        gene_sets,
        length
    )
)

print(gene_set_sizes)

write.csv(
    gene_set_sizes,
    "results/enrichment/patient_pathway_gene_set_sizes.csv",
    row.names = FALSE
)


# ------------------------------------------------------------
# 9. Compute ssGSEA scores
#
# Current GSVA versions use an ssgseaParam object followed by
# gsva().
#
# ssGSEA transforms the gene × sample matrix into:
#
#     pathway × sample
#
# Each value summarizes relative representation of a gene set
# within one individual sample.
# ------------------------------------------------------------

if (!exists(
    "ssgseaParam",
    where = asNamespace("GSVA"),
    inherits = FALSE
)) {
    stop(
        "Your installed GSVA package does not provide ssgseaParam(). ",
        "\nPlease update GSVA through Bioconductor."
    )
}

ssgsea_parameter <- GSVA::ssgseaParam(
    exprData = expr,
    geneSets = gene_sets,
    minSize = 5,
    maxSize = 1000,
    normalize = TRUE
)

ssgsea_result <- GSVA::gsva(
    ssgsea_parameter,
    verbose = TRUE
)

score_matrix <- as.matrix(
    ssgsea_result
)


# ------------------------------------------------------------
# 10. Save pathway score matrix
# ------------------------------------------------------------

write.csv(
    score_matrix,
    "results/enrichment/patient_pathway_ssGSEA_scores.csv"
)

saveRDS(
    score_matrix,
    "r_objects/patient_pathway_ssGSEA_scores.rds"
)


# ------------------------------------------------------------
# 11. Sample metadata
# ------------------------------------------------------------

sample_info <- as.data.frame(
    colData(dds)
)

sample_info$sample_id <- rownames(
    sample_info
)

sample_info$patient_id <- as.character(
    sample_info$patient_id
)

sample_info$timepoint <- as.character(
    sample_info$timepoint
)


# ------------------------------------------------------------
# 12. Convert pathway matrix into long format
# ------------------------------------------------------------

score_df <- as.data.frame(
    t(score_matrix),
    check.names = FALSE
)

score_df$sample_id <- rownames(
    score_df
)

score_long <- score_df %>%
    left_join(
        sample_info[
            ,
            c(
                "sample_id",
                "patient_id",
                "timepoint"
            )
        ],
        by = "sample_id"
    ) %>%
    pivot_longer(
        cols = all_of(names(gene_sets)),
        names_to = "signature",
        values_to = "ssGSEA_score"
    )

write.csv(
    score_long,
    "results/enrichment/patient_pathway_scores_long.csv",
    row.names = FALSE
)


# ------------------------------------------------------------
# 13. Calculate paired T3 - T1 differences
# ------------------------------------------------------------

paired_scores <- score_long %>%
    select(
        patient_id,
        timepoint,
        signature,
        ssGSEA_score
    ) %>%
    pivot_wider(
        names_from = timepoint,
        values_from = ssGSEA_score
    )

if (!all(c("T1", "T3") %in% colnames(paired_scores))) {
    stop(
        "Expected both T1 and T3 columns in pathway score data."
    )
}

paired_scores$delta_T3_minus_T1 <-
    paired_scores$T3 - paired_scores$T1


# ------------------------------------------------------------
# 14. Paired Wilcoxon testing
#
# This asks whether each pathway score systematically changes
# from T1 to T3 across the 24 matched patients.
# ------------------------------------------------------------

pathway_test <- lapply(
    split(
        paired_scores,
        paired_scores$signature
    ),
    function(df) {
        test <- wilcox.test(
            df$T3,
            df$T1,
            paired = TRUE,
            exact = FALSE
        )

        expected_direction <- case_when(
            grepl(
                "^T3_",
                unique(df$signature)
            ) ~ "T3 higher",
            grepl(
                "^T1_",
                unique(df$signature)
            ) ~ "T1 higher",
            TRUE ~ "unspecified"
        )

        data.frame(
            signature = unique(
                df$signature
            ),
            direction = expected_direction,
            n_pairs = sum(
                complete.cases(
                    df$T1,
                    df$T3
                )
            ),
            median_T1 = median(
                df$T1,
                na.rm = TRUE
            ),
            median_T3 = median(
                df$T3,
                na.rm = TRUE
            ),
            median_delta_T3_minus_T1 =
                median(
                    df$delta_T3_minus_T1,
                    na.rm = TRUE
                ),
            mean_delta_T3_minus_T1 =
                mean(
                    df$delta_T3_minus_T1,
                    na.rm = TRUE
                ),
            p_value = test$p.value
        )
    }
) %>%
    bind_rows()


# ------------------------------------------------------------
# 15. Multiple-testing correction
# ------------------------------------------------------------

pathway_test$p_adj <- p.adjust(
    pathway_test$p_value,
    method = "BH"
)

pathway_test <- pathway_test[
    order(
        pathway_test$p_adj,
        pathway_test$p_value
    ),
]

write.csv(
    pathway_test,
    "results/enrichment/patient_pathway_paired_tests.csv",
    row.names = FALSE
)


# ------------------------------------------------------------
# 16. Concordance across patients
#
# For each signature, calculate how many patients move in the
# expected direction.
# ------------------------------------------------------------

expected_direction_function <- function(signature) {
    if (grepl("^T3_", signature)) {
        return(1)
    }

    if (grepl("^T1_", signature)) {
        return(-1)
    }

    return(NA)
}


concordance_table <- paired_scores %>%
    group_by(signature) %>%
    summarise(
        n_pairs = sum(
            complete.cases(
                T1,
                T3
            )
        ),
        median_delta = median(
            delta_T3_minus_T1,
            na.rm = TRUE
        ),
        expected_direction =
            expected_direction_function(
                unique(signature)
            ),
        patients_expected_direction =
            sum(
                expected_direction *
                    delta_T3_minus_T1 > 0,
                na.rm = TRUE
            ),
        percent_expected_direction =
            100 *
                patients_expected_direction /
                n_pairs,
        .groups = "drop"
    )

write.csv(
    concordance_table,
    "results/enrichment/patient_pathway_concordance.csv",
    row.names = FALSE
)


# ------------------------------------------------------------
# 17. Final pathway-score figure
#
# Each line connects the T1 and T3 score for the same patient.
#
# This is the key visualization:
#
#       T1 ─────────→ T3
#
# repeated across patients.
# ------------------------------------------------------------

plot_df <- score_long

p <- ggplot(
    plot_df,
    aes(
        x = timepoint,
        y = ssGSEA_score,
        group = patient_id
    )
) +
    geom_line(
        alpha = 0.35
    ) +
    geom_point(
        size = 2.3
    ) +
    facet_wrap(
        ~signature,
        scales = "free_y"
    ) +
    labs(
        title = "Patient-level pathway activity",
        subtitle =
            "Paired ssGSEA scores for selected biological programs",
        x = NULL,
        y = "ssGSEA score"
    ) +
    theme_classic() +
    theme(
        strip.text = element_text(
            face = "bold"
        ),
        plot.title = element_text(
            face = "bold"
        )
    )


ggsave(
    "results/figures/Figure7_patient_level_pathway_scores.pdf",
    p,
    width = 11,
    height = 8,
    units = "in"
)


# ------------------------------------------------------------
# 18. Print results
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PATIENT-LEVEL PATHWAY RESULTS\n")
cat("============================================\n\n")

print(
    pathway_test
)

cat("\nConcordance:\n\n")

print(
    concordance_table
)

cat("\nFiles written:\n")

cat(
    "  results/enrichment/patient_pathway_gene_set_sizes.csv\n"
)

cat(
    "  results/enrichment/patient_pathway_ssGSEA_scores.csv\n"
)

cat(
    "  results/enrichment/patient_pathway_scores_long.csv\n"
)

cat(
    "  results/enrichment/patient_pathway_paired_tests.csv\n"
)

cat(
    "  results/enrichment/patient_pathway_concordance.csv\n"
)

cat(
    "  results/figures/Figure7_patient_level_pathway_scores.pdf\n"
)

cat(
    "  r_objects/patient_pathway_ssGSEA_scores.rds\n"
)

cat("\n============================================\n")
cat("SCRIPT 10 COMPLETED\n")
cat("============================================\n")
