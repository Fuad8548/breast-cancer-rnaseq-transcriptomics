# Stage 10: Patient-level validation of the candidate genes

Now we're asking a different question from DESeq2/GSEA:
> Does the T1→T3 shift appear consistently across individual matched patients, or is the pathway signal being driven by only a few patients?

For visualization, we'll use VST expression from our fitted DESeq2 object. VST stabilizes the relationship between expression magnitude and variance, making genes easier to visualize across samples.

We'll examine four representative genes from each side:
- T3: KNL1, CDK1, CCNB1, RAD51
- T1: GYS1, SCD, INSR, TNC
This gives us coverage of kinetochore/cell-cycle/DNA-repair biology on the T3 side and metabolic/signaling/ECM biology on the T1 side.

1. **Load the fitted DESeq2 object**
```r
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)

dds <- readRDS(
    "r_objects/dds_T1_T3_fitted.rds"
)
```

2. **Generate VST expression**
```r
vsd <- vst(
    dds,
    blind = FALSE
)

vst_mat <- assay(vsd)
```

The matrix is:
> genes × 48 samples

Now removing any Ensembl version suffix just in case:
```r
rownames(vst_mat) <- sub(
    "\\..*$",
    "",
    rownames(vst_mat)
)
```
3. **Mapping our candidate genes**
```r
gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = rownames(vst_mat),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
)

# Defining the candidates
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
```

4. **Building the patient-level expression table**
```r
meta <- as.data.frame(
    colData(dds)
)

sample_ids <- colnames(vst_mat)

plot_df <- data.frame()

for (sym in candidate_info$SYMBOL) {

    ens <- names(gene_symbols)[
        gene_symbols == sym
    ][1]

    if (is.na(ens)) {
        warning(
            paste("Gene not found:", sym)
        )
        next
    }

    tmp <- data.frame(
        sample_id = sample_ids,
        SYMBOL = sym,
        expression = as.numeric(
            vst_mat[ens, sample_ids]
        ),
        patient_id = meta[
            sample_ids,
            "patient_id"
        ],
        timepoint = meta[
            sample_ids,
            "timepoint"
        ],
        stringsAsFactors = FALSE
    )

    plot_df <- rbind(
        plot_df,
        tmp
    )
}

plot_df$timepoint <- factor(
    plot_df$timepoint,
    levels = c("T1", "T3")
)
```

**Output:**


5. **Making the paired plot**
```r
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
        subtitle = "24 matched patients, VST-transformed expression",
        x = "Timepoint",
        y = "VST expression"
    ) +
    theme_bw()

print(p)
```

Each line represents one patient:
> T1 ●────────● T3

- If the gene is generally T3-associated, we'd expect many lines to move upward.
- If T1-associated, we'd expect many to move downward.

6. **The plot is useful, but let's quantify what our eyes see**

We don't want to say:
> "Looks like most patients go upward."

Let's calculate:
```r
summary_list <- list()

for (sym in candidate_symbols) {

    sub <- plot_df[
        plot_df$SYMBOL == sym,
    ]

    wide <- reshape(
        sub[
            ,
            c(
                "patient_id",
                "timepoint",
                "expression"
            )
        ],
        idvar = "patient_id",
        timevar = "timepoint",
        direction = "wide"
    )

    delta <- wide$expression.T3 -
             wide$expression.T1

    expected <- candidate_info[
        candidate_info$SYMBOL == sym,
        "expected_direction"
    ]

    concordant <- if (
        expected == "T3"
    ) {
        delta > 0
    } else {
        delta < 0
    }

    summary_list[[sym]] <- data.frame(
        SYMBOL = sym,
        expected_direction = expected,
        n_pairs = sum(!is.na(delta)),
        median_delta = median(delta, na.rm = TRUE),
        mean_delta = mean(delta, na.rm = TRUE),
        n_concordant = sum(concordant, na.rm = TRUE),
        concordance_percent =
            mean(concordant, na.rm = TRUE) * 100
    )
}

candidate_validation <- do.call(
    rbind,
    summary_list
)

candidate_validation
```

7. **What does concordance_percent mean?**
Suppose KNL1 gives:
```text
n_pairs = 24
n_concordant = 20
concordance = 83.3%
```

That means:
> In 20 of the 24 matched patients, KNL1's expression moved in the expected T3 direction.
That's very different from simply reporting the average expression change.

And if we found:
> 12 / 24 = 50%

that would tell us the gene's direction isn't particularly consistent across patients.
This is why patient-level validation matters.

8. **A very important distinction**
We are not replacing DESeq2 with this analysis.
DESeq2 remains our formal differential-expression model:
> ~ patient_id + timepoint

This new analysis is asking:
> How consistent is the observed direction across individual paired patients?

So the two analyses answer complementary questions:
```text
DESeq2
→ Is there evidence for a systematic T1/T3 effect
  under the statistical model?

Paired plot
→ How does that effect behave in individual patients?
```

### Reproducibility
```r
dir.create(
    "results/candidates",
    recursive = TRUE,
    showWarnings = FALSE
)

# saving the table
write.csv(
    candidate_validation,
    "results/candidates/paired_candidate_validation.csv",
    row.names = FALSE
)

# saving the figure
ggsave(
    "results/candidates/paired_candidate_expression.pdf",
    plot = p,
    width = 12,
    height = 7
)
```

4. **Now let's interpret the validation correctly**
The output will give us:
```text
n_pairs
median_delta
mean_delta
n_concordant
concordance_percent
```

For example, conceptually:
> KNL1   T3   24   +...   +...   18   75%

would mean:
> 18 of the 24 patients showed an increase in KNL1 from T1 to T3.

That is within-patient consistency. This is a different piece of evidence from the population-level GSEA result.

Our evidence hierarchy is now:
```text
DESeq2
   ↓
overall modeled T3 vs T1 effect

GSEA
   ↓
coordinated pathway signal

GO + KEGG convergence
   ↓
same biological programs highlighted

Patient-level paired plot
   ↓
does the direction occur across individuals?
```

5. **Wilcoxon signed-rank test**
For these eight genes, we'll also calculate a paired Wilcoxon signed-rank test as a descriptive check of the within-patient shifts.
This is not replacing DESeq2, and it won't change our main conclusion. It simply tells us whether the paired expression differences tend systematically away from zero for each selected candidate.

```r
wilcox_results <- lapply(
    candidate_symbols,
    function(sym) {

        sub <- plot_df[
            plot_df$SYMBOL == sym,
        ]

        wide <- reshape(
            sub[
                ,
                c(
                    "patient_id",
                    "timepoint",
                    "expression"
                )
            ],
            idvar = "patient_id",
            timevar = "timepoint",
            direction = "wide"
        )

        test <- wilcox.test(
            wide$expression.T3,
            wide$expression.T1,
            paired = TRUE,
            exact = FALSE
        )

        data.frame(
            SYMBOL = sym,
            paired_wilcox_p = test$p.value
        )
    }
)

wilcox_results <- do.call(
    rbind,
    wilcox_results
)

wilcox_results
```

### What our results tell us
Three candidates show a nominally significant paired Wilcoxon result:
| Gene     | Direction | Concordance | Paired Wilcoxon p |
| -------- | --------- | ----------: | ----------------: |
| **KNL1** | T3        |       70.8% |        **0.0135** |
| **GYS1** | T1        |       70.8% |       **0.00894** |
| **INSR** | T1        |       70.8% |        **0.0135** |

The remaining five have paired Wilcoxon p-values **above 0.05**.

So, descriptively:

- KNL1 shows a fairly consistent T1→T3 increase across patients.
- GYS1 shows a fairly consistent T1-associated direction.
- INSR likewise shows a fairly consistent T1-associated direction.
- CDK1, CCNB1 and RAD51 have the expected T3 direction overall, but the individual-patient consistency is weaker.
- SCD and TNC show T1-associated direction, but their paired Wilcoxon evidence is weaker.

However, these are still exploratory candidate-gene results, not validated biomarkers. And importantly, **all eight still have DESeq2 FDR ≈ 0.869**, so our original conclusion of no individually FDR-significant DEGs remains unchanged.

**Inconsistency**
Our table says:
> KNL1   n_pait = 7    concordance = 70.83%

But 70.83% of 24 patients is:
> 17 / 24 = 70.83%

Similarly:
CDK1   54.17% = 13 / 24
CCNB1  62.50% = 15 / 24
RAD51  66.67% = 16 / 24

So the `n_pait` column doesn't agree with the concordance percentage.
That means we should not save this particular table as our final validation table.
The likely problem is that the earlier calculation of the `n_pairs` field was not actually counting the same observations used for the concordance calculation.

### Fix
```r
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
            concordance_percent =
                100 * mean(concordant),
            paired_wilcox_p = wt$p.value
        )
    }
)

candidate_validation <- do.call(
    rbind,
    validation_results
)

candidate_validation
```

**Output**:
```bash
SYMBOL expected_direction n_pairs median_delta mean_delta n_concordant
1   KNL1   T3     24    0.5692191  0.6028285   17
2   CDK1   T3     24    0.3536958  0.7429317   13
3  CCNB1   T3     24    0.4382296  0.5774797   15
4  RAD51   T3     24    0.2755488  0.4721754   16
5   GYS1   T1     24   -0.2817543 -0.4025785   17
6    SCD   T1     24   -0.3912577 -0.6163015   15
7   INSR   T1     24   -0.2723577 -0.3202184   17
8    TNC   T1     24   -0.3636894 -0.4244384   15
  concordance_percent paired_wilcox_p
1            70.83333     0.013457443
2            54.16667     0.100412494
3            62.50000     0.100412494
4            66.66667     0.089130926
5            70.83333     0.008941423
6            62.50000     0.149060718
7            70.83333     0.013457443
8            62.50000     0.106465467
```

We want every number to have a clear provenance:
```text
24 matched patients
        ↓
24 T1/T3 pairs
        ↓
T3 − T1 VST difference
        ↓
directional concordance
        ↓
paired Wilcoxon test
```

**And now we have an interesting biological distinction**
Our pathway analysis suggested:
**T3**
```text
mitosis
kinetochore
chromosome segregation
DNA replication
DNA repair
```
with candidate genes including **KNL1, CDK1, CCNB1, RAD51**. Yet patient-level behavior tells us that these genes aren't equally consistent.

For example, if after correction we confirm:
> KNL1 → 17/24 patients → 70.8%
> CDK1 → 13/24 → 54.2%

Then KNL1 appears to represent the T3-associated program more consistently across this cohort than CDK1, descriptively.
That doesn't make KNL1 a validated biomarker. It simply makes it a more useful representative candidate for our visualization and discussion.
Likewise, on the T1 side, GYS1 and INSR currently show the strongest paired evidence among our selected candidates.

**One methodological point**
The paired Wilcoxon p-values are not replacing DESeq2. Our primary genome-wide inference remains:
> No individual gene survives FDR correction.

The paired Wilcoxon analysis is a candidate-level descriptive validation performed after pathway-based selection.

So the hierarchy remains:
Genome-wide DESeq2
        ↓
no FDR-significant individual genes

GSEA
        ↓
strong coordinated pathway signals

GO + KEGG
        ↓
convergent biological programs

Candidate genes
        ↓
paired patient-level validation
        ↓
exploratory representatives

**T3** => **KNL1 is the clearest representative among our four T3 candidate**

| Gene     | State | Median T3−T1 VST change | Concordant pairs | Concordance | Paired Wilcoxon p |
| -------- | ----- | ----------------------: | ---------------: | ----------: | ----------------: |
| **KNL1** | T3    |                  +0.569 |            17/24 |       70.8% |        **0.0135** |

It has:
  - 17/24 patients moving in the expected direction
  - 70.8% concordance
  - median VST increase of +0.569
  - paired Wilcoxon p = 0.0135
And biologically it sits right in the center of the pathway signal: kinetochore/chromosome segregation.
RAD51 is also useful because it represents the **DNA-repair/genome-maintenance** side of the T3 program, although its patient-level test is weaker.
CDK1 is biologically important to the cell-cycle story, but only 13/24 patients show the expected direction, so it is less useful as our representative patient-level marker.

**T1:**
**GYS1 and INSR are the clearest representatives among the four T1 candidates.**

GYS1:
17/24 concordant, 70.8%, median change −0.282, p = 0.00894

INSR:
17/24 concordant, 70.8%, median change −0.272, p = 0.0135

They also fit different parts of the T1 program:

GYS1
 ↓
carbohydrate/glycogen metabolism

INSR
 ↓
insulin/metabolic signaling

That gives us complementary biological representation rather than picking four genes from the same pathway.

They also fit different parts of the T1 program:
```text
GYS1
 ↓
carbohydrate/glycogen metabolism

INSR
 ↓
insulin/metabolic signaling
```

### We've now completed the discovery/validation pipeline
Our Project analysis is essentially:
```text
Processed counts
      ↓
48 matched T1/T3 samples
      ↓
Filtering + normalization
      ↓
VST + PCA
      ↓
Paired DESeq2
      ↓
0 FDR-significant individual DEGs
      ↓
Genome-wide ranked gene list
      ↓
GO-BP GSEA
      ↓
GO redundancy reduction
      ↓
KEGG GSEA
      ↓
GO/KEGG core-gene convergence
      ↓
Candidate modules
      ↓
24-patient paired validation
      ↓
KNL1 / RAD51
GYS1 / INSR
```

And the resulting biological interpretation is now quite coherent:
- **T3-associated**: proliferative/cell-cycle and genome-maintenance machinery, particularly kinetochore, chromosome segregation, DNA replication and DNA repair.
- **T1-associated**: metabolic and cell/tissue-interaction-associated biology, particularly carbohydrate/metabolic signaling and extracellular/tissue organization.





















