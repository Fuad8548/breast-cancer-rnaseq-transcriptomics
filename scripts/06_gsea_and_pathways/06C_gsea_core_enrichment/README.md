## What our core_enrichment result is telling us
For example:
```bash
GO:0008608
attachment of spindle microtubules to kinetochore
NES = +2.76
FDR = 7.68 × 10⁻⁹
```

The `core_enrichment` field contains Entrez IDs such as:

```text
57082/151648/83540/221150/983/...
```
These are the genes contributing to the enrichment signal for that GO term. Because we used:

```r
keyType = "ENTREZID"
```
those numbers are **Entrez Gene identifiers**.

**One subtle point we should remember**
Look at:
> negative regulation of chromosome segregation
with `NES = +2.69`

Someone might incorrectly say:
> "T3 negatively regulates chromosome segregation."

That's not what the result means.
The GO term is simply the biological category assigned to those genes.

The positive NES means:
> Genes annotated to the GO process "negative regulation of chromosome segregation" are enriched toward the T3 side of our ranked gene list.

It does not directly tell us that the biological process itself is inhibited.

## Now: which genes are repeatedly driving the T3 signal?
This is what we want to calculate.
Notice that our top 10 pathways are highly related:
```text
spindle/kinetochore
      ↕
nuclear division
      ↕
chromosome segregation
      ↕
cell cycle
      ↕
DNA replication
```
Therefore many genes will appear in several `core_enrichment` lists.

A gene appearing in 6 related pathways is interesting, but remember:
> those pathways are themselves related, so recurrence does not mean six independent discoveries.

Still, it is an excellent way to identify candidate central genes/modules for our biological interpretation.

## Extract the leading-edge genes
```r
library(AnnotationDbi)
library(org.Hs.eg.db)

t3_top10 <- go_simple_df[
    go_simple_df$NES > 0 & go_simple_df$p.adjust < 0.05,
]

t3_top10 <- t3_top10[
    order(t3_top10$p.adjust),
]

t3_top10 <- t3_top10[
    1:min(10, nrow(t3_top10)),
]

t3_top10[
    ,
    c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]
```

We already have essentially that output, but now let's **split the long strings into individual genes**.

```r
# 1. Start fresh from the top 10 data
t3_long <- t3_top10

# 2. Split the genes into a clean, flat vector first
gene_list <- strsplit(t3_long$core_enrichment, "/")
expanded_genes <- unlist(gene_list)

# 3. Repeat the rows of the dataframe based on the original split counts
t3_long <- t3_long[rep(seq_len(nrow(t3_long)), lengths(gene_list)), ]

# 4. Now assigning matches perfectly because both are exactly 59,846 items long
t3_long$ENTREZID <- expanded_genes

# 5. Safely drop the old column
t3_long$core_enrichment <- NULL
```

Now `t3_long` contains something conceptually like:

```bash
GO term       NES     Entrez
--------------------------------
kinetochore   2.76    57082
kinetochore   2.76    151648
kinetochore   2.76    83540
...
```

## Now translate those mysterious numbers into gene names

```r
gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(t3_long$ENTREZID),
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
)

t3_long$ENTREZID <- as.character(t3_long$ENTREZID)

t3_long <- merge(
    t3_long,
    gene_map,
    by = "ENTREZID",
    all.x = TRUE
)

head(
    t3_long[
        ,
        c("ENTREZID", "SYMBOL", "GENENAME",
          "Description", "NES")
    ],
    30
)
```

Let's ask:
> Which genes occur repeatedly across our top T3-enriched pathways?

```r
gene_frequency <- t3_long[
    !duplicated(t3_long[, c("Description", "ENTREZID")]),
]

gene_frequency <- aggregate(
    Description ~ ENTREZID + SYMBOL + GENENAME,
    data = gene_frequency,
    FUN = length
)

names(gene_frequency)[4] <- "n_top_GO_terms"

gene_frequency <- gene_frequency[
    order(-gene_frequency$n_top_GO_terms),
]

head(gene_frequency, 30)
```

## Why this is better than simply taking the biggest log2FC

We are following:
```text
DESeq2
   ↓
Wald statistic
   ↓
ranked transcriptome
   ↓
GSEA
   ↓
significant biological programs
   ↓
core/leading-edge genes
   ↓
recurrent genes
   ↓
candidate biological players
```

We have:
```text
GO:0008608   spindle → kinetochore
GO:0000280   nuclear division
GO:0051985   chromosome segregation
GO:0050000   chromosome localization
GO:0006261   DNA replication
GO:0007094   spindle checkpoint
GO:0071174   spindle checkpoint
GO:0090068   cell cycle
GO:0007059   chromosome segregation
GO:0007088   mitotic nuclear division
```

that gives us a biological model:
```text
                 T3-associated state
                         │
          ┌──────────────┼──────────────┐
          ↓              ↓              ↓
      Cell cycle    Chromosome       Genome
                     segregation     replication
          │              │              │
          ↓              ↓              ↓
      mitosis       kinetochore       DNA replication
                     spindle          DNA repair
                     checkpoint
```

**Now the pathway signal has turned into actual molecular biology**

1. **Look at the genes themselves**

A. **Kinetochore machinery**
We have:
```text
KNL1
NDC80
NUF2
SPC24
SPC25
SKA1
SKA3
ZW10
ZWINT
ZWILCH
SPDL1
```
These are heavily concentrated around the **kinetochore and spindle attachment machinery.** Think of a dividing chromosome like this:

```text
chromosome
    │
centromere
    │
kinetochore
    │
spindle microtubules
    │
opposite cell pole
```

The kinetochore is the protein assembly that connects chromosomes to spindle microtubules during mitosis.
So having **KNL1 + NDC80 + NUF2 + SPC24 + SPC25 + SKA1 + SKA3 + ZW10**-family genes repeatedly appearing across your enriched GO terms is much more convincing than seeing one isolated cell-cycle gene.

2. **Then there's the spindle checkpoint machinery**
We also have:
```text
BUB1
BUB1B
MAD2L1
MAD2L1BP
TTK
TRIP13
CDC20
```
These are associated with the machinery that monitors chromosome-spindle attachment and controls progression through mitosis.

So our result is becoming:
```text
             T3-associated signal
                     │
          ┌──────────┴──────────┐
          │                     │
    Kinetochore            Spindle checkpoint
          │                     │
   KNL1, NDC80             BUB1, BUB1B
   NUF2, SPC24             MAD2L1, TTK
   SPC25, SKA1             TRIP13, CDC20
   SKA3, ZW10
```

3. **And we see direct cell-cycle regulators**
Our list also contains:
```text
CCNB1
CDK1
CDC20
BIRC5
CENPF
CDCA8
```
These are consistent with active mitotic/cell-cycle transcriptional programs
So we're seeing **multiple layers of the same biological system**:
```text
Cell-cycle control
       ↓
Mitotic entry/progression
       ↓
Spindle formation & checkpoint
       ↓
Kinetochore attachment
       ↓
Chromosome segregation
```
That's much more informative than simply saying:
> "GO analysis found cell-cycle pathways."

We can now say **which molecular machinery is repeatedly driving the enrichment**.

4. **The frequency result is particularly useful**
Look at our top three:
```text
KNL1   → 9 top GO terms
NDC80  → 9
SKA1   → 9
SKA3   → 9
```
Then:
```text
BIRC5  → 8
CCNB1  → 8
NUF2   → 8
ZW10   → 8
```
and:
```text
CDCA8
GEN1
SPC24
SPC25
TTK
ZNF207
ZWINT
ZWILCH
```
This tells us that the **same molecular genes repeatedly contribute to several related enriched processes**.

However, there's an important statistical caveat:
> `n_top_GO_terms` is not an independent evidence score.

For example, KNL1 may occur in:
```text
chromosome segregation
nuclear division
kinetochore attachment
mitotic nuclear division
```
because those biological concepts are themselves related.
So we should use this frequency as a way to identify **candidate central genes**, not as proof that KNL1 is "9 times more important."

5. **We can now formulate a much stronger biological hypothesis**
At this point, our evidence chain is:
```text
DESeq2
│
├── no individual gene: FDR < 0.05
│
└── ranked Wald statistics
         ↓
       GSEA
         ↓
T3-associated GO programs
         ↓
mitosis
chromosome segregation
spindle checkpoint
DNA replication
DNA repair
         ↓
core-enrichment genes
         ↓
KNL1, NDC80, NUF2, SKA1,
SKA3, SPC24, SPC25,
BUB1, BUB1B, MAD2L1,
TTK, CCNB1, CDK1, CDC20...
```

6. **There's an especially nice finding here**

Our original DESeq2 result gave us:
> 0 FDR-significant genes

Yet GSEA found an extremely strong pathway-level signal.

7. **Before calling these "candidate biomarkers", we need one more check**

A recurrent GSEA gene is not automatically a biomarker.
For a portfolio-quality analysis, we should connect each candidate back to the original DESeq2 statistics:
```text
Gene
 ↓
log2 fold change
 ↓
Wald statistic
 ↓
p-value
 ↓
FDR
 ↓
GSEA core enrichment
 ↓
biological role
```
That tells us whether, for example, KNL1 is:

- strongly T3-associated,
- moderately T3-associated,
- highly variable,
- or merely included in the pathway because of the gene-set structure.

8. Let's take these recurrent T3 genes and merge them with the original DESeq2 results.

```r
library(DESeq2)

candidate_t3 <- gene_frequency

de_results <- as.data.frame(
    results(
        readRDS("data/processed/dds_T1_T3_fitted.rds"),
        contrast = c("timepoint", "T3", "T1")
    )
)

de_results$ENTREZID <- mapIds(
    org.Hs.eg.db,
    keys = rownames(de_results),
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
)

candidate_t3$ENTREZID <- as.character(candidate_t3$ENTREZID)
de_results$ENTREZID <- as.character(de_results$ENTREZID)

candidate_t3 <- merge(
    candidate_t3,
    de_results[, c(
        "ENTREZID",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj"
    )],
    by = "ENTREZID",
    all.x = TRUE
)

candidate_t3 <- candidate_t3[
    order(
        -candidate_t3$n_top_GO_terms,
        -abs(candidate_t3$log2FoldChange)
    ),
]

candidate_t3[
    ,
    c(
        "SYMBOL",
        "GENENAME",
        "n_top_GO_terms",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj"
    )
]
```

| Gene      | GO-term recurrence | log2FC T3 vs T1 | nominal p | DESeq2 FDR |
| --------- | -----------------: | --------------: | --------: | ---------: |
| **KNL1**  |                  8 |          +1.048 |   0.00365 |      0.869 |
| **NUF2**  |                  7 |          +1.065 |   0.00401 |      0.869 |
| **CENPF** |                  6 |          +0.940 |   0.01195 |      0.869 |
| **SKA3**  |                  8 |          +0.862 |   0.00899 |      0.869 |
| **SPC25** |                  8 |          +0.800 |   0.03236 |      0.869 |
| **CDK1**  |                  6 |          +0.755 |   0.02038 |      0.869 |
| **NDC80** |                  8 |          +0.735 |   0.03516 |      0.869 |
| **BUB1**  |                  6 |          +0.732 |   0.04017 |      0.869 |
| **CCNB1** |                  8 |          +0.682 |   0.02527 |      0.869 |
| **SPC24** |                  7 |          +0.679 |   0.02288 |      0.869 |
| **SPDL1** |                  6 |          +0.538 |   0.00232 |      0.869 |
| **BUB1B** |                  5 |          +0.704 |   0.04267 |      0.869 |





























































