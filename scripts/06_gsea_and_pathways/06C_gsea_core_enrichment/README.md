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
gene_list <- strsplit(
    t3_long$core_enrichment, 
    "/"
)
expanded_genes <- unlist(gene_list)

# 3. Repeat the rows of the dataframe based on the original split counts
t3_long <- t3_long[
    rep(
        seq_len(nrow(t3_long)), 
        lengths(gene_list)
    ), 
]

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
        readRDS("r_objects/dds_T1_T3_fitted.rds"),
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

The key thing is that the GSEA result and the DESeq2 result are telling two different layers of the same story.

1. **Look at the strongest recurrent genes**
From your merged table, these stand out:

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

The important pattern is that virtually all of these genes have **positive log2FC**.
So they are not merely appearing in the same GO pathways by coincidence. They are also positioned on the **T3-associated side of our DESeq2 ranking**.
For example:
```bash
KNL1    +1.05
NUF2    +1.07
CENPF   +0.94
SKA3    +0.86
SPC25   +0.80
CDK1    +0.76
NDC80   +0.73
BUB1    +0.73
CCNB1   +0.68
```

2. **Now notice the apparent contradiction**
Take KNL1:
```text
log2FC = +1.048
p = 0.00365
padj = 0.869
```
p = 0.0036 doesn't mean KNL1 is significant. We performed thousands of gene-level tests. After multiple-testing correction:
> FDR = 0.869

So KNL1 is **not statistically significant as an individual DEG** under our DESeq2 criterion.
The same applies to NUF2, CENPF, CDK1, etc.
This is precisely why our earlier statement:
> 0 FDR-significant individual genes
remains completely correct.

3. **Yet these genes collectively make a very strong GSEA signal**
We have:
```text
                 T3-associated ranking
                         ↑
                         │
        KNL1 ────────────┤
        NUF2 ────────────┤
        CENPF ───────────┤
        SKA3 ────────────┤
        CDK1 ────────────┤
        NDC80 ───────────┤
        BUB1 ────────────┤
        CCNB1 ───────────┤
        SPC24 ───────────┤
        SPC25 ───────────┤
                         │
                         ↓
                 T1-associated
```
A whole collection of genes involved in the same machinery is shifted toward T3.

Therefore GSEA can detect the **coordinated movement of the gene set**, even though none of the individual genes survives genome-wide correction.
That is the central statistical lesson of this project.

4. **We can now describe the T3 program much more precisely**
Previously we said:
> "T3 is associated with cell-cycle pathways."

Now we can be more specific => our data support a T3-associated program involving:
**Kinetochore assembly and chromosome attachment**

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
```

**Mitotic checkpoint:**
```text
BUB1
BUB1B
MAD2L1
TTK
TRIP13
SPDL1
```

**Mitotic progression**
```text
CDK1
CCNB1
CDC20
PLK1
CENPF
BIRC5
CDCA8
```

**DNA replication / genome maintenance**
We also see genes such as:
```text
CDT1
CDC6
MCM-family genes
RAD51
BRCA2
FANCD2
ATR
...
```
So the signal is not one isolated process.
It spans several connected aspects of proliferative and genome-maintenance biology.

So the appropriate language for this project is:
> exploratory candidate genes associated with the T3-enriched mitotic program
rather than:
> confirmed biomarkers.

Before we rank candidates, we should ensure:
> one Entrez ID = one gene = one recurrence count per GO term.

The GSEA itself isn't invalidated by this, because the underlying ranked list had already been deduplicated at the Entrez level. But our post-GSEA candidate summary should be cleaned.


5. **Our current biological story is now**

```text
T1
│
├── epithelial/tissue organization
├── extracellular structure
├── metabolic processes
└── vascular/developmental processes
│
│
╞═══════════════════════════════╡
│         T1 → T3
│
└── T3
    ├── mitotic cell cycle
    ├── chromosome segregation
    ├── kinetochore/spindle machinery
    ├── spindle checkpoint
    ├── DNA replication
    ├── DNA repair
    └── immune-associated programs
```

And within the T3 program:
```text
KNL1 ─ NDC80 ─ NUF2 ─ SPC24/SPC25
             │
          kinetochore
             │
    SKA1 ─ SKA3 ─ ZW10
             │
       spindle attachment
             │
   BUB1 ─ BUB1B ─ MAD2L1 ─ TTK
             │
       spindle checkpoint
             │
      CDK1 ─ CCNB1 ─ CDC20
             │
        mitotic progression
```

6. T1 core enrichment
We have properly characterized the T3 side.
We should now perform exactly the same process for the negative-NES pathways:
```text
T1-associated GO
       ↓
core enrichment genes
       ↓
remove duplicate mappings
       ↓
merge with DESeq2
       ↓
identify recurring T1 genes
```
We want to know not only:
> "What becomes more T3-associated?"

but also:
> "What molecular programs characterize the T1 side?"

## T1 associated genes
We now ask the mirror-image question:
> Which genes are driving the pathways enriched toward the T1 side of the transcriptome?

We already saw T1-associated GO themes such as **cell junction assembly, synapse organization, epithelial morphogenesis, lipid/carbohydrate metabolism, vasculature development, and small-molecule catabolism**.

1. **Select the strongest T1 pathways**
```r
t1_top <- go_simple_df[
    go_simple_df$NES < 0 & go_simple_df$p.adjust < 0.05,
]

t1_top <- t1_top[
    order(t1_top$p.adjust),
]

t1_top <- t1_top[
    1:min(10, nrow(t1_top)),
]

t1_top[
    ,
    c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]
```

What this does => We keep only:
> NES < 0
because negative NES means the gene set is enriched toward the **T1-associated end** of our DESeq2 ranking.
Then we keep the 10 strongest pathways.

2. **Break `core_enrichment` into individual genes**
The long string:
> 79834/8934/1400/79414/...
is actually many Entrez IDs separated by `/`.
We'll convert it into rows:

```r
t1_long <- t1_top

gene_list <- strsplit(
    t1_long$core_enrichment, 
    "/"
)
expanded_genes <- unlist(gene_list)

t1_long <- t1_long[
    rep(
        seq_len(nrow(t1_long)), 
        lengths(gene_list)
    ), 
]

t1_long$ENTREZID <- expanded_genes

t1_long$core_enrichment <- NULL
```

Conceptually, we have converted:
```text
Pathway A → 10/20/30/40
Pathway B → 20/30/50/60
```
into:
```text
Pathway A   10
Pathway A   20
Pathway A   30
Pathway A   40
Pathway B   20
Pathway B   30
Pathway B   50
Pathway B   60
```
Now we can count which genes repeatedly occur.

3. **Remove duplicate gene-pathway pairs**
This is the cleanup we identified from the T3 analysis.
```r
t1_long <- unique(
    t1_long[
        ,
        c(
            "Description",
            "NES",
            "p.adjust",
            "ENTREZID"
        )
    ]
)
```
This means a particular gene can contribute once per pathway, rather than being accidentally counted twice because of **annotation duplication**.

4. **Annotate the Entrez IDs**
```r
gene_map_t1 <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(t1_long$ENTREZID),
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
)
```

Then:
```r
gene_map_t1 <- gene_map_t1[
    !duplicated(gene_map_t1$ENTREZID),
]

t1_long <- merge(
    t1_long,
    gene_map_t1,
    by = "ENTREZID",
    all.x = TRUE
)
```
Now our mysterious Entrez numbers become interpretable gene symbols.

5. **Find the recurring T1 genes**
```r
t1_frequency <- aggregate(
    Description ~ ENTREZID + SYMBOL + GENENAME,
    data = t1_long,
    FUN = length
)

names(t1_frequency)[4] <- "n_top_GO_terms"

t1_frequency <- t1_frequency[
    order(-t1_frequency$n_top_GO_terms),
]

head(t1_frequency, 30)
```
This is exactly analogous to what we did for T3.

Output:
```bash
ENTREZID   SYMBOL          GENENAME     
207         AKT1    AKT serine/threonine kinase 1
25          ABL1    ABL proto-oncogene 1, non-receptor tyrosine kinase
351         APP     Amyloid beta precursor protein
n_top_GO_terms
6
5
5
.........
```
Again, this is not a statistical importance score. It simply tells us which genes recur across several of the strongest T1-associated GO terms.

6. **Now connect those genes back to DESeq2**
We don't want to say:
> "This gene is important because it appears in many GO terms."

We want to know:
> "Does this gene actually have T1-associated expression in the original DESeq2 analysis?"

**Step 1: Load the original DESeq2 results**
```r
de_results <- read.csv(
    "results/differential_expression/DESeq2_T3_vs_T1_all_genes.csv",
    check.names = FALSE
)
```
Then inspect the columns:
```r
names(de_results)
```
We want `gene_id`, `log2FoldChange`, `pvalue`, `padj`.

**Step 2: Clean the Ensembl IDs**
Our original DESeq2 result contains Ensembl IDs.
```r
de_results$ENSEMBL <- sub(
    "\\..*$",
    "",
    de_results$gene_id
)
```
This removes an Ensembl version suffix if one exists.
For example:
> ENSG00000123456.7
becomes:
> ENSG00000123456
This makes annotation more robust.

**Step 3: Map Ensembl → Entrez**
```r
de_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(de_results$ENSEMBL),
    columns = "ENTREZID",
    keytype = "ENSEMBL"
)
```
Now clean duplicate mappings:
```r
de_map <- de_map[
    !is.na(de_map$ENTREZID),
]

de_map <- de_map[
    !duplicated(de_map$ENSEMBL),
]
```

**Step 4: Attach Entrez IDs to DESeq2 results**
```r
de_results <- merge(
    de_results,
    de_map,
    by = "ENSEMBL",
    all.x = TRUE
)
```
Now our DESeq2 table has:
```bash
ENSEMBL
ENTREZID
baseMean
log2FoldChange
stat
pvalue
padj
```

**Step 5: Make one DESeq2 record per Entrez gene**
This is important because we want the same identifier system as our GSEA ranking.

```r
de_entrez <- de_results[
    !is.na(de_results$ENTREZID),
]

de_entrez <- de_entrez[
    order(
        de_entrez$ENTREZID,
        -abs(de_entrez$stat)
    ),
]

de_entrez <- de_entrez[
    !duplicated(de_entrez$ENTREZID),
]
```
**Why are we doing this?**
Suppose annotation produces:
```bash
Ensembl A → Entrez 123
Ensembl B → Entrez 123
```
We don't want Entrez 123 appearing twice in our candidate table.
And because our GSEA ranking is based on the DESeq2 Wald statistic, we retain the record with the largest absolute `stat`.

**Step 6: Merge with our T1 GSEA genes**
We already created:
```r
t1_frequency
```
which contains:
```text
ENTREZID
SYMBOL
GENENAME
n_top_GO_terms
```
Now:
```r
t1_candidates <- merge(
    t1_frequency,
    de_entrez[
        ,
        c(
            "ENTREZID",
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "padj",
            "stat"
        )
    ],
    by = "ENTREZID",
    all.x = TRUE
)
```
Then:
```r
t1_candidates <- t1_candidates[
    order(
        -t1_candidates$n_top_GO_terms,
        t1_candidates$log2FoldChange
    ),
]
```
Finally:
```r
t1_candidates[
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

Our earlier T3 candidate table had duplicate symbols such as SKA1 and SPC25 because the post-GSEA annotation contained duplicate mapping records. We're fixing that here with:

```r
!duplicated(de_map$ENSEMBL)
```

and:
```r
!duplicated(de_entrez$ENTREZID)
```

**Output**
Because these genes came from **negative-NES T1-associated pathways**, we're interested in seeing negative values such as:

```bash
SYMBOL n_top_GO_terms log2FoldChange pvalue padj  DVL1      6            -0.503   0.00024 0.869 
EPHB2     6            -0.367   0.08080 0.869 
DAG1      6            -0.235   0.19421 0.869 
AKT1      6            -0.232   0.12393 0.869 
ABL1      6            -0.177   0.29830 0.886
....
```
1. **The T1 genes are pointing toward a different biological module**
Our recurrent T1-associated genes include:
**DVL1, EPHB2, DAG1, AKT1, ABL1, SRPX2, PLXND1, SNCA, THBS2, NECTIN1, NRP1, APP, SEMA4C, VLDLR, SYNDIG1, RTN4RL2, NLGN3, CLSTN2, LEP, AMOT**.

And importantly, their log2 fold changes are predominantly negative, for example:
```bash
DVL1      -0.503
EPHB2     -0.367
DAG1      -0.235
AKT1      -0.232
VLDLR     -0.783
SYNDIG1   -0.736
RTN4RL2   -0.703
NLGN3     -0.597
AMOT      -0.496
```
So these aren't merely genes that happened to occur in negatively enriched pathways. Their individual expression estimates also point toward the T1 side.

2. **There are three major T1 themes**
- Theme 1: Cell-cell interaction and adhesion
```text
EPHB2
DAG1
NECTIN1
AMOT
PLXND1
```

These fit naturally with the earlier T1 GO findings involving:
> cell junction assembly and epithelial morphogenesis.
That gives us a coherent tissue-organization signal rather than a single isolated gene.

- Theme 2: Neural/synaptic-associated annotation
```text
SNCA
APP
VLDLR
SYNDIG1
RTN4RL2
NLGN3
CLSTN2
SEMA4C
```
And remember that one of our strongest T1 GO terms was:
> synapse organization, NES ≈ −1.95

So there really is a reproducible synaptic/neural-associated annotation signal on the T1 side.

It does not establish that the breast tumors are becoming neuronal. A safer interpretation is:
> The T1 transcriptional state is enriched for genes involved in cell-cell interaction, neuronal/synaptic-associated processes, and membrane/cell-adhesion biology.

- Theme 3: Extracellular/vascular biology
We also have:
```text
THBS2
NRP1
PLXND1
```
which fits nicely with the earlier T1 enrichment for:
```text
vasculature development
blood-vessel development
extracellular structure organization
extracellular matrix-related biology.
```
So the T1 side isn't just "synapse." It looks more like a broader **tissue architecture / cell interaction / extracellular and vascular-associated program**, with several genes also carrying synaptic-associated GO annotations.

3. **Compare this with T3**
**T3**
We found:
```text
KNL1
NDC80
NUF2
SKA1
SKA3
SPC24
SPC25
BUB1
BUB1B
TTK
MAD2L1
CDK1
CCNB1
CDC20
CENPF
...
```
These genes converge on:
**mitosis → kinetochore → spindle checkpoint → chromosome segregation → DNA replication/repair**

**T1**
We found:
```text
DVL1
EPHB2
DAG1
PLXND1
THBS2
NECTIN1
NRP1
APP
NLGN3
VLDLR
...
```

These converge on:
**cell interaction → adhesion/junctions → tissue organization → extracellular/vascular biology → synaptic-associated processes**

So we have something much closer to a state transition model:
```text
T1                                      T3
│                                        │
│ tissue architecture                   │ proliferative machinery
│ cell junctions                        │ mitotic machinery
│ cell-cell interaction                 │ kinetochore
│ extracellular/vascular biology        │ spindle checkpoint
│ metabolic processes                   │ chromosome segregation
│ synaptic-associated annotation        │ DNA replication/repair
│                                        │
└────────────── transcriptional shift ───┘
```

4. **FDR Column**
Take DVL1:
```bash
log2FC = -0.503
p = 0.000243
padj = 0.869
```

And VLDLR:
```bash
log2FC = -0.783
p = 0.00490
padj = 0.869
```

So, just like KNL1 on the T3 side, these are not statistically significant individual DEGs after multiple-testing correction.

Therefore:
**DVL1 is not a validated differentially expressed biomarker in this analysis.**

But it can be considered an exploratory candidate gene contributing to the **T1-associated pathway signal**.

5. **Our whole Project result is now much clearer**
We started with:
> "Are there individual genes that significantly differ between matched T1 and T3 samples?"

Answer:
> None after FDR correction.

Then we asked:
> "Are there coordinated biological programs?"

Answer:
> Yes

And now:
> "Which genes contribute to those coordinated programs?"

We have actual molecular candidates on both sides.

**T3 candidate module**
> KNL1 / NDC80 / NUF2 / SKA1 / SKA3 / SPC24 / SPC25 / BUB1 / TTK / CDK1 / CCNB1 / CENPF...

**T1 candidate module**
> DVL1 / EPHB2 / DAG1 / PLXND1 / THBS2 / NECTIN1 / NRP1 / VLDLR / NLGN3 / CLSTN2...

6. One methodological point before KEGG pathway
Our `n_top_GO_terms` values are useful for identifying recurring genes, but they **should not be treated as a gene importance score**.

For example:
```text
DVL1 → 6 GO terms
EPHB2 → 6 GO terms
```
doesn't mean DVL1 is automatically more biologically important than a gene appearing in only 3 terms.

The correct interpretation is:
> These genes recur across multiple related enriched GO gene sets and therefore are useful exploratory representatives of the underlying transcriptional program.

**What we have done so far:**
```text
DESeq2
   ↓
ranked gene list
   ↓
GO-BP GSEA
   ↓
GO redundancy reduction
   ↓
T3 core genes
   ↓
T1 core genes
```
The next question is:
> Do we see the same biological signal using an independent pathway database?

That's what KEGG gives us.

- GO organizes genes into biological processes.
- KEGG organizes them into named molecular/cellular pathways.

[KEGG GSEA](/scripts/06_gsea_and_pathways/06D_kegg_gsea/README.md)










