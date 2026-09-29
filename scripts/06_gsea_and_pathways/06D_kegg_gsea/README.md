# Module 06D: KEGG GSEA

Now we move from Gene Ontology to KEGG.
This is useful because GO and KEGG organize biological knowledge differently. We want to ask whether the strong signal we found with GO is also visible when the genes are grouped into **named molecular pathways**.

Our question is:
> Which KEGG pathways are shifted toward the T3 side or T1 side of our ranked transcriptome?

1. **First, understanding what KEGG is adding**
So far:
```text
DESeq2
   ↓
rank all genes
   ↓
GO-BP GSEA
   ↓
biological processes
```

GO told us things like:
```text
chromosome segregation
DNA replication
cell junction assembly
vasculature development
```

KEGG asks a slightly different question:
```text
Which recognized molecular pathways
contain genes concentrated toward T3 or T1?
```

So we might encounter pathways such as:
```text
Cell cycle
DNA replication
Homologous recombination
...
```
The important point is that KEGG is complementary evidence, not a replacement for GO.

2. **Reconstruct the ranked list**
Start from our permanent checkpoint:
```r
library(clusterProfiler)
library(org.Hs.eg.db)

gene_list <- readRDS(
    "data/processed/gene_list_T3_vs_T1_Entrez.rds"
)

length(gene_list)
```
**Output:**
> 12665
These are Entrez IDs ranked by the **DESeq2 Wald statistic**.

3. **Run KEGG GSEA**
```r
gsea_kegg <- gseKEGG(
    geneList = gene_list,
    organism = "hsa",
    keyType = "ncbi-geneid",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE
)
```
**Explanation**:
- `organism = "hsa"` => hsa means Homo sapiens: we're therefore asking KEGG to use human pathways.
- `keyType = "ncbi-geneid"`: Our ranked list contains Entrez/NCBI Gene IDs, so we tell `gseKEGG()` what identifier system we're using.
- `minGSSize` / `maxGSSize` => Same reasoning as GSEA before: 
> 10 ≤ pathway size ≤ 500
We don't want extremely tiny or excessively broad pathways dominating the analysis.
- `pAdjustMethod = "BH"`: Again, we're testing many pathways, so we need multiple-testing correction.

4. **Saving the result**

```r
dir.create(
    "results/enrichment",
    recursive = TRUE,
    showWarnings = FALSE
)

# save gsea_kegg result
saveRDS(
    gsea_kegg,
    "results/enrichment/gsea_KEGG_T3_vs_T1.rds"
)

# Load the saved result
gsea_kegg_res <- readRDS("results/enrichment/gsea_KEGG_T3_vs_T1.rds")

kegg_df <- as.data.frame(gsea_kegg_res)

write.csv(
    kegg_df,
    "results/enrichment/KEGG_GSEA_T3_vs_T1_all_pathways.csv",
    row.names = FALSE
)
```

5. **How do we read the result?**
The important columns are again:
```text
ID
Description
NES
pvalue
p.adjust
qvalue
core_enrichment
```
The central one is **NES**.

Remember:
```text
NES > 0  → pathway enriched toward T3
NES < 0  → pathway enriched toward T1
```
And:
```bash
p.adjust < 0.05
```
means the pathway survives multiple-testing correction.

6. **Strongest overall KEGG pathways**
```r
head(
    kegg_df[
        order(kegg_df$p.adjust),
        c(
            "ID",
            "Description",
            "NES",
            "pvalue",
            "p.adjust"
        )
    ],
    20
)
```

7. **Then separating T3 and T1**
**T3**

```r
kegg_t3 <- kegg_df[
    kegg_df$NES > 0 & kegg_df$p.adjust < 0.05,
]

kegg_t3 <- kegg_t3[
    order(kegg_t3$p.adjust),
]

kegg_t3[
    1:min(15, nrow(kegg_t3)),
    c(
        "ID",
        "Description",
        "NES",
        "p.adjust"
    )
]
```

**T1**
```r
kegg_t1 <- kegg_df[
    kegg_df$NES < 0 & kegg_df$p.adjust < 0.05,
]

kegg_t1 <- kegg_t1[
    order(kegg_t1$p.adjust),
]

kegg_t1[
    1:min(15, nrow(kegg_t1)),
    c(
        "ID",
        "Description",
        "NES",
        "p.adjust"
    )
]
```

8. **What are we hoping to learn?**
We already have a hypothesis from GO.

**GO suggested T3:**
```text
mitosis
chromosome segregation
kinetochore
spindle checkpoint
DNA replication
DNA repair
```

So if KEGG independently produces pathways related to:
```text
cell cycle
DNA replication
homologous recombination
```
then the two analyses are telling a consistent story.

Likewise, GO suggested T1:
```text
cell junction
epithelial organization
metabolism
vascular/extracellular biology
```
If KEGG identifies corresponding metabolic or signaling pathways, that gives us another layer of support.
Suppose KEGG finds something different. That's still valuable. We would ask why:
> Are the databases defining pathways differently? Is one signal more specific to GO annotations? Are there pathways represented poorly in KEGG?

That's proper biological analysis.

9. **Why KEGG does not invalidate our earlier GO result**

```text
GO map                         KEGG map

biological process             named pathway
       │                              │
chromosome segregation         cell cycle
DNA replication                DNA replication
cell junction                  signaling pathway
vascular development           metabolic pathway
       │                              │
       └──────── same transcriptome ──┘
```
They're two different ways of organizing the same genes. If the same broad biological phenomenon appears in both, our interpretation becomes stronger.

After running the following KEGG command, we get:

```r
head(
    kegg_df[
        order(kegg_df$p.adjust),
        c("ID","Description","NES","pvalue","p.adjust")
    ],
    20
)
```

**Output:** 
```bash
ID             Description 
hsa03010       Ribosome 
hsa04110       Cell cycle 
hsa05322       Systemic lupus erythematosus
hsa03030       DNA replication 
hsa00190       Oxidative phosphorylation 
hsa05171       Coronavirus disease 
hsa03040       Spliceosome

         NES      pvalue       p.adjust 
hsa03010 2.748215 1.000000e-10 1.750000e-08 
hsa04110 2.468838 1.000000e-10 1.750000e-08 
hsa05322 2.444788 1.914121e-08 2.233141e-06 
hsa03030 2.484317 1.290855e-07 1.129498e-05 
hsa00190 2.031015 5.140934e-06 2.998878e-04 
hsa05171 1.940099 4.835630e-06 2.998878e-04 
hsa03040 1.919139 8.515894e-06 3.977932e-04
```

This KEGG result is giving us a second, fairly coherent view of the same biology we saw with GO. And there are a couple of things here that we absolutely need to understand before we treat the pathway names as biological conclusions.

1. **The first major result: KEGG strongly confirms the T3 cell-cycle signal**
Our top KEGG pathways include:
| KEGG pathway                 |   NES |       FDR |
| ---------------------------- | ----: | --------: |
| **Ribosome**                 | +2.75 | 1.75×10⁻⁸ |
| **Cell cycle**               | +2.47 | 1.75×10⁻⁸ |
| **DNA replication**          | +2.48 | 1.13×10⁻⁵ |
| **Homologous recombination** | +2.17 | 5.94×10⁻⁴ |
| **Spliceosome**              | +1.92 | 3.98×10⁻⁴ |

So remember what GO independently told us:
```text
GO:
T3 → mitosis
     chromosome segregation
     kinetochore
     spindle checkpoint
     DNA replication
     DNA repair
```

Now KEGG says:
```text
KEGG:
T3 → Cell cycle
      DNA replication
      Homologous recombination
      Ribosome
      Spliceosome
```

That's important because two different pathway annotation systems are converging on a T3-associated proliferative/genome-maintenance program.

1. **Why "Ribosome" ranks 1?**
```text
Ribosome
NES = +2.75
FDR = 1.75 × 10⁻⁸
```
It means ribosome-associated genes are disproportionately concentrated toward the T3 end of our ranked transcriptome.
Biologically, that is compatible with increased emphasis on protein synthesis / translational machinery in the T3 transcriptional state.

And that fits reasonably well with the broader picture:
```text
T3
│
├── cell cycle
├── DNA replication
├── chromosome segregation
├── DNA repair
└── protein-production machinery
```

3. **The Cell Cycle + DNA Replication result is especially valuable**
Our GO analysis gave:
> DNA-templated DNA replication, mitotic division, chromosome segregation, spindle-checkpoint processes.

KEGG gives:
> Cell cycle + DNA replication + homologous recombination.

4. **Now the immune-related pathways**
Our KEGG results also contain:
```bash
Systemic lupus erythematosus        NES +2.44
Graft-versus-host disease           NES +2.13
Hematopoietic cell lineage          NES +2.01
Antigen processing/presentation     NES +2.06
Allograft rejection                 NES +2.04
T-cell receptor signaling           NES +1.81
```

This fits surprisingly well with our GO findings:
```text
adaptive immune response
lymphocyte differentiation
T-cell activation
immune effector process
leukocyte activation
```
So there appears to be a second T3-associated immune-related program.

But **"Systemic lupus erythematosus" does NOT mean our patients have lupus.**

Likewise:
> "Coronavirus disease" does not mean coronavirus infection.

Why? because KEGG disease pathways are collections of genes associated with those diseases. Many genes are shared with generic processes such as:

- immune signaling
- cytokine signaling
- antigen presentation
- inflammation
- cell proliferation

So if our transcriptome contains a strong immune signal, several apparently unrelated disease pathways can become enriched because they share parts of that molecular machinery.

> KEGG enrichment included several immune-associated disease pathways, consistent with the immune-related GO enrichment observed in the T3 transcriptional state. These pathway labels should be interpreted as shared molecular signatures rather than evidence of the corresponding diseases.

**"Why is coronavirus appearing in my breast-cancer dataset?"**
But KEGG is effectively saying:
> Some genes that participate in the molecular network represented by this pathway are collectively enriched toward T3.
It is not diagnosing infection.

6. **We also have some T1-associated KEGG pathways**
Two of your top 20 have negative NES:

```text
Other glycan degradation
NES = -2.21
Cytoskeleton in muscle cells
NES = -1.73
```
So these are enriched toward the **T1 side**.

That is interesting because GO had already shown a T1-associated signal involving:
```text
cell junctions
tissue organization
extracellular biology
metabolism
```

7. **Let's compare GO and KEGG directly**

| Biological theme               | GO       | KEGG                                                                  |
| ------------------------------ | -------- | --------------------------------------------------------------------- |
| Cell division                  | ✅        | **Cell cycle ✅**                                                      |
| DNA replication                | **✅**    | **✅**                                                                 |
| Genome repair                  | **✅**    | **Homologous recombination ✅**                                        |
| Mitosis/chromosome segregation | **✅**    | **Cell cycle ✅**                                                      |
| Immune response                | **✅**    | **T-cell signaling / antigen presentation / hematopoietic lineage ✅** |
| Metabolic biology              | **✅ T1** | Some metabolic-related pathways including glycan degradation          |
| Tissue/cell interaction        | **✅ T1** | Not yet obvious among top 20                                          |

9. **KEGG disease names are not the endpoints we want to highlight**
So, now run:
**T3-associated KEGG**

```r
kegg_t3 <- kegg_df[
    kegg_df$NES > 0 & kegg_df$p.adjust < 0.05,
]

kegg_t3 <- kegg_t3[
    order(kegg_t3$p.adjust),
]

kegg_t3[
    1:min(15, nrow(kegg_t3)),
    c("ID", "Description", "NES", "p.adjust")
]
```

**T1-associated KEGG**
```r
kegg_t1 <- kegg_df[
    kegg_df$NES < 0 & kegg_df$p.adjust < 0.05,
]

kegg_t1 <- kegg_t1[
    order(kegg_t1$p.adjust),
]

kegg_t1[
    1:min(15, nrow(kegg_t1)),
    c("ID", "Description", "NES", "p.adjust")
]
```

1. **T3: KEGG strongly reinforces the proliferative/genome-maintenance program**
Our strongest T3-associated pathways are:
```bash
Ribosome                  NES +2.75
Cell cycle                NES +2.47
DNA replication           NES +2.48
Homologous recombination  NES +2.17
Spliceosome               NES +1.92
```

These match our GO findings extremely well. We already had:
```text
GO → T3
     chromosome segregation
     mitosis
     kinetochore
     spindle checkpoint
     DNA replication
     DNA repair
```

Now:
```text
KEGG → T3
       cell cycle
       DNA replication
       homologous recombination
       ribosome
       spliceosome
```

So we have **convergence between two independent pathway annotation systems**.

2. **T3 looks broader than just "cell proliferation"**
```text
Ribosome
Cell cycle
DNA replication
Homologous recombination
Spliceosome
Oxidative phosphorylation
Nucleocytoplasmic transport
Ribosome biogenesis
```

So the T3-associated state appears to involve several coordinated programs:
- Genome duplication and division
> Cell cycle + DNA replication + homologous recombination

- Chromosome organization
> From GO: kinetochore + spindle + chromosome segregation

- Gene-expression/protein-production machinery
> Spliceosome + ribosome + ribosome biogenesis

That's why we would describe the T3 result as something like:
> a coordinated proliferative and genome-maintenance transcriptional program, accompanied by enrichment of RNA-processing and protein-production machinery.

3. **T3 also has a strong immune-associated signal**
KEGG gives:
```text
Systemic lupus erythematosus
Graft-versus-host disease
Hematopoietic cell lineage
Antigen processing and presentation
Viral protein interaction with cytokine/cytokine receptor
Allograft rejection
T-cell receptor signaling
```

And GO independently gave:
```text
adaptive immune response
lymphocyte differentiation
T-cell activation
immune effector process
leukocyte activation
```

So there is a second coherent T3 theme:
> immune-associated transcriptional activity.

These are pathway/gene modules containing genes that can also participate in general immune, inflammatory, signaling, or cellular processes.

4. **Now look at the T1 side**
The T1-associated KEGG pathways are actually quite informative:
```bash
Other glycan degradation       NES -2.21
Butanoate metabolism           NES -1.99
beta-Alanine metabolism        NES -1.96
AMPK signaling                 NES -1.72
Hormone signaling              NES -1.69
Integrin signaling             NES -1.62
Gap junction                   NES -1.67
Focal adhesion                 NES -1.51
Cadherin signaling             NES -1.46
Insulin signaling              NES -1.63
Proteoglycans in cancer        NES -1.59
```

Now compare that with our GO result:
```text
cell junction assembly
epithelial morphogenesis
extracellular matrix organization
vasculature development
carbohydrate metabolism
lipid metabolic regulation
small-molecule catabolism
```
The two analyses are converging again.

5. **So the T1 program can now be characterized more clearly**
Rather than calling it simply "metabolism", we'd divide it into two related components.

**T1 tissue-interaction/signaling program**
```text
Integrin signaling
Gap junction
Focal adhesion
Cadherin signaling
Proteoglycans
```

This fits beautifully with:
```text
GO:
cell junction assembly
epithelial morphogenesis
extracellular structure organization
```
This suggests a **cell-cell / cell-matrix / tissue-architecture-associated transcriptional state**.

**T1 metabolic/signaling program**
```text
glycan degradation
butanoate metabolism
beta-alanine metabolism
AMPK signaling
insulin signaling
hormone signaling
```
This fits the earlier GO metabolic enrichment.

6. **We can now put the whole project together**
```text
                  T1 ↔ T3 transcriptional shift
                            │
             ┌──────────────┴──────────────┐
             │                             │
            T1                            T3
             │                             │
   Tissue architecture        Proliferative program
   Cell-cell interaction      Cell cycle
   Cell-matrix interaction    Mitosis
   Integrin/cadherin          Chromosome segregation
   Gap junction               Kinetochore/spindle
   ECM-related biology        DNA replication
                              DNA repair
             │                RNA/protein machinery
   Metabolic/signaling        Immune-associated prog
```

7. **The most important scientific finding**
We started from:
```text
DESeq2:
0 genes with FDR < 0.05
```

It would conclude:
> "There is no difference between T1 and T3."

But our analysis found:
```text
Individual genes
       ↓
no genome-wide significant DEGs

BUT

Ranked transcriptome
       ↓
GSEA
       ↓
strong coordinated pathway changes
       ↓
GO + KEGG convergence
```
So the more accurate conclusion is:
> The analysis did not identify individually FDR-significant DEGs, but it did identify coordinated pathway-level differences between the matched T1 and T3 transcriptional states.

8. **Now we need to connect KEGG back to genes**
We found T3 core genes such as:
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
CDK1
CCNB1
CENPF
```

And T1-associated candidates such as:
```text
DVL1
EPHB2
DAG1
PLXND1
THBS2
NECTIN1
NRP1
VLDLR
NLGN3
CLSTN2
```
Now KEGG lets us ask:
> Are the same genes driving both the GO and KEGG signals?

9. **Let's extract KEGG core genes**
First, the T3 pathways:
```r
kegg_t3_top <- kegg_df[
    kegg_df$NES > 0 & kegg_df$p.adjust < 0.05,
]

kegg_t3_top <- kegg_t3_top[
    order(kegg_t3_top$p.adjust),
]

kegg_t3_top <- kegg_t3_top[
    1:min(10, nrow(kegg_t3_top)),
]

kegg_t3_top[
    ,
    c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]
```

and T1:
```r
kegg_t1_top <- kegg_df[
    kegg_df$NES < 0 & kegg_df$p.adjust < 0.05,
]

kegg_t1_top <- kegg_t1_top[
    order(kegg_t1_top$p.adjust),
]

kegg_t1_top <- kegg_t1_top[
    1:min(10, nrow(kegg_t1_top)),
]

kegg_t1_top[
    ,
    c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]
```

## GO + KEGG convergence
Until now we treated the analyses separately:
```text
GO GSEA
   ↓
biological processes

KEGG GSEA
   ↓
named pathways
```

Now we ask:
> Do GO and KEGG point to the same genes and the same biological programs?

1. **First, let's make the top GO and KEGG gene sets**
We'll use the strongest 10 pathways on each side.
**T3 GO:**
```r
go_t3_top <- go_simple_df[
    go_simple_df$NES > 0 & go_simple_df$p.adjust < 0.05,
]

go_t3_top <- go_t3_top[
    order(go_t3_top$p.adjust),
]

go_t3_top <- go_t3_top[
    1:min(10, nrow(go_t3_top)),
]
```

**T3 KEGG:**
```r
kegg_t3_top <- kegg_df[
    kegg_df$NES > 0 & kegg_df$p.adjust < 0.05,
]

kegg_t3_top <- kegg_t3_top[
    order(kegg_t3_top$p.adjust),
]

kegg_t3_top <- kegg_t3_top[
    1:min(10, nrow(kegg_t3_top)),
]
```

2. **Extract the actual genes**
We had `core_enrichment` that contained Entrez IDs separated by `/`

We'll turn them into unique gene sets.
```r
extract_core_genes <- function(df) {
    ids <- unlist(
        strsplit(
            df$core_enrichment,
            "/"
        )
    )
    
    unique(as.character(ids))
}

go_t3_genes <- extract_core_genes(go_t3_top)
kegg_t3_genes <- extract_core_genes(kegg_t3_top)
```

Now:
```r
length(go_t3_genes)
length(kegg_t3_genes)
```

These numbers tell us how many unique core-enrichment genes are represented by the top 10 GO and KEGG pathways.

3. **Find the genes shared by GO and KEGG**

```r
shared_t3 <- intersect(
    go_t3_genes,
    kegg_t3_genes
)

length(shared_t3)
```

We're asking:
> How many genes appear in the core enrichment of both the strongest GO and KEGG T3 pathways?

4. **Converting those numbers into gene names**

```r
shared_t3_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = shared_t3,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
)

shared_t3_map <- shared_t3_map[
    !duplicated(shared_t3_map$ENTREZID),
]

head(
    shared_t3_map[
        ,
        c("ENTREZID", "SYMBOL", "GENENAME")
    ],
    30
)
```

5. **Let's make the result even more useful**
We're not interested only in how many genes overlap.
We want to know:
> Which genes are supported by both GO and KEGG AND are actually shifted toward T3 in DESeq2?

Run:
```r
t3_convergence <- merge(
    shared_t3_map,
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


t3_convergence <- t3_convergence[
    order(
        -t3_convergence$log2FoldChange
    ),
]


head(
    t3_convergence[
        ,
        c(
            "SYMBOL",
            "GENENAME",
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "padj"
        )
    ],
    30
)
```

6. **This is powerful, because**
Imagine we find:
```text
KNL1
NDC80
NUF2
CDK1
CCNB1
BUB1
...
```

Then the evidence becomes:
```text
                     T3
                      │
             ┌────────┴────────┐
             │                 │
            GO               KEGG
             │                 │
       mitosis /           cell cycle /
       kinetochore         DNA replication
             │                 │
             └────────┬────────┘
                      ↓
               shared genes
                      ↓
              KNL1, NDC80,
              NUF2, CDK1...
```

That's much stronger than saying:
> "GO found cell cycle and KEGG found cell cycle."

We can show that the molecular genes contributing to both analyses overlap.

7. **We need to do the exact same thing for T1**
**T1 GO:**
```r
go_t1_top <- go_simple_df[
    go_simple_df$NES < 0 & go_simple_df$p.adjust < 0.05,
]

go_t1_top <- go_t1_top[
    order(go_t1_top$p.adjust),
]

go_t1_top <- go_t1_top[
    1:min(10, nrow(go_t1_top)),
]
```

**T1 KEGG**
```r
kegg_t1_top <- kegg_df[
    kegg_df$NES < 0 & kegg_df$p.adjust < 0.05,
]

kegg_t1_top <- kegg_t1_top[
    order(kegg_t1_top$p.adjust),
]

kegg_t1_top <- kegg_t1_top[
    1:min(10, nrow(kegg_t1_top)),
]

go_t1_genes <- extract_core_genes(go_t1_top)
kegg_t1_genes <- extract_core_genes(kegg_t1_top)

shared_t1 <- intersect(
    go_t1_genes,
    kegg_t1_genes
)

length(shared_t1)
```

And annotate:
```r
shared_t1_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = shared_t1,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
)

shared_t1_map <- shared_t1_map[
    !duplicated(shared_t1_map$ENTREZID),
]

head(
    shared_t1_map[
        ,
        c("ENTREZID", "SYMBOL", "GENENAME")
    ],
    30
)
```

8. **What we're ultimately building**
Our final evidence chain will look like:
```text
                      T1 ↔ T3
                         │
             ┌───────────┴───────────┐
             │                       │
            T1                      T3
             │                       │
       GO + KEGG                 GO + KEGG
             │                       │
             ↓                       ↓
   tissue/signaling/             cell cycle/
   metabolic programs            genome-maintenance
             │                       │
             ↓                       ↓
       shared genes              shared genes
             │                       │
             └───────────┬───────────┘
                         ↓
                 candidate genes
                         ↓
               biological model
```

Now, running the T3 and T1 convergence part through
```r
length(shared_t3)
head(shared_t3_map[, c("ENTREZID","SYMBOL","GENENAME")], 30)

length(shared_t1)
head(shared_t1_map[, c("ENTREZID","SYMBOL","GENENAME")], 30)
```

## GO and KEGG are converging on the same molecular genes?

We got:
```text
T3 GO ↔ T3 KEGG shared core genes = 88
T1 GO ↔ T1 KEGG shared core genes = 123
```
These are genes that occur in the core-enrichment sets of both the top GO and top KEGG pathways.
This does not mean 88 and 123 genes were independently validated. Both analyses use the same RNA-seq data and the same DESeq2-ranked gene list.

What it means is:
> Two different biological annotation frameworks are highlighting overlapping molecular components of the same transcriptional states.

1. **T3 convergence is extremely coherent**
Our first 30 T3 shared genes include:
```text
KNL1
SGO1
CDK1
CCNB1
NDC80
RB1
CDT1
CDC25C
RAD51
BUB1
BUB1B
BRCA2
RAD50
WAPL
BRIP1
TTK
PTTG1
CHEK2
TRIP13
...
```

Notice what happens when we group them:
**Kinetochore/chromosome segregation**
```text
KNL1
NDC80
SGO1
BUB1
BUB1B
TTK
WAPL
PTTG1
```

**Cell-cycle control**
```text
CDK1
CCNB1
CDC25C
RB1
```

**DNA replication/ repair**
```text
CDT1
RAD51
BRCA2
RAD50
BRIP1
RAD54L
RAD54B
CHEK2
```

Those give us:
```text
T3
 │
 ├── Cell-cycle progression
 │
 ├── Kinetochore / spindle checkpoint
 │
 ├── Chromosome segregation
 │
 └── DNA replication / repair
```

And that agrees with both our GO and KEGG analyses.

2. **T1 convergence is also remarkably coherent**
Our T1 shared genes begin:
```text
RAF1
ITGB3
LAMA5
ACTB
LAMB1
ITGA3
THBS2
ACTN1
CAMK2B
SORT1
DAG1
AGRN
PDLIM5
FZD5
TLN2
AKT1
MYH10
INSR
AGT
CAV1
ITGB4
FN1
JUP
TLN1
VCL
PLEC
TGFB1
MMP2
FLNA
...
```
Now the biological pattern becomes obvious.
**Integrin / cell-matrix interaction**
```text
ITGA3
ITGB3
ITGB4
LAMA5
LAMB1
FN1
TLN1
TLN2
VCL
DAG1
```

**Cytoskeleton / adhesion**
```text
ACTB
ACTN1
MYH10
PLEC
FLNA
JUP
```

**Extracellular matrix / remodeling**
```text
THBS2
FN1
TGFB1
MMP2
LAMA5
LAMB1
```

**Signalling**
```text
RAF1
AKT1
INSR
FZD5
TGFB1
```
So T1 is developing into a very recognizable:
> cell-adhesion / extracellular-matrix / cytoskeletal / signaling-associated transcriptional program.
And that fits our earlier T1 GO findings beautifully.

3. **This is stronger than simply counting GO terms**
```text
                    T1 ↔ T3
                       │
               DESeq2 ranked genes
                       │
             ┌─────────┴─────────┐
             ↓                   ↓
            GO                  KEGG
             ↓                   ↓
       enriched processes    pathways
             ↓                   ↓
             └─────────┬─────────┘
                       ↓
               shared core genes
                       ↓
          molecular interpretation
```

For T3:
```text
GO: chromosome segregation
             +
KEGG: cell cycle / DNA replication
             ↓
KNL1, NDC80, CDK1, BUB1,
CCNB1, RAD51, BRCA2...
```

For T1:
```text
GO: junction / epithelial / ECM
             +
KEGG: integrin / focal adhesion /
      cadherin / gap junction
             ↓
ITGA3, ITGB3, FN1, THBS2,
TLN1, VCL, DAG1, AKT1...
```

4. **Now the question is**
> Which genes are the most useful representatives of these two programs?

Not "which genes are significant?", because we already know that **none survived DESeq2 FDR < 0.05**.
Instead we want an **exploratory candidate set** supported by several pieces of evidence.

We can use four pieces:
   1. Present in GO core enrichment
   2. Present in KEGG core enrichment
   3. Direction agrees with the pathway side
   4. DESeq2 effect/statistic gives additional gene-level context

5. **Let's build the integrated candidate table**
We already have `shared_t3` and `shared_t1`.
For T3:

```r
t3_convergence <- merge(
    shared_t3_map,
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

t3_convergence <- t3_convergence[
    t3_convergence$log2FoldChange > 0,
]

t3_convergence <- t3_convergence[
    order(-abs(t3_convergence$stat)),
]

head(
    t3_convergence[
        ,
        c(
            "SYMBOL",
            "GENENAME",
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "padj",
            "stat"
        )
    ],
    20
)
```

For T1:
```r
t1_convergence <- merge(
    shared_t1_map,
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

t1_convergence <- t1_convergence[
    t1_convergence$log2FoldChange < 0,
]

t1_convergence <- t1_convergence[
    order(-abs(t1_convergence$stat)),
]

head(
    t1_convergence[
        ,
        c(
            "SYMBOL",
            "GENENAME",
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "padj",
            "stat"
        )
    ],
    20
)
```

6. **Why are we sorting by `abs(stat)`?**
The DESeq2 Wald statistic combines:
- estimated effect
- uncertainty around the estimate

So a larger absolute Wald statistic means the gene has stronger evidence for departure from zero within this model, even though it may still fail the final genome-wide FDR threshold.
We're therefore not simply saying:
> "largest fold change = most important."

We're incorporating the statistical evidence from the original DESeq2 model.


7. **The final candidate table will represent**
**T3 candidate module**
```text
KNL1
NDC80
CDK1
CCNB1
BUB1
BUB1B
TTK
RAD51
BRCA2
...
```

**T1 candidate module**
```text
ITGA3
ITGB3
FN1
THBS2
LAMA5
TLN1
DAG1
AKT1
VCL
...
```

8. **Comparing the approximate effect sizes we have seen**
T3:
```text
KNL1   +1.05
NUF2   +1.07
CENPF  +0.94
CDK1   +0.76
CCNB1  +0.68
```

T1:
```text
DVL1   -0.50
VLDLR  -0.78
NLGN3  -0.60
AMOT   -0.50
```

So the T3 proliferative genes appear to have a stronger and more uniformly positive gene-level shift than many of the T1 genes we've examined.
But we should not turn that observation into a formal comparison of biological strength without a prespecified statistical analysis. For now, it's descriptive.

Then we can move into:
```text
Candidate genes
      ↓
gene-level visualization
      ↓
pathway/gene network
      ↓
final biological interpretation
      ↓
limitations
      ↓
portfolio figures + README
```

1. **T3: a very coherent proliferative/genome-maintenance module**
Our top T3 GO+KEGG-convergent genes include:
**KNL1, SGO1, CCNE2, CDC25C, DBF4, IL1B, CDK1, CCNB1, CDC7, RAD51, RBBP8, NDC80, BUB1, NBN, BUB1B, RB1, MAD2L1BP, BLM, DNA2, BRCA2**. Their DESeq2 log2FC values are all positive, ranging from about +0.26 to +1.09 in these 20 rows.

We can group them into modules:

**Mitosis / kinetochore**
```text
KNL1
SGO1
NDC80
BUB1
BUB1B
MAD2L1BP
```

**Cell-cycle progression**
```text
CCNE2
CDC25C
CDK1
CCNB1
CDC7
DBF4
```

**DNA replication / repair**
```text
RAD51
RBBP8
NBN
BLM
DNA2
BRCA2
```

So this is no longer just "cell-cycle enrichment." We have **multiple molecular layers of a coordinated T3-associated proliferative/genome-maintenance program**.

2. **T1: the convergence gives us a better interpretation than GO alone**
Our T1 GO+KEGG-convergent genes are:
**GYS1, OXCT2, HEXA, GAA, INSR, HEXD, HRAS, PDGFB, PFKL, AKT2, VEGFB, COL2A1, CPT1C, SCD, AVPR1A, TNC, FHL2, PLEC, MAN2C1, ACADS.**

And all 20 have negative log2FC, from roughly −0.26 to −0.85.
This lets us refine the T1 story into two major components.

**Metabolic program**
> GYS1, OXCT2, HEXA, GAA, HEXD, PFKL, CPT1C, SCD, ACADS
which fits very nicely with our T1 GO/KEGG observations involving carbohydrate metabolism, small-molecule catabolism, glycan degradation, butanoate metabolism and related metabolic processes.

**Signaling / extracellular / tissue-organization program**
> INSR, HRAS, PDGFB, VEGFB, TNC, FHL2, PLEC
which fits the T1-associated signaling, extracellular and tissue-interaction picture.

And this is why we would not make the earlier "synaptic" interpretation a central conclusion. Synapse-associated GO terms were present, but the GO+KEGG convergence is much more clearly pointing toward a broader metabolic + extracellular/tissue/signaling state.
That's a useful example of why we integrate multiple analyses instead of interpreting one GO term in isolation.

3. **The most important statistical caveat**
Every one of these candidate genes still has:
> DESeq2 padj ≈ 0.869

For example:
```text
KNL1     log2FC +1.048   p 0.00365   padj 0.869
CDK1     log2FC +0.755   p 0.02038   padj 0.869

GYS1     log2FC -0.452   p 0.00343   padj 0.869
OXCT2    log2FC -0.569   p 0.00826   padj 0.869
```
So our conclusion must not be:
> "KNL1 and GYS1 are statistically significant DEGs."
They aren't.

Our conclusion is:
> These genes are exploratory candidates because they contribute to GSEA signals that are supported by both GO and KEGG, with gene-level direction consistent with the corresponding T3 or T1 state.

4. **Let's define our final candidate modules**
**T3 candidate module**
A representative set could be:
```text
KNL1
NDC80
CDK1
CCNB1
BUB1
BUB1B
CDC25C
RAD51
BRCA2
NBN
```
This covers:
```text
kinetochore
+
mitosis
+
cell-cycle control
+
DNA repair
```

**T1 candidate module**
```text
GYS1
HEXA
GAA
PFKL
SCD
INSR
PDGFB
VEGFB
TNC
FHL2
```
this covers:
```text
metabolism
+
metabolic regulation
+
signaling
+
extracellular/tissue biology
```

5. **Now making this reproducible**
```r
saveRDS(
    shared_t3,
    "results/candidates/shared_T3_GO_KEGG_genes.rds"
)

saveRDS(
    shared_t1,
    "results/candidates/shared_T1_GO_KEGG_genes.rds"
)

dir.create(
    "results/candidates",
    recursive = TRUE,
    showWarnings = FALSE
)

write.csv(
    t3_convergence,
    "results/candidates/T3_GO_KEGG_convergence_candidates.csv",
    row.names = FALSE
)

write.csv(
    t1_convergence,
    "results/candidates/T1_GO_KEGG_convergence_candidates.csv",
    row.names = FALSE
)
```

6. **But we should NOT immediately call these "final biomarkers"**
We have established so far:
```text
T1-associated genes
T3-associated genes
```
but we haven't yet visualized their sample-level behavior.

For a candidate such as KNL1 or GYS1, we want to know:
> Does the gene generally show the expected T1 → T3 shift across individual paired patients, or is the apparent effect being driven by a small subset of samples?
That's important because our PCA already suggested patient heterogeneity.

So now we're moving from:
> population-level statistics

to:
> individual patient-level behavior

7. **The next figure is a paired expression plot**
For selected representative candidates, we'll plot:
```text
Patient 3     T1 ───── T3
Patient 11    T1 ───── T3
Patient 13    T1 ───── T3
...
Patient 51    T1 ───── T3
```
with one line per patient.

For a strongly T3-associated gene:
```text
T1 ●
    ╲
     ● T3
```
for most patients.

For highly heterogenous gene:
```text
T1 ●──● T3
T1 ●──● T3
T1 ●
      ╲
       ● T3
...
```
This is particularly valuable because our DESeq2 model is **paired**:
> ~ patient_id + timepoint

















