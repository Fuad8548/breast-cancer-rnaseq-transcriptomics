We are not starting the RNA-seq analysis again. We already have the ranked gene list saved from Stage 5. We are simply rebuilding the GSEA object because R sessions are temporary.

Our DESeq2 analysis asked:
> Which individual genes differ between T3 and T1?

GSEA asks a different question:
> Do groups of biologically related genes tend to occur toward the T3 side or T1 side of the entire ranked transcriptome?

1. **Loading the tools and our saved ranked genes**
```r
library(clusterProfiler)
library(org.Hs.eg.db)

gene_list <- readRDS(
    "r_objects/gene_list_T3_vs_T1_Entrez.rds"
)

length(gene_list)
head(gene_list)
```

We previously ended up with approximately 12,665 unique Entrez genes. 
The important point is that `gene_list` is already ranked:

```text
positive values → T3-associated
negative values → T1-associated
```
and we are using the **DESeq2 Wald statistic** for that ranking.

2. **Run GO Biological Process GSEA again**

```r
gsea_go_bp <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE
)
```

Conceptually, this does:
```text
12,665 ranked genes
        ↓
GO database
        ↓
groups genes into biological processes
        ↓
asks whether each process is concentrated
toward T3 or T1
        ↓
calculates enrichment statistics
        ↓
corrects for many GO-term tests
```

The result is stored in:
```r
gsea_go_bp
```

3. **Verify that we recovered our previous analysis**

```r
go_df <- as.data.frame(gsea_go_bp)

nrow(go_df)
```

Then:
```r
sum(go_df$p.adjust < 0.05, na.rm = TRUE)
```
We previously obtained: **863**
significant GO Biological Process terms.

## Save the GSEA result

```r
dir.create(
    "results/enrichment",
    recursive = TRUE,
    showWarnings = FALSE
)

saveRDS(
    gsea_go_bp,
    "results/enrichment/gsea_GO_BP_T3_vs_T1.rds"
)

write.csv(
    go_df,
    "results/enrichment/GO_BP_GSEA_T3_vs_T1_all_terms.csv",
    row.names = FALSE
)
```







