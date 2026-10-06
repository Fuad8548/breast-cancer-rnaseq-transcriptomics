## 1. Install the annotation/enrichment packages

```r
BiocManager::install(
    c("org.Hs.eg.db", "clusterProfiler")
)
```
then:
```r
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
```

## 2. Reload our fitted DESeq2 model
```r
dds <- readRDS(
    "r_objects/dds_T1_T3_fitted.rds"
)
```

then recreate exactly our comparison:
```r
res <- results(
    dds,
    contrast = c("timepoint", "T3", "T1")
)
```
Add the gene IDs:
```r
res_df <- as.data.frame(res)

res_df$gene_id <- rownames(res_df)
```


## 3. Annotate the Ensembl IDs
```r
annotation <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = rownames(res),
    columns = c(
        "SYMBOL",
        "ENTREZID",
        "GENENAME"
    ),
    keytype = "ENSEMBL"
)
```
Now inspect:
```r
head(annotation)
```

Output at: [Enrichment](/results/enrichment/ensembl_annotation_map.csv)


## 4. Check our annotation success rate

```r
sum(!is.na(annotation$SYMBOL))
```

and:
```r
sum(!is.na(annotation$ENTREZID))
```

Also:
```r
length(unique(annotation$ENSEMBL))
```

The numbers tell us how much of our 13,903-gene result can be annotated.

Save the mapping:
```r
dir.create(
    "results/enrichment",
    recursive = TRUE,
    showWarnings = FALSE
)

write.csv(
    annotation,
    "results/enrichment/ensembl_annotation_map.csv",
    row.names = FALSE
)
```

## 🧠 Why do we care about Entrez IDs?
> Different biological databases speak different identifier languages.

- Our expression matrix: Ensembl
- DESeq2: Ensembl
- Some enrichment resources: Entrez
So part of bioinformatics is **identifier translation**.

`clusterProfiler` explicitly provides ID-conversion tools, and its GO/KEGG GSEA functions work with ranked gene lists.


## 5. why GSEA?
Think of this:
```text
Strongest T3 evidence
        ↑
Gene
Gene
Gene
Gene
...
Gene
Gene
Gene
        ↓
Strongest T1 evidence
```

GSEA then asks:
> Do genes belonging to a particular biological pathway tend to occur disproportionately toward one end of the ranking?

The current `clusterProfiler` documentation describes GSEA as operating on an ordered ranked gene list.
This is exactly the kind of analysis that makes sense for our current result.

## 6. What should we use for the ranking?
We'll use the **DESeq2 Wald statistic**, `stat`, because it is signed. (Not `p-value`, not `log2FC` alone)

So: 
```text
positive statistic
        ↓
evidence toward T3
```

while:
```text
negative statistic
        ↓
evidence toward T1
```
And its magnitude incorporates both the estimated effect and its uncertainty.
That's much more informative for ranking than simply sorting by fold change.


## 7. Build our ranked list

First combine the statistics with the annotation:
```r
rank_df <- data.frame(
    gene_id = rownames(res),
    stat = res$stat,
    stringsAsFactors = FALSE
)

rank_df <- merge(
    rank_df,
    annotation[, c("ENSEMBL", "ENTREZID")],
    by.x = "gene_id",
    by.y = "ENSEMBL"
)
```

Remove genes without a usable statistic or Entrez ID:
```r 
rank_df <- rank_df[
    !is.na(rank_df$stat) &
    !is.na(rank_df$ENTREZID),
]
```

## 8. Deal with duplicate Entrez mappings
One Entrez ID can sometimes correspond to more than one Ensembl record in an annotation mapping.
GSEA expects a ranked vector with unique identifiers, so we need a transparent rule.
We'll retain, for each Entrez ID, the Ensembl record with the strongest absolute Wald statistic.

```r
rank_df <- rank_df[
    order(
        abs(rank_df$stat),
        decreasing = TRUE
    ),
]

rank_df <- rank_df[
    !duplicated(rank_df$ENTREZID),
]
```

Then construct the actual ranked vector:
```r
gene_list <- rank_df$stat

names(gene_list) <- as.character(
    rank_df$ENTREZID
)

gene_list <- sort(
    gene_list,
    decreasing = TRUE
)
```

Check:
```r
head(gene_list)
```
and:
```r
tail(gene_list)
```

Output at: [Ranked Gene List](/results/enrichment/ranked_gene_list_for_GSEA.csv)

















