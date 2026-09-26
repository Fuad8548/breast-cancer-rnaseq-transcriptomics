## Step 1: Step 1: Why did we make a ranked gene list?
This was the purpose of:
```r
gene_list <- readRDS(
    "data/processed/gene_list_T3_vs_T1_Entrez.rds"
)
```

Our ranked list looks conceptually like:

```bash
Gene        Wald statistic
──────────────────────────
Gene A          +3.67
Gene B          +3.61
Gene C          +3.54
Gene D          +3.42
...
Gene X           0.01
Gene Y          -0.02
...
Gene Z          -3.67
```

We decided earlier that:
- positive Wald statistic = T3-associated
- negative Wald statistic = T1-associated

So we've transformed our DESeq2 result into a **continuous ranking of the entire transcriptome**.

## Step 2: What is GSEA actually asking?
Suppose we have this ranked list:

```r
T3 side
 ↑
 │  Gene A
 │  Gene B
 │  Gene C
 │  Gene D
 │
 │  Gene E
 │  Gene F
 │
 │  Gene G
 │
 │  Gene H
 │
 │  Gene I
 │
 │  Gene J
 ↓
T1 side
```
Now suppose we have a biological pathway: "**DNA replication**"
and its genes are:
```bash
DNA replication genes:
A
C
D
G
```

Look where they appear:
```bash
A  ← pathway gene
B
C  ← pathway gene
D  ← pathway gene
E
F
G  ← pathway gene
H
I
J
```
They're disproportionately concentrated toward the **T3 end**.

GSEA asks:
> Are the genes belonging to this biological pathway unusually concentrated toward one end of our ranked transcriptome?

## Step 3: What does gseGO() mean?

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

**Explanation:**

1. `gseGO()` means:
> Perform Gene Set Enrichment Analysis using Gene Ontology.

We are asking:
> "Which GO biological processes have their genes disproportionately concentrated toward the T3 or T1 end of my ranked transcriptome?"

2. `geneList = gene_list`
This tells GSEA:
> "Here is my ranked transcriptome."

Our ranking comes from the DESeq2 Wald statistic. So, we are connecting:

```text
DESeq2
   ↓
Wald statistics
   ↓
ranked genes
   ↓
GSEA
```

3. `OrgDb = org.Hs.eg.db`
This is basically our **human gene annotation database**.

GSEA needs to know things like:
```text
ENTREZ 8115 → which gene?
                  ↓
             which GO terms?
                  ↓
          which biological processes?
```
Without an annotation database, the software wouldn't know which genes belong to which biological processes.

4. `keyType = "ENTREZID"`

This tells clusterProfiler:
> "The identifiers in my ranked list are Entrez Gene IDs."

Remember that our original data used Ensembl IDs. We previously performed:

```text
Ensembl ID
    ↓
gene annotation
    ↓
Entrez ID
    ↓
ranked list
```
because the pathway databases we're using need compatible gene identifiers.

5. `ont = "BP"`
GO has three major branches:

```text
Gene Ontology
│
├── BP = Biological Process
├── MF = Molecular Function
└── CC = Cellular Component
```

We're currently asking specifically about:
> Biological Processes

So we're interested in things like:
- DNA replication
- chromosome segregation
- cell division
- extracellular matrix organization
- carbohydrate metabolism

## Step 4: Why minGSSize and maxGSSize?

```r
minGSSize = 10
maxGSSize = 500
```
A gene set is simply a group of genes belonging to a biological category.

For example:
```text
DNA replication
    ↓
Gene A
Gene B
Gene C
...
Gene 80
```
We are saying:
> Don't analyze extremely tiny gene sets containing fewer than 10 genes.

and:
> Don't analyze enormous gene sets containing more than 500 genes.

This prevents very tiny or extremely broad categories from dominating the analysis.

## Step 5: Why pAdjustMethod = "BH"?
This connects directly to our earlier DESeq2 problem.
We're testing **many GO terms**.
Imagine:
```text
GO term 1
GO term 2
GO term 3
...
GO term 5,000
```
If we use ordinary p-values, some pathways will appear significant purely by chance. So we apply Benjamini-Hochberg multiple-testing correction.
The resulting: `p.adjust` is the value we care about when determining whether a pathway survives multiple-testing correction.

## Step 6: What did GSEA actually tell us?
For example:
```bash
mitotic sister chromatid segregation

NES = +2.898
p.adjust = 8.04 × 10⁻⁹
```
This means: **Positive NES**

The genes belonging to this biological process are disproportionately concentrated toward the **T3-associated end** of the ranked list.

**Large absolute NES** => the enrichment is relatively strong

**Extremely small adjusted p-value**
The enrichment is unlikely to be explained simply by random ranking, even after correcting for testing many GO terms.
Therefore:
> The transcriptome contains a coordinated T3-associated signal involving mitotic sister chromatid segregation.

## Step 7: Why did we get 863 significant GO terms?

GO is hierarchical. For example, these aren't necessarily 20 independent discoveries:
```text
             Cell division
                  │
          Nuclear division
                  │
        Chromosome segregation
                  │
       Sister chromatid segregation
                  │
       Mitotic sister chromatid
             segregation
```
The same genes can appear in many related GO terms.

Therefore:
> 863 significant GO terms does not mean 863 independent biological discoveries.
This is why we introduced `simplify()`.

## Step 8: What does simplify() actually do?

```r
gsea_go_simplified <- simplify(
    gsea_go_bp,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
)
```

It is essentially asking:
> "Among highly similar GO terms, can we keep representative terms rather than showing hundreds of redundant descriptions?"

For example, imagine:
```text
chromosome separation
sister chromatid separation
chromosome segregation
nuclear chromosome segregation
mitotic chromosome segregation
```
They may contain many of the same genes. A redundancy-reduction step makes the final biological interpretation easier.

## Step 9: And this is where our analysis becomes scientifically interesting

**In our entire project, individual-gene analysis**

```bash
0 genes
FDR < 0.05
```
In pathway level analysis, we found strong coordinated signals:

T3:
```text
mitosis
chromosome segregation
kinetochore
spindle checkpoint
DNA replication
```

T1:
```text
cell junction
epithelial organization
extracellular matrix
vascular development
metabolism
```

So our scientific story is evolving from:
> "There are no significant genes."

to something much more nuanced:
> "Although no individual gene passed genome-wide FDR < 0.05 in the paired T3-vs-T1 analysis, ranked gene-set analysis revealed coordinated transcriptional programs associated with the two timepoints."

## Final Codeblock:

```r
library(clusterProfiler)
library(org.Hs.eg.db)

gene_list <- readRDS(
    "data/processed/gene_list_T3_vs_T1_Entrez.rds"
)

length(gene_list)
head(gene_list)
```

Output:
```bash
12665
```
because we previously had 12,665 unique Entrez genes in the final ranked list.

Then run:
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

Then verify:
```r
nrow(as.data.frame(gsea_go_bp))
```

and:
```r
sum(as.data.frame(gsea_go_bp)$p.adjust < 0.05, na.rm = TRUE)
```
We should recover approximately the **863 significant GO BP terms** you obtained previously.

Then run:
```r
gsea_go_simplified <- simplify(
    gsea_go_bp,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
)
```

then:
```r
go_simple_df <- as.data.frame(gsea_go_simplified)

nrow(go_simple_df)
```






































