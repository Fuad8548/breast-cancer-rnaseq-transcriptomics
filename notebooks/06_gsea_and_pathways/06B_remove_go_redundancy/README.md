# Remove GO redundancy

`simplify()`asks:
> Which highly similar GO terms are redundant enough that we can retain a representative term?

```r
gsea_go_simplified <- simplify(
    gsea_go_bp,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
)
```

Then:
```r
go_simple_df <- as.data.frame(gsea_go_simplified)

nrow(go_simple_df)
```

### What does cutoff = 0.7 mean?
It is a similarity threshold.

GO term A ─────────────── GO term B
             similarity

When terms are sufficiently similar, `simplify()` considers them redundant and retains a representative according to the criterion we specify. This is not changing the underlying expression data or DESeq2 statistics. It is only making the biological interpretation cleaner.

### Why by = "`p.adjust`" and `select_fun = min`?
```r
by = "p.adjust"
```
to compare terms using their adjusted significance.

And:
```r
select_fun = min
```
means that among redundant terms, the representative is selected using the smallest adjusted p-value.
So we're basically saying:
> "When several highly similar GO descriptions are telling essentially the same story, retain the statistically strongest representative."

### Save the simplified result too
```r
saveRDS(
    gsea_go_simplified,
    "results/enrichment/gsea_GO_BP_T3_vs_T1_simplified.rds"
)

write.csv(
    go_simple_df,
    "results/enrichment/GO_BP_GSEA_T3_vs_T1_simplified.csv",
    row.names = FALSE
)
```

Now we have both:
```text
Original GSEA
       ↓
all GO terms
       ↓
redundancy reduction
       ↓
simplified GO terms
```

## Inspect what survived after redundancy reduction
```r
head(
    go_simple_df[
        order(go_simple_df$p.adjust),
        c("ID", "Description", "NES", "p.adjust")
    ],
    30
)
```

### 1. T3-associated programs
The strongest group is still overwhelmingly **cell division** and **genome maintenance**:

| Theme                  | Your GO terms                                                                                  |                 NES |
| ---------------------- | ---------------------------------------------------------------------------------------------- | ------------------: |
| Chromosome/mitosis     | spindle microtubule attachment to kinetochore                                                  |               +2.76 |
| Nuclear division       | nuclear division                                                                               |               +2.70 |
| Chromosome segregation | negative regulation of chromosome segregation; chromosome localization; chromosome segregation | +2.69, +2.66, +2.63 |
| DNA replication        | DNA-templated DNA replication; regulation of DNA replication                                   |        +2.65, +2.48 |
| Spindle checkpoint     | mitotic spindle assembly/checkpoint signaling                                                  |               +2.65 |
| Cell cycle             | positive regulation of cell-cycle process; mitotic cell-cycle transition                       |        +2.64, +2.23 |
| DNA repair             | homologous recombination double-strand-break repair                                            |               +2.49 |
| Chromatin              | nucleosome organization                                                                        |               +2.42 |


So we can reasonably describe the dominant T3-associated transcriptional program as:
> cell-cycle progression, chromosome segregation, mitotic spindle/kinetochore activity, DNA replication, DNA repair, and chromatin organization.

Notice how much more informative this is than our original "0 significant DEGs."


### 2. There is also a T3-associated immune signal

```text
adaptive immune response                NES +2.284
lymphocyte differentiation              NES +2.101
immune receptor signaling               NES +2.100
regulation of lymphocyte activation     NES +1.981
T cell activation                       NES +1.978
immune effector process                 NES +1.913
leukocyte activation                    NES +2.114
```

So among the simplified significant terms, there is a second major T3-associated theme:
> immune activation / lymphocyte-related biology.

But notice that these are pathway-level associations, not evidence that a specific immune cell population increased. We would need additional analyses, such as cell-type deconvolution or marker analysis, before making that stronger claim.

### 3. What about T1?
The most obvious negative-NES terms in our simplified top 30 are:

```text
GO:0034329  cell junction assembly       NES -1.956
GO:0050808  synapse organization         NES -1.946
```

And from the unsimplified result you previously showed, T1 also had:

```text
carbohydrate metabolism
lipid metabolic regulation
epithelial morphogenesis
extracellular matrix organization
vasculature development
blood vessel development
```

So the T1 side appears to contain processes related to:
- **cellular/tissue organization + metabolism + extracellular/tissue structure.**

However, we should be careful here because we haven't yet inspected the genes responsible for each pathway.

### 4. There is an important subtlety in our result

```text
T3 → cell junction assembly? No
T1 → cell junction assembly → NES -1.956
```
and:
```text
T1 → synapse organization → NES -1.946
```
These "synapse" terms shouldn't immediately make us say:
> "The breast tissue became neuronal at T1."

GO pathways can contain genes that participate in broader processes such as:
- cell adhesion
- cytoskeletal organization
- membrane organization
- cell-cell interaction

which may also be annotated under synaptic processes. The next step therefore has to be **gene-level inspection**.

### 5. This is where core_enrichment becomes important
A GSEA result contains a field called:
```r
core_enrichment
```
This is essentially the set of genes that are contributing most directly to the enrichment signal for that gene set, often called the **leading-edge/core-enrichment genes**.

For example, suppose:
```text
DNA replication
core_enrichment:
GENE_A/GENE_B/GENE_C/GENE_D/...
```

We can now ask:
> Which actual genes are producing the strong +2.65 NES?

### 6. Let's inspect it for the major T3 program
Run this:
```r
t3_top <- go_simple_df[
    go_simple_df$NES > 0 & go_simple_df$p.adjust < 0.05,
]

t3_top <- t3_top[
    order(t3_top$p.adjust),
]

t3_top[
    1:min(10, nrow(t3_top)),
    c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]
```

and separately for T1:
```r
t1_top <- go_simple_df[
    go_simple_df$NES < 0 & go_simple_df$p.adjust < 0.05,
]

t1_top <- t1_top[
    order(t1_top$p.adjust),
]

t1_top[
    1:min(10, nrow(t1_top)),
    c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]
```

This takes us from:
```text
GO term
   ↓
biological interpretation
```

to:
```text
GO term
   ↓
core-enrichment genes
   ↓
specific molecular players
```

### 7. And this will connect back to our DESeq2 results
Remember:
```text
DESeq2
→ 0 genes at FDR < 0.05
```

Yet GSEA tells us:
```text
T3
→ mitosis
→ chromosome segregation
→ DNA replication
→ DNA repair
→ immune programs
```

Now we'll inspect the core genes.
Some may have:
```text
moderate log2FC
moderate p-value
```
But collectively they generate a powerful pathway-level signal.

That gives us a scientifically defensible explanation for why **gene-level significance** and **pathway-level significance** don't have to agree.
















