# Stage 6C: Remove GO redundancy

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

### Why by = "p.adjust" and select_fun = min?
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

### Inspect what survived after redundancy reduction
```r
head(
    go_simple_df[
        order(go_simple_df$p.adjust),
        c("ID", "Description", "NES", "p.adjust")
    ],
    30
)
```

## T3-associated programs
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
























