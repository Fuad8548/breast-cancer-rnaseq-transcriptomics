# Building the volcano plot
A volcano plot simultaneously shows:
- x-axis: `log2FoldChange`
- y-axis: `-log10(p-value)`

So:
- A point far to the right means a large positive effect.
- A point far to the left means a large negative effect.
- A point high up means a small p-value.

First create a plotting dataframe:
```r
volcano_df <- as.data.frame(res)

volcano_df$gene_id <- rownames(volcano_df)

volcano_df$neg_log10_p <- -log10(volcano_df$pvalue)

volcano_df$significant <- volcano_df$padj < 0.05
```

Now make the plot:
```r 
library(ggplot2)

volcano <- ggplot(
    volcano_df,
    aes(
        x = log2FoldChange,
        y = neg_log10_p
    )
) +
    geom_point(
        alpha = 0.6
    ) +
    geom_vline(
        xintercept = c(-1, 1),
        linetype = "dashed"
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    ) +
    labs(
        title = "Volcano Plot: T3 vs T1",
        x = "log2 fold change",
        y = "-log10(p-value)"
    ) +
    theme_minimal()
```

Display:
```r
volcano
```

Save:
```r
ggsave(
    "results/differential_expression/volcano_T3_vs_T1.pdf",
    volcano,
    width = 8,
    height = 6
)
```

**One thing to notice:**
We are deliberately using the raw p-value on the y-axis for the classic volcano plot.
That lets us see the underlying signal:
```text
718 genes with p < 0.05
115 with p < 0.01
9 with p < 0.001
```
But our interpretation remains:
FDR < 0.05
    ↓
0 genes



# How to read our volcano plot

**Left vs Right**
← T3 lower              T3 higher →
   expression          expression

     -log2FC     0   +log2FC

So:
- **left side** = genes estimated to be lower at T3 than T1
- **right side** = genes estimated to be higher at T3 than T1

**Up vs down**

Higher on the graph means:
    −log10​(p)
is larger, which means the raw p-value is smaller.

So the upper corners are where you'd expect visually compelling candidates:

                  stronger evidence
                         ↑
                 ●       │       ●
              ●          │          ●
                         │
 lower T3  ←─────────────┼──────────────→ higher T3
                         │


## What our plot actually shows
### 1. Most genes are tightly concentrated around log2FC = 0
There's a huge dense "V/funnel" centered around:
> log2FC ≈ 0

That means most genes have relatively modest estimated changes between T1 and T3.
That agrees with what we already found numerically:
> maximum absolute log2FC ≈ 2.09.

So there aren't enormous transcriptome-wide shifts across thousands of genes.

### 2. There IS some nominal signal
We can clearly see a substantial population of points rising above the horizontal line:
> p=0.05

That's consistent with your earlier calculation:
```text
p < 0.05       718 genes
p < 0.01       115 genes
p < 0.001        9 genes
```

### 3. There are genes with appreciable effect sizes
Our vertical dashed lines are:
```text
log2FC = -1
log2FC = +1
```
which correspond roughly to:
```text
-1 → half as much expression
+1 → twice as much expression
```
Our plot has some points beyond these boundaries, especially on the positive side, and our numerical analysis found a maximum absolute log2FC of **2.093**.

### 4. And this is the main lesson from your volcano
> 718 genes showed nominal p < 0.05, but no genes reached FDR < 0.05 in the paired T3-vs-T1 analysis.

### Why the horizontal line can be misleading
Our dashed horizontal line corresponds to:
> p=0.05

But we tested 13,903 genes.
If we use 0.05 independently for thousands of hypotheses, some genes will cross that line just by chance.

That's why we use the adjusted p-value:
```text
raw p-value
     ↓
multiple-testing correction
     ↓
padj / FDR
```
And our actual result was:
> FDR < 0.05 → 0 genes
So the genes above that horizontal line are nominal candidates, not confirmed DEGs.

























