# Project 1 — Methods

## Study design

This analysis used processed transcriptomic data from GEO dataset GSE122630.

The selected analysis cohort consisted of 24 patients with matched T1 and T3
samples, producing 48 samples in total.

Because each patient contributed both timepoints, the comparison was treated as
a paired design.

The primary statistical model was:

```text
~ patient_id + timepoint
```

with T1 as the reference level and the coefficient
`timepoint_T3_vs_T1` representing the T3-versus-T1 comparison after accounting
for patient-specific effects.

## Expression preprocessing

The processed expression matrix initially contained approximately 19,503
rows.

After removal of an invalid/duplicate feature, 19,502 genes remained.

Low-expression filtering removed 5,599 genes, leaving 13,903 genes for the
DESeq2 analysis.

Filtering was used to reduce noise and avoid spending statistical power on
features with insufficient expression evidence.

## Normalization and exploratory analysis

DESeq2 size-factor normalization was applied to account for differences in
library size and composition.

The estimated normalization factors ranged from approximately 0.555 to 2.075,
with a median near 1.

Variance-stabilizing transformation was subsequently used for exploratory
visualization.

Principal component analysis was performed to inspect dominant sources of
variation and evaluate the sample structure.

## Differential expression

Differential expression was performed with DESeq2 using its negative-binomial
model.

The primary contrast was:

```text
T3 vs T1
```

The significance of individual genes was assessed using Benjamini-Hochberg
false-discovery-rate adjustment.

The analysis included 13,903 genes. No gene reached an adjusted p-value below
0.10 or 0.05.

## Ranked gene-set analysis

Because no individual genes passed genome-wide FDR correction, pathway-level
analysis was performed using the complete ranked transcriptome.

Genes were ranked by the signed DESeq2 Wald statistic.

Positive values indicate association with T3; negative values indicate
association with T1.

This approach avoids imposing an arbitrary individual-gene significance cutoff
before pathway analysis.

## GO enrichment

Gene Ontology Biological Process gene-set enrichment analysis was performed
using the ranked gene list.

Multiple-testing correction was applied to pathway-level results.

Because GO terms are hierarchical and overlapping, redundant terms were
reduced using semantic similarity-based simplification.

## KEGG enrichment

KEGG pathway enrichment was performed using the same ranked transcriptomic
signal.

GO and KEGG were interpreted as complementary representations of the same
underlying expression data rather than independent biological experiments.

## Candidate selection

Representative candidate genes were selected from recurrent pathway-level
signals.

T3-associated representatives included:

* KNL1
* RAD51
* CDK1
* CCNB1

T1-associated representatives included:

* GYS1
* SCD
* INSR
* TNC

Candidate status does not imply genome-wide significance.

## Patient-level validation

For each candidate, variance-stabilized expression was summarized within each
patient using:

```text
T3 - T1
```

paired differences.

The analysis examined:

* median paired difference,
* mean paired difference,
* number and percentage of patients showing the expected direction,
* paired Wilcoxon signed-rank testing.

These analyses were treated as exploratory because the genes were selected after
inspection of the pathway-level results.

## Reproducibility

The complete analysis is divided into executable scripts, notebooks, generated
results and documentation.

The scripts contain the reproducible computational workflow, whereas the
notebooks provide the reasoning and interpretation behind the analyses.
Generated figures and tables are stored in the results directories.
