# Baseline Transcriptomic Profiling of Breast Cancer Response to Neoadjuvant Chemotherapy

## Dataset
GSE122630

## Organism
Homo sapiens

## Platform
Ion Torrent Proton RNA-seq / transcriptome sequencing.

## Samples we'll use
50 baseline T1 samples
   - 16 responders
   - 34 non-responders
(source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122630)

## Main question
How does chemotherapy alter the breast-cancer transcriptome from before treatment (T1) to mid-chemotherapy (T3)?

## Final pipeline

GSE122630
    ↓
Download expression data + metadata
    ↓
Select baseline T1 samples
    ↓
Verify 16 R / 34 NR
    ↓
Explore count matrix
    ↓
Filter low-expression genes
    ↓
Normalization
    ↓
PCA
    ↓
Differential expression
    ├── log2FC
    ├── p-value
    └── adjusted p-value
    ↓
Volcano plot
    ↓
Heatmap
    ↓
GO enrichment
    ↓
KEGG/pathway analysis
    ↓
Candidate response-associated genes
    ↓
Biological interpretation
    ↓
GitHub README + figures



## PCA and Exploratory QC

Variance-stabilized expression values were used for principal component
analysis.

PC1 explained 25% of the total variance and PC2 explained 12%, giving
37% cumulative variance across the first two principal components.

The T1 and T3 samples showed substantial overlap in PC1-PC2 space rather
than a strong global separation by timepoint. The paired PCA also showed
heterogeneous movement between T1 and T3 across patients, indicating
substantial inter-patient transcriptomic variation.

PCA was used as an exploratory QC and visualization method and was not
used to determine differential-expression significance or to exclude
samples solely on the basis of their position in PCA space.


## Reproducibility

The analysis is organized into numbered Bash and R scripts under `scripts/`.

The scripts are designed to recreate the processed data, quality-control
results, normalization, PCA, and differential-expression analysis from
the project input files.

The R session itself is not treated as the permanent record of the
analysis. Intermediate DESeq2 objects are saved as `.rds` files so that
later stages can be resumed without manually reconstructing the session.


### Differential Expression

A paired DESeq2 model was fitted using:

```text
~ patient_id + timepoint
```

























