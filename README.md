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



# Paired Transcriptomic Analysis of T1 vs T3

## Overview

This project analyzes paired transcriptomic profiles from **24 patients**, each contributing a matched **T1 and T3 sample** (**48 samples total**).

The primary biological question is:

> **What transcriptomic programs distinguish T3 from T1 after accounting for patient-specific variation?**

The analysis follows a reproducible RNA-seq workflow from expression preprocessing and exploratory analysis through paired differential expression, gene-set enrichment, pathway-level patient analysis, and candidate-gene validation.

---

## Main finding

The paired DESeq2 analysis identified **no individual genes reaching genome-wide false-discovery-rate significance**.

Rather than selecting nominally significant genes arbitrarily, the analysis therefore examined the **complete ranked transcriptome** using gene-set enrichment analysis.

The pathway-level results revealed coordinated biological structure:

### T3-associated programs

* Cell-cycle progression
* Chromosome segregation
* Kinetochore and spindle-checkpoint biology
* DNA replication
* DNA repair / homologous recombination
* Ribosome and protein-synthesis-related processes
* Selected immune-associated signaling

### T1-associated programs

* Carbohydrate and small-molecule metabolism
* Metabolic signaling
* Extracellular-matrix organization
* Cell-junction and adhesion biology
* Tissue/epithelial interaction programs

These findings suggest that the principal transcriptomic signal is distributed across **coordinated biological programs rather than a small set of individually genome-wide-significant genes**.

---

## Analytical workflow

```text
Matched T1/T3 samples
        ↓
Quality control
        ↓
Low-expression filtering
        ↓
DESeq2 normalization
        ↓
Variance-stabilizing transformation
        ↓
PCA / sample-structure analysis
        ↓
Paired differential expression
        ↓
Gene annotation
        ↓
Signed Wald-statistic ranking
        ↓
GO Biological Process GSEA
        ↓
GO redundancy reduction
        ↓
KEGG pathway GSEA
        ↓
GO–KEGG biological convergence
        ↓
Leading-edge gene analysis
        ↓
Representative candidate genes
        ↓
Patient-level paired validation
        ↓
Patient-level pathway scores
        ↓
Final figures + reproducibility audit
```

---

## Study design

The analysis uses a paired design because each patient contributes measurements at both T1 and T3.

The primary DESeq2 model is:

```text
~ patient_id + timepoint
```

with T1 as the reference level.

The primary contrast is:

```text
T3 vs T1
```

Including `patient_id` allows the analysis to account for patient-specific baseline differences.

---

## Key analysis results

### Preprocessing

| Metric                                    |  Result |
| ----------------------------------------- | ------: |
| Initial expression rows                   | ~19,503 |
| Rows after cleanup                        |  19,502 |
| Genes removed by low-expression filtering |   5,599 |
| Genes entering DESeq2                     |  13,903 |
| Patient pairs                             |      24 |
| Total samples                             |      48 |

### Exploratory structure

PCA showed substantial inter-patient heterogeneity:

* **PC1:** ~25% variance
* **PC2:** ~12% variance

T1 and T3 samples overlap substantially, reinforcing the importance of the paired model.

### Differential expression

| Criterion         | Genes |
| ----------------- | ----: |
| Nominal p < 0.05  |   718 |
| Nominal p < 0.01  |   115 |
| Nominal p < 0.001 |     9 |
| padj < 0.10       |     0 |
| padj < 0.05       |     0 |

The nominally significant genes were therefore **not treated as genome-wide significant discoveries**.

---

## Representative candidate genes

Candidate genes were selected from recurrent pathway-level signals and subsequently examined across individual patients.

### T3-associated

* **KNL1** — kinetochore / chromosome-segregation program
* **RAD51** — homologous recombination / DNA-repair program
* **CDK1** — cell-cycle regulation
* **CCNB1** — cell-cycle progression

### T1-associated

* **GYS1** — glycogen/carbohydrate metabolism
* **INSR** — metabolic signaling
* **SCD** — lipid metabolism
* **TNC** — extracellular/tissue-associated biology

These genes are **exploratory candidates**, not validated biomarkers and not genome-wide-significant DEGs.

---

## Patient-level validation

Selected candidates were examined using paired variance-stabilized expression differences:

```text
T3 - T1
```

Examples:

| Gene  | Expected direction | Patients showing expected direction | Paired Wilcoxon p |
| ----- | ------------------ | ----------------------------------: | ----------------: |
| KNL1  | T3 ↑               |                               17/24 |            0.0135 |
| RAD51 | T3 ↑               |                               16/24 |            0.0891 |
| GYS1  | T1 ↑               |                               17/24 |           0.00894 |
| INSR  | T1 ↑               |                               17/24 |            0.0135 |

These tests are exploratory because candidate selection was informed by the same dataset.

---

## Repository structure

```text
breast-cancer-rnaseq-transcriptomics/
│
├── data/
│   ├── raw/
│   └── processed/
│
├── scripts/
│   ├── 01_extract_T1_T3.sh
│   ├── 02_qc_filter_normalize.R
│   ├── 03_deseq2.R
│   ├── 04_diagnose_DE.R
│   ├── 05_annotation_and_ranking.R
│   ├── 06_gsea_and_pathways.R
│   ├── 07_final_candidate_validation_and_figures.R
│   ├── 08_final_figures.R
│   ├── 09_final_audit.R
│   ├── 10_patient_pathway_scores.R
│   └── 11_leading_edge_heatmap.R
│
├── notebooks/
│   ├── 01_project_walkthrough.ipynb
│   └── 02_final_results.ipynb
│
├── results/
│   ├── differential_expression/
│   ├── enrichment/
│   ├── candidates/
│   └── figures/
│
├── docs/
│   ├── PROJECT_REPORT.md
│   ├── METHODS.md
│   └── LIMITATIONS.md
│
└── README.md
```

### Directory roles

**`scripts/`**
Contains the executable analysis pipeline. Each script performs a defined computational step and writes its outputs to the repository.

**`notebooks/`**
Contains the human-readable reasoning, selected code, visualizations and interpretation behind the analysis.

**`results/`**
Contains generated tables, pathway results, candidate results and publication-style figures.

**`docs/`**
Contains the formal scientific documentation, including methods, final interpretation and limitations.

**`data/`**
Contains the source and processed data required by the pipeline.

---

## Reproducibility

The analysis is designed so that the major results can be regenerated from the saved inputs and scripts.

The main computational sequence is:

```bash
01_extract_T1_T3.sh
02_qc_filter_normalize.R
03_deseq2.R
04_diagnose_DE.R
05_annotation_and_ranking.R
06_gsea_and_pathways.R
07_final_candidate_validation_and_figures.R
08_final_figures.R
10_patient_pathway_scores.R
11_leading_edge_heatmap.R
09_final_audit.R
```

The final audit is intentionally run last so that missing files, unexpected sample structure, and inconsistencies in the primary statistical conclusion can be detected before the project is considered complete.

---

## Main figures

| Figure   | Description                                   |
| -------- | --------------------------------------------- |
| Figure 1 | PCA of paired T1/T3 samples                   |
| Figure 2 | T3 vs T1 differential-expression volcano plot |
| Figure 3 | GO Biological Process GSEA                    |
| Figure 4 | KEGG pathway GSEA                             |
| Figure 5 | Candidate-gene paired expression              |
| Figure 6 | Candidate-gene heatmap                        |
| Figure 7 | Patient-level pathway activity scores         |
| Figure 8 | T3 leading-edge gene heatmap                  |

---

## Statistical interpretation

A central principle of this project is the distinction between:

```text
Individual-gene evidence
        versus
Pathway-level evidence
```

No individual gene survived genome-wide FDR correction.

However, coordinated shifts across groups of biologically related genes were detected by ranked gene-set analysis.

Therefore, the principal conclusions are expressed at the **pathway/program level**, while individual genes are treated as exploratory representatives.

---

## Important limitations

This project is hypothesis-generating.

In particular:

1. The absence of genome-wide-significant individual genes limits claims about specific transcripts.
2. GSEA identifies coordinated statistical enrichment but does not prove biological causality.
3. GO terms are redundant and should not be interpreted as independent discoveries.
4. Disease-labelled KEGG pathways do not imply that the corresponding disease is present in the study subjects.
5. Candidate genes were selected from the same dataset used for validation.
6. GO/KEGG convergence represents internal concordance rather than independent replication.
7. RNA expression does not directly establish protein abundance or functional activity.
8. External cohort or experimental validation would be required for definitive biomarker or mechanistic claims.

See [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md) for the full discussion.

---

## Documentation

For a detailed description of the analysis:

* [`docs/METHODS.md`](docs/METHODS.md)
* [`docs/PROJECT_REPORT.md`](docs/PROJECT_REPORT.md)
* [`docs/LIMITATIONS.md`](docs/LIMITATIONS.md)

For the analysis narrative:

* [`notebooks/01_project_walkthrough.ipynb`](notebooks/01_project_walkthrough.ipynb)
* [`notebooks/02_final_results.ipynb`](notebooks/02_final_results.ipynb)

---

## Tools

The project uses R/Bioconductor-based transcriptomic analysis, including:

* DESeq2
* GSVA / ssGSEA
* clusterProfiler
* org.Hs.eg.db
* AnnotationDbi
* ggplot2
* pheatmap

---

## Project status

**Analysis:** Complete
**Figures:** Complete
**Documentation:** Complete
**Reproducibility audit:** Final stage
**Biological interpretation:** Hypothesis-generating

---

## Final takeaway

The major result of this project is not a list of significant genes.

It is the identification of a coordinated transcriptomic shift from a T1-associated
metabolic/tissue-interaction state toward a T3-associated proliferative and
genome-maintenance program.

The project is designed to demonstrate not only the ability to run
bioinformatics tools, but also the ability to:

* choose an appropriate statistical design,
* recognize the limitations of individual-gene testing,
* move to pathway-level analysis when justified,
* connect pathway results back to genes and individual patients,
* distinguish exploratory evidence from validated findings,
* and maintain a reproducible computational workflow.


























