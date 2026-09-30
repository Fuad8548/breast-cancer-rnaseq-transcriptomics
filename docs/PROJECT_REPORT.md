# Project 1 — Paired Transcriptomic Analysis of T1 vs T3 Samples

## 1. Overview

This project investigates transcriptomic differences between matched T1 and T3
samples from 24 patients using a paired RNA-seq analysis framework.

The central question was:

> What molecular programs distinguish T3 from T1 after accounting for
> patient-specific baseline variation?

The analysis progressed from preprocessing and exploratory visualization to
paired differential expression, ranked gene-set enrichment, GO/KEGG
convergence and patient-level candidate validation.

## 2. Study design

The final cohort contained:

* 24 patients
* 24 T1 samples
* 24 T3 samples
* 48 total samples

Because T1 and T3 samples were matched within patients, the primary model
included patient identity as a blocking factor:

```text
~ patient_id + timepoint
```

The principal contrast was T3 versus T1.

## 3. Expression preprocessing

The processed matrix contained approximately 19,503 genes before filtering.

After removal of an invalid/duplicate feature, 19,502 rows remained.

Low-expression filtering removed 5,599 genes, leaving 13,903 genes for
differential-expression analysis.

DESeq2 normalization produced size factors ranging approximately from 0.555 to
2.075.

## 4. Exploratory structure

Variance-stabilized PCA showed that PC1 explained approximately 25% of the total
variance and PC2 approximately 12%.

T1 and T3 samples overlapped substantially, while patient-level heterogeneity
remained apparent.

This observation supported the decision to model patient identity explicitly.

## 5. Individual-gene differential expression

The paired DESeq2 analysis included 13,903 genes.

No genes reached:

```text
FDR < 0.10
```

or:

```text
FDR < 0.05
```

Thus, there is no evidence in this dataset for a broad set of individually
genome-wide-significant T3-vs-T1 genes under the chosen model.

This negative result is an important component of the analysis rather than a
reason to discard the dataset.

## 6. Pathway-level analysis

The complete ranked transcriptome was subsequently analyzed using GSEA.

Genes were ranked by their signed Wald statistics.

This revealed strong coordinated enrichment on the T3 side for processes
involving:

* mitotic chromosome segregation,
* kinetochore function,
* spindle checkpoint signaling,
* cell-cycle progression,
* DNA replication,
* DNA repair,
* ribosome-related processes.

The T1 side was characterized more strongly by:

* carbohydrate metabolism,
* lipid metabolism,
* extracellular-matrix organization,
* cell junctions,
* epithelial/tissue morphogenesis,
* integrin/focal-adhesion signaling.

## 7. GO results

GO Biological Process enrichment produced a large number of statistically
significant terms.

After redundancy reduction, the most prominent T3-associated themes included
spindle microtubule attachment to the kinetochore, nuclear division, chromosome
segregation, DNA replication and mitotic checkpoint signaling.

The T1-associated themes included cell-junction assembly, epithelial
morphogenesis, metabolic regulation and tissue-associated processes.

Because GO terms overlap extensively, the results were interpreted as groups of
related biological themes rather than as hundreds of independent findings.

## 8. KEGG results

KEGG analysis provided similar high-level structure.

T3-associated pathways included:

* Ribosome
* Cell cycle
* DNA replication
* Homologous recombination
* Oxidative phosphorylation
* antigen processing and presentation
* T-cell receptor signaling

T1-associated pathways included:

* other glycan degradation
* butanoate metabolism
* AMPK signaling
* integrin signaling
* focal adhesion
* cadherin signaling
* insulin signaling

The concordance between GO and KEGG provides internal consistency for the
pathway interpretation, while not constituting independent replication.

## 9. Representative candidate genes

Representative genes were chosen to connect pathway-level observations with
individual molecular features.

### T3-associated

**KNL1** — associated with kinetochore organization and chromosome segregation.

**RAD51** — associated with homologous recombination and DNA repair.

### T1-associated

**GYS1** — associated with glycogen/carbohydrate metabolism.

**INSR** — central to insulin-associated metabolic signaling.

Additional exploratory candidates included CDK1, CCNB1, SCD and TNC.

None of these candidates should be described as genome-wide significant DEGs.

## 10. Paired patient-level validation

Candidate expression changes were examined within each patient using VST
differences between T3 and T1.

Examples include:

| Gene  | Direction | Patients showing expected direction | Paired Wilcoxon p-value |
| ----- | --------- | ----------------------------------: | ----------------------: |
| KNL1  | T3 ↑      |                               17/24 |                  0.0135 |
| CDK1  | T3 ↑      |                               13/24 |                  0.1004 |
| CCNB1 | T3 ↑      |                               15/24 |                  0.1004 |
| RAD51 | T3 ↑      |                               16/24 |                  0.0891 |
| GYS1  | T1 ↑      |                               17/24 |                 0.00894 |
| SCD   | T1 ↑      |                               15/24 |                  0.1491 |
| INSR  | T1 ↑      |                               17/24 |                  0.0135 |
| TNC   | T1 ↑      |                               15/24 |                  0.1065 |

These results show that some pathway-derived candidates reproduce their
expected direction across a substantial fraction of patients.

Because candidate selection was informed by the same dataset, these tests
remain exploratory.

## 11. Integrated biological interpretation

The combined analysis suggests a transition from a T1 state enriched for
metabolic and tissue-interaction programs toward a T3 state with stronger
proliferative and genome-maintenance signatures.

The most coherent T3 signal involves the machinery required for chromosome
segregation, DNA replication, cell-cycle progression and repair.

The strongest T1 themes involve metabolic regulation together with extracellular
and adhesion-associated biology.

The transcriptomic evidence therefore appears to be distributed across
coordinated biological programs rather than concentrated in a small set of
dramatically changing individual genes.

## 12. Conclusion

This project demonstrates a complete paired transcriptomic analysis from
preprocessing through biological interpretation.

The principal statistical result is that no individual gene survives
genome-wide FDR correction.

The principal biological result is that pathway-level analysis nevertheless
reveals coordinated T3-associated proliferation/genome-maintenance programs and
T1-associated metabolic/tissue-interaction programs.

The candidate-gene analysis provides patient-level support for selected
directions but remains exploratory.

The resulting hypotheses should be tested in an independent cohort or through
experimental validation before being interpreted as definitive biomarkers or
mechanistic findings.
