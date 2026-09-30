# Project 1 — Limitations

## 1. No genome-wide significant individual genes

The paired DESeq2 analysis identified no genes with adjusted p-value below
0.05, and none below 0.10.

Therefore, claims about specific genes must remain exploratory.

Moderate fold changes or nominal p-values cannot be interpreted as definitive
differential expression after multiple-testing correction.

## 2. Pathway enrichment does not prove mechanism

GSEA identifies coordinated statistical enrichment of genes belonging to
annotated biological sets.

It does not establish causal mechanisms.

For example, enrichment of cell-cycle or DNA-repair pathways indicates that genes
associated with those processes are shifted in the ranked expression profile;
it does not by itself prove increased cell division, DNA damage or a causal
biological process.

## 3. GO redundancy

The GO database contains many related and hierarchically nested terms.

Consequently, hundreds of significant terms may represent a smaller number of
underlying biological themes.

Redundancy reduction was therefore used before interpretation.

## 4. Disease-labelled pathway names require caution

Some KEGG pathways are named for diseases or clinical syndromes.

Enrichment of such a pathway does not mean that the study subjects have that
disease.

The appropriate interpretation is that genes represented in the corresponding
pathway contribute to the observed ranked expression signal.

## 5. Candidate selection is exploratory

Candidate genes were selected after examining pathway-level results.

Subsequent paired testing of these candidates is therefore not equivalent to a
fully independent validation experiment.

The candidate results should be interpreted as hypothesis-generating evidence.

## 6. Internal rather than external validation

GO and KEGG analyses use the same transcriptomic dataset.

Their convergence increases confidence that the observed signal is internally
consistent, but it does not constitute independent replication.

External validation in another patient cohort would provide substantially
stronger evidence.

## 7. Biological interpretation of broad annotations

Some annotations, particularly broad terms involving synapses, membrane
organization, adhesion or cytoskeletal processes, can arise from genes with
roles in multiple tissue types.

Such annotations should therefore be interpreted in the context of the actual
genes contributing to enrichment rather than taken literally from the pathway
name alone.

## 8. Transcriptomic association is not protein-level evidence

RNA abundance does not necessarily correspond directly to protein abundance,
protein activity or phenotype.

Experimental validation would therefore be necessary for mechanistic
conclusions.

## Overall interpretation

The strongest conclusion supported by this analysis is the presence of
coordinated pathway-level transcriptomic structure between T1 and T3, rather
than a set of individually validated differentially expressed genes.

The project should consequently be presented as a rigorous,
hypothesis-generating transcriptomic analysis rather than as a biomarker
validation study.
