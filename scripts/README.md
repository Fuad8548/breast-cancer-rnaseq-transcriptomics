# Project 1 Script Guide

## Purpose

These scripts preserve the computational workflow so the analysis can be reproduced after restarting R, reused in a future project, and audited by another reader.

## Files

### `01_extract_T1_T3.sh`
Creates the paired T1/T3 sample list and metadata, then extracts the 48 selected sample columns from the original 95-sample matrix.

Outputs:
- `r_objects/sample_ids.txt`
- `r_objects/metadata_T1_T3.csv`
- `r_objects/counts_T1_T3.tsv`

### `02_qc_filter_normalize.R`
Loads the selected count matrix and metadata, performs basic sample QC, filters low-expression genes, constructs the paired DESeq2 object, estimates size factors, creates normalized counts, performs VST, and generates PCA plots.

Outputs include:
- library-size QC
- detected-gene QC
- filtering summary
- normalization size factors
- normalized count matrix
- standard and paired PCA PDFs
- `r_objects/dds_normalized.rds`

### `03_deseq2.R`
Loads the prepared DESeq2 object, fits the negative-binomial model, extracts the explicit T3-vs-T1 contrast, and saves complete and sorted differential-expression results.

Outputs include:
- `results/differential_expression/results_names.txt`
- `results/differential_expression/DESeq2_T3_vs_T1_all_genes.csv`
- sorted result tables
- DESeq2 summary
- `r_objects/dds_fitted.rds`

### `run_pipeline.sh`
Runs the first three stages sequentially from the project root.

## Reproducibility principle

The R session is temporary. These scripts are the permanent record of the computational procedure. Intermediate `.rds` files are saved so later sessions can resume without rebuilding every object manually.

## Important

The scripts should be run from the project root, not from inside the `scripts/` directory.
