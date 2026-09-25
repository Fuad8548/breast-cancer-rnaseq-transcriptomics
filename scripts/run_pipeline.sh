#!/usr/bin/env bash
set -euo pipefail

# Run the reproducible Project 1 workflow from the project root.

bash scripts/01_extract_T1_T3.sh
Rscript scripts/02_qc_filter_normalize.R
Rscript scripts/03_deseq2.R

echo "Project 1 pipeline completed through differential expression."
