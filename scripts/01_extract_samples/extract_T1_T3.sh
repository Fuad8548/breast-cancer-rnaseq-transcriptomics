#!/usr/bin/env bash
set -euo pipefail

# Project 1: Breast-cancer RNA-seq transcriptomics
# Stage 01: Extract the paired T1/T3 samples from the original GEO matrix.

RAW="data/raw/GSE122630_Processed_data.txt.gz"
PROCESSED="data/processed"

mkdir -p "$PROCESSED"

# Extract sample IDs for the T1/T3 subset.
gzip -dc "$RAW" \
  | head -1 \
  | tr '\t' '\n' \
  | tail -n +2 \
  | grep -E '_(T1|T3)$' \
  > "$PROCESSED/sample_ids.txt"

# Build metadata from the sample names.
{
    echo "sample_id,patient_id,timepoint"
    while IFS= read -r sample; do
        patient="${sample#Patient_}"
        patient="${patient%_*}"
        timepoint="${sample##*_}"
        printf '%s,%s,%s\n' "$sample" "$patient" "$timepoint"
    done < "$PROCESSED/sample_ids.txt"
} > "$PROCESSED/metadata_T1_T3.csv"

# Extract only the selected sample columns, keeping ID_REF as column 1.
awk -F'\t' -v OFS='\t' '
NR==FNR {
    wanted[$1]=1
    next
}
FNR==1 {
    for (i=1; i<=NF; i++) {
        if (i==1 || ($i in wanted)) idx[++n]=i
    }
    for (j=1; j<=n; j++) {
        if (j>1) printf OFS
        printf "%s", $idx[j]
    }
    print ""
    next
}
{
    for (j=1; j<=n; j++) {
        if (j>1) printf OFS
        printf "%s", $idx[j]
    }
    print ""
}' "$PROCESSED/sample_ids.txt" \
   <(gzip -dc "$RAW") \
   > "$PROCESSED/counts_T1_T3.tsv"

# Basic validation.
rows=$(wc -l < "$PROCESSED/counts_T1_T3.tsv")
cols=$(head -1 "$PROCESSED/counts_T1_T3.tsv" | awk -F'\t' '{print NF}')
meta_rows=$(($(wc -l < "$PROCESSED/metadata_T1_T3.csv") - 1))

printf 'Count matrix: %s rows, %s columns\n' "$rows" "$cols"
printf 'Metadata samples: %s\n' "$meta_rows"

if [[ "$cols" -ne 49 ]]; then
    echo "ERROR: expected 49 columns (1 gene ID + 48 samples)." >&2
    exit 1
fi

if [[ "$meta_rows" -ne 48 ]]; then
    echo "ERROR: expected 48 selected samples." >&2
    exit 1
fi

echo "Stage 01 completed successfully."
