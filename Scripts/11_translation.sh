#!/usr/bin/env bash
set -euo pipefail

# --- Conda setup ---
# Function to check if a command exists
command_exists () {
    command -v "$1" &> /dev/null
}

if ! command_exists conda; then
  echo "[ERROR] Conda not found. Please ensure conda is installed and initialized."
  exit 1
fi
eval "$(conda shell.bash hook)"
conda activate fasta3_env || { echo "[ERROR] Failed to activate 'fasta3_env'. Please check your Conda environments."; exit 1; }

echo "--- Nucleotide to Amino Acid Alignment Script (fasty36) ---"
echo "----------------------------------------------------------"

# === REQUIREMENTS CHECK ===
# Check if fasty36 is installed
if ! command_exists fasty36; then
  echo "[ERROR] 'fasty36' command not found. Please ensure the FASTA package is installed and its binaries are in your PATH."
  echo "Refer to the setup_environment.sh script for installation instructions."
  exit 1
fi

# === USER INPUT ===
read -p "Enter path to query nucleotide file (e.g., denovo_genes_seq.fasta): " QUERY_NUC
read -p "Enter path to reference protein file (e.g., reference_proteins.fasta): " REF_PROTEIN
read -p "Enter path for readable output file (e.g., fasty_output_readable.txt): " OUTPUT_FORMAT0
read -p "Enter path for tabular output file (e.g., fasty_output_tabular.tsv): " OUTPUT_FORMAT6
read -p "Enter path for summary output file (e.g., fasty_best_hits_summary.tsv): " OUTPUT_SUMMARY

# === INPUT VALIDATION ===
if [[ ! -f "$QUERY_NUC" ]]; then
    echo "[ERROR] Query nucleotide file not found: $QUERY_NUC. Exiting."
    exit 1
fi
if [[ ! -f "$REF_PROTEIN" ]]; then
    echo "[ERROR] Reference protein file not found: $REF_PROTEIN. Exiting."
    exit 1
fi

# Create parent directories for output files
mkdir -p "$(dirname "$OUTPUT_FORMAT0")" "$(dirname "$OUTPUT_FORMAT6")" "$(dirname "$OUTPUT_SUMMARY")" || { echo "[ERROR] Failed to create output directories."; exit 1; }


# === MAIN EXECUTION ===

echo -e "\n[1/3] Running fasty36 in readable format (-m 0)..."
fasty36 -E 1e-5 -m 0 "$QUERY_NUC" "$REF_PROTEIN" > "$OUTPUT_FORMAT0" || {
    echo "[ERROR] fasty36 failed for readable output. Check log for details."; exit 1;
}
echo "[INFO] Readable output saved to: $OUTPUT_FORMAT0"

echo -e "\n[2/3] Running fasty36 in tabular format (-m 8)..."
fasty36 -E 1e-5 -m 8 "$QUERY_NUC" "$REF_PROTEIN" > "$OUTPUT_FORMAT6" || {
    echo "[ERROR] fasty36 failed for tabular output. Check log for details."; exit 1;
}
echo "[INFO] Tabular output saved to: $OUTPUT_FORMAT6"

echo -e "\n[3/3] Extracting best hits from tabular output and summarizing..."
if [[ ! -f "$OUTPUT_FORMAT6" || ! -s "$OUTPUT_FORMAT6" ]]; then
    echo "[ERROR] Tabular output file is missing or empty: $OUTPUT_FORMAT6. Cannot extract best hits."
    exit 1
fi

# AWK script to extract best hits
awk '
BEGIN {
    FS="\t"; # Tab separated file
    OFS="\t"; # Output tab separated
    print "Query", "Target", "Score", "Expect", "Identity", "Overlap", "Query_start", "Query_end", "Target_start", "Target_end", "Query_len", "Target_len"; # Header
}
{
    query = $1;
    target = $2;
    score = $3;
    expect = $4;
    identity = $5;
    overlap = $6;
    q_start = $7;
    q_end = $8;
    t_start = $9;
    t_end = $10;
    q_len = $11;
    t_len = $12;

    # Check if this query has been seen before or if this hit has a better score
    if (!(query in best_hits) || score > best_hits[query]["score"]) {
        best_hits[query]["target"] = target;
        best_hits[query]["score"] = score;
        best_hits[query]["expect"] = expect;
        best_hits[query]["identity"] = identity;
        best_hits[query]["overlap"] = overlap;
        best_hits[query]["q_start"] = q_start;
        best_hits[query]["q_end"] = q_end;
        best_hits[query]["t_start"] = t_start;
        best_hits[query]["t_end"] = t_end;
        best_hits[query]["q_len"] = q_len;
        best_hits[query]["t_len"] = t_len;
    }
}
END {
    # Print all best hits
    for (q in best_hits) {
        print q, best_hits[q]["target"], best_hits[q]["score"], best_hits[q]["expect"], \
              best_hits[q]["identity"], best_hits[q]["overlap"], \
              best_hits[q]["q_start"], best_hits[q]["q_end"], \
              best_hits[q]["t_start"], best_hits[q]["t_end"], \
              best_hits[q]["q_len"], best_hits[q]["t_len"];
    }
}' "$OUTPUT_FORMAT6" > "$OUTPUT_SUMMARY" || {
    echo "[ERROR] Failed to extract best hits using awk. Check $OUTPUT_FORMAT6"; exit 1;
}
echo "[INFO] Best hits summary saved to: $OUTPUT_SUMMARY"


echo -e "\n----------------------------------------------------------"
echo "Nucleotide to Amino Acid alignment pipeline completed successfully!"
echo "Outputs:"
echo "  - Readable format: $OUTPUT_FORMAT0"
echo "  - Tabular format: $OUTPUT_FORMAT6"
echo "  - Best hits summary: $OUTPUT_SUMMARY"
