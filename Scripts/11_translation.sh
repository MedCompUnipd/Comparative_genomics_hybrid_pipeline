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
fasty36 -E 1e-5 -m 0 -s P250 "$QUERY_NUC" "$REF_PROTEIN" > "$OUTPUT_FORMAT0" || {
    echo "[ERROR] fasty36 failed for readable output. Check log for details."; exit 1;
}
echo "[INFO] Readable output saved to: $OUTPUT_FORMAT0"

echo -e "\n[2/3] Running fasty36 in tabular format (-m 8)..."
fasty36 -E 1e-5 -m 8 -s P250 "$QUERY_NUC" "$REF_PROTEIN" > "$OUTPUT_FORMAT6" || {
    echo "[ERROR] fasty36 failed for tabular output. Check log for details."; exit 1;
}
echo "[INFO] Tabular output saved to: $OUTPUT_FORMAT6"

echo -e "\n[3/3] Extracting best hits from tabular output and summarizing..."
if [[ ! -f "$OUTPUT_FORMAT6" || ! -s "$OUTPUT_FORMAT6" ]]; then
    echo "[ERROR] Tabular output file is missing or empty: $OUTPUT_FORMAT6. Cannot extract best hits."
    exit 1
fi

awk -F'\t' '
BEGIN {
    OFS = "\t";
    print "Query", "Target", "Identity", "Align_length", "Mismatches", "Gaps", \
          "Query_start", "Query_end", "Target_start", "Target_end", \
          "E-value", "Score", "Query_coverage(%)";
}

function abs(x) { return x < 0 ? -x : x }

{
    query = $1;
    target = $2;
    identity = $3;
    align_len = $4;
    mismatches = $5;
    gaps = $6;
    q_start = $7 + 0;
    q_end = $8 + 0;
    t_start = $9 + 0;
    t_end = $10 + 0;
    evalue = $11;
    score = $12;

    # Estrai prot_start e prot_end da "location_XXXX..YYYY"
    match(query, /location_([0-9]+)\.\.([0-9]+)/, loc);
    if (loc[1] && loc[2]) {
        prot_start = loc[1] + 0;
        prot_end = loc[2] + 0;
        prot_len = abs(prot_end - prot_start) + 1;

        covered_nt = abs(q_end - q_start) + 1;
        covered_aa = covered_nt / 3;
        query_coverage = (covered_aa / prot_len) * 100;

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%.1f\n", \
               query, target, identity, align_len, mismatches, gaps, \
               q_start, q_end, t_start, t_end, evalue, score, query_coverage;
    } else {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\tNA\n", \
               query, target, identity, align_len, mismatches, gaps, \
               q_start, q_end, t_start, t_end, evalue, score;
    }
}
' "$OUTPUT_FORMAT6" > "$OUTPUT_SUMMARY"


echo -e "\n----------------------------------------------------------"
echo "Nucleotide to Amino Acid alignment pipeline completed successfully!"
echo "Outputs:"
echo "  - Readable format: $OUTPUT_FORMAT0"
echo "  - Tabular format: $OUTPUT_FORMAT6"
echo "  - Best hits summary: $OUTPUT_SUMMARY"
