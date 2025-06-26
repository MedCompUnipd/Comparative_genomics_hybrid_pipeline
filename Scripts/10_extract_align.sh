#!/bin/bash
set -euo pipefail

# --- Conda setup ---
command_exists () {
    command -v "$1" &> /dev/null
}

if ! command_exists conda; then
  echo "[ERROR] Conda not found. Please ensure conda is installed and initialized."
  exit 1
fi
eval "$(conda shell.bash hook)"
conda activate fasta3_env || { echo "[ERROR] Failed to activate 'fasta3_env'. Please check your Conda environments."; exit 1; }

echo "--- Masked Region Extraction and Alignment Script ---"
echo "---------------------------------------------------"

# === REQUIREMENTS CHECK ===
if ! command_exists glsearch36; then
  echo "[ERROR] 'glsearch36' command not found. Please ensure the FASTA package is installed and its binaries are in your PATH."
  echo "Refer to the setup_environment.sh script for installation instructions."
  exit 1
fi

# === USER INPUT ===
read -p "Enter the path to the input masked FASTA file (e.g., from script 9): " INPUT_MASKED_FASTA
read -p "Enter the path to the output alignment directory: " ALIGNMENT_DIR
read -p "Enter the path to the output unaligned regions FASTA file (e.g., unaligned_regions.fasta): " OUTPUT_UNALIGNED_FASTA
read -p "Enter the path to the reference database FASTA file: " DB_FILE

# === INPUT VALIDATION ===
if [[ ! -f "$INPUT_MASKED_FASTA" ]]; then
  echo "[ERROR] Masked input FASTA file not found: $INPUT_MASKED_FASTA. Exiting."
  exit 1
fi

if [[ ! -f "$DB_FILE" ]]; then
  echo "[ERROR] Reference database FASTA file not found: $DB_FILE. Exiting."
  exit 1
fi

# === SETUP ===
mkdir -p "$ALIGNMENT_DIR" || { echo "[ERROR] Failed to create alignment directory: $ALIGNMENT_DIR"; exit 1; }
> "$OUTPUT_UNALIGNED_FASTA" || { echo "[ERROR] Failed to create/clear unaligned FASTA file: $OUTPUT_UNALIGNED_FASTA"; exit 1; }

TMP_FASTA_DIR="${ALIGNMENT_DIR}/tmp_region_fastas"
mkdir -p "$TMP_FASTA_DIR" || { echo "[ERROR] Failed to create temporary FASTA directory: $TMP_FASTA_DIR"; exit 1; }

# --- FUNCTION: Process a continuous region of non-N nucleotides ---
process_region() {
  local header="$1"
  local region_count="$2"
  local start="$3"
  local end="$4"
  local subseq="$5"

  local region_id="${header}_region_${region_count}_${start}_${end}"

  echo "[INFO] Processing region: $region_id"

  echo ">${region_id}" >> "$OUTPUT_UNALIGNED_FASTA"
  echo "$subseq" | fold -w 60 >> "$OUTPUT_UNALIGNED_FASTA"

  local region_file="${TMP_FASTA_DIR}/${region_id}.fasta"
  echo ">${region_id}" > "$region_file"
  echo "$subseq" | fold -w 60 >> "$region_file"

  local glsearch_output="${ALIGNMENT_DIR}/${region_id}.glsearch.out"
  echo "  → Running glsearch36 on $region_id against $DB_FILE. Output to $glsearch_output"
  glsearch36 -m 0 "$region_file" "$DB_FILE" > "$glsearch_output" || { echo "[WARNING] glsearch36 failed for region $region_id. Check $glsearch_output"; }
}

# === MAIN LOGIC ===
echo "[INFO] Reading input masked FASTA file: $INPUT_MASKED_FASTA"

CURRENT_HEADER=""
SEQUENCE=""
REGION_COUNT=1
IN_REGION=0
START_POS=0

while IFS= read -r line || [[ -n "$line" ]]; do
  line=$(echo "$line" | tr -d '\r')

  if [[ "$line" =~ ^">" ]]; then
    if [[ -n "$SEQUENCE" ]]; then
      SEQ_LEN=${#SEQUENCE}
      if [[ $IN_REGION -eq 1 ]]; then
        END_POS="$SEQ_LEN"
        SUBSEQ="${SEQUENCE:$((START_POS-1)):$(($END_POS - $START_POS + 1))}"
        process_region "$CURRENT_HEADER" "$REGION_COUNT" "$START_POS" "$END_POS" "$SUBSEQ"
        REGION_COUNT=$((REGION_COUNT + 1))
      fi
    fi
    CURRENT_HEADER=$(echo "$line" | cut -d' ' -f1 | sed 's/^>//' | sed 's/[^a-zA-Z0-9._-]/_/g')
    SEQUENCE=""
    REGION_COUNT=1
    IN_REGION=0
    START_POS=0
    echo "[INFO] Reading sequence for header: $CURRENT_HEADER"
  else
    SEQUENCE+="$line"
  fi
done < "$INPUT_MASKED_FASTA"

# Process last sequence
if [[ -n "$SEQUENCE" ]]; then
  SEQ_LEN=${#SEQUENCE}
  IN_REGION=0
  START_POS=0
  for (( i=0; i<SEQ_LEN; i++ )); do
    CHAR="${SEQUENCE:$i:1}"
    CURRENT_1_BASED_POS=$((i+1))

    if [[ "$CHAR" != "N" && $IN_REGION -eq 0 ]]; then
      START_POS="$CURRENT_1_BASED_POS"
      IN_REGION=1
    elif [[ ( "$CHAR" == "N" || "$CURRENT_1_BASED_POS" -eq "$SEQ_LEN" ) && $IN_REGION -eq 1 ]]; then
      if [[ "$CURRENT_1_BASED_POS" -eq "$SEQ_LEN" && "$CHAR" != "N" ]]; then
        END_POS="$CURRENT_1_BASED_POS"
      else
        END_POS=$((CURRENT_1_BASED_POS - 1))
      fi

      if [[ "$START_POS" -le "$END_POS" ]]; then
        SUBSEQ="${SEQUENCE:$((START_POS-1)):$(($END_POS - $START_POS + 1))}"
        process_region "$CURRENT_HEADER" "$REGION_COUNT" "$START_POS" "$END_POS" "$SUBSEQ"
        REGION_COUNT=$((REGION_COUNT + 1))
      fi
      IN_REGION=0
    fi
  done
fi

rm -rf "$TMP_FASTA_DIR" || echo "[WARNING] Failed to remove temporary directory: $TMP_FASTA_DIR"

echo ""
echo "---------------------------------------------------"
echo "Masked region extraction and alignment completed!"
echo "Unaligned regions are in: $OUTPUT_UNALIGNED_FASTA"
echo "Individual glsearch results are in: $ALIGNMENT_DIR"
