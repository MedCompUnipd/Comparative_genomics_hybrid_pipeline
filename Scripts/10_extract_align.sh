#!/bin/bash
set -euo pipefail

echo ""
echo "=== Masked Region Extraction & glsearch36 Alignment ==="

# --- Check conda ---
command -v conda >/dev/null || { echo "[ERROR] conda not found in PATH"; exit 1; }
eval "$(conda shell.bash hook)"
conda activate fasta3_env || { echo "[ERROR] Failed to activate fasta3_env"; exit 1; }

# --- Check glsearch36 ---
command -v glsearch36 >/dev/null || { echo "[ERROR] glsearch36 not found in PATH"; exit 1; }

# === User Input ===
read -p "Masked FASTA input path: " INPUT_MASKED_FASTA
read -p "Output alignment directory: " ALIGNMENT_DIR
read -p "Output unaligned regions FASTA: " OUTPUT_UNALIGNED_FASTA
read -p "Reference DB FASTA: " DB_FILE

[[ -f "$INPUT_MASKED_FASTA" ]] || { echo "[ERROR] Input FASTA not found"; exit 1; }
[[ -f "$DB_FILE" ]] || { echo "[ERROR] DB FASTA not found"; exit 1; }

# --- Setup ---
mkdir -p "$ALIGNMENT_DIR"
TMP_FASTA_DIR="${ALIGNMENT_DIR}/tmp"
mkdir -p "$TMP_FASTA_DIR"
> "$OUTPUT_UNALIGNED_FASTA"

echo ""
echo "[INFO] Processing $INPUT_MASKED_FASTA"

CURRENT_HEADER=""
SEQUENCE=""
REGION_COUNT=1

clean_header() {
  echo "$1" | cut -d' ' -f1 | sed 's/^>//; s/[^A-Za-z0-9_.-]/_/g'
}

write_region() {
  local header="$1"
  local count="$2"
  local start="$3"
  local end="$4"
  local subseq="$5"

  local region_id="${header}_region_${count}_${start}_${end}"
  echo "[INFO] Found region: $region_id"

  if (( ${#subseq} < 20 )); then
    echo "[INFO] Skipping too-short region ($region_id length=${#subseq})"
    return
  fi

  if ! grep -Eq '^[ACGTacgt]+$' <<< "$subseq"; then
    echo "[WARNING] Region $region_id has invalid bases; skipping"
    return
  fi

  {
    echo ">${region_id}"
    echo "$subseq" | fold -w 60
  } >> "$OUTPUT_UNALIGNED_FASTA"

  local region_fasta="${TMP_FASTA_DIR}/${region_id}.fasta"
  {
    echo ">${region_id}"
    echo "$subseq" | fold -w 60
  } > "$region_fasta"

  if [[ ! -s "$region_fasta" ]]; then
    echo "[WARNING] Region FASTA not created or empty: $region_fasta"
    return
  fi

  ls -lh "$region_fasta"

  echo "    [RUN] glsearch36 - $region_id"
  local gl_out="${ALIGNMENT_DIR}/${region_id}.glsearch.out"
  if ! glsearch36 -m 0 "$region_fasta" "$DB_FILE" > "$gl_out"; then
    echo "[WARNING] glsearch36 failed for $region_id"
  fi
}

# === MAIN LOOP ===
while IFS= read -r line || [[ -n "$line" ]]; do
  line=$(echo "$line" | tr -d '\r')
  if [[ "$line" == ">"* ]]; then
    if [[ -n "$SEQUENCE" ]]; then
      SEQ_LEN=${#SEQUENCE}
      IN_REGION=0
      for (( i=0; i<SEQ_LEN; i++ )); do
        CHAR="${SEQUENCE:i:1}"
        POS=$((i+1))
        if [[ "$CHAR" != "N" && $IN_REGION -eq 0 ]]; then
          START_POS=$POS
          IN_REGION=1
        elif [[ ( "$CHAR" == "N" || "$POS" -eq "$SEQ_LEN" ) && $IN_REGION -eq 1 ]]; then
          if [[ "$CHAR" == "N" ]]; then
            END_POS=$((POS - 1))
          else
            END_POS=$POS
          fi
          if (( START_POS <= END_POS )); then
            SUBSEQ="${SEQUENCE:START_POS-1:END_POS - START_POS + 1}"
            write_region "$CURRENT_HEADER" "$REGION_COUNT" "$START_POS" "$END_POS" "$SUBSEQ"
            REGION_COUNT=$((REGION_COUNT + 1))
          fi
          IN_REGION=0
        fi
      done
    fi
    CURRENT_HEADER=$(clean_header "$line")
    SEQUENCE=""
    REGION_COUNT=1
  else
    SEQUENCE+="$line"
  fi
done < "$INPUT_MASKED_FASTA"

# Last record
if [[ -n "$SEQUENCE" ]]; then
  SEQ_LEN=${#SEQUENCE}
  IN_REGION=0
  for (( i=0; i<SEQ_LEN; i++ )); do
    CHAR="${SEQUENCE:i:1}"
    POS=$((i+1))
    if [[ "$CHAR" != "N" && $IN_REGION -eq 0 ]]; then
      START_POS=$POS
      IN_REGION=1
    elif [[ ( "$CHAR" == "N" || "$POS" -eq "$SEQ_LEN" ) && $IN_REGION -eq 1 ]]; then
      if [[ "$CHAR" == "N" ]]; then
        END_POS=$((POS - 1))
      else
        END_POS=$POS
      fi
      if (( START_POS <= END_POS )); then
        SUBSEQ="${SEQUENCE:START_POS-1:END_POS - START_POS + 1}"
        write_region "$CURRENT_HEADER" "$REGION_COUNT" "$START_POS" "$END_POS" "$SUBSEQ"
        REGION_COUNT=$((REGION_COUNT + 1))
      fi
      IN_REGION=0
    fi
  done
fi

echo ""
echo "[INFO] Extraction and glsearch36 done"

# === SUMMARY PARSING ===
SUMMARY_TSV="${ALIGNMENT_DIR}/glsearch_summary.tsv"
echo -e "Region_ID\tContig_Start\tContig_End\tGene_ID\tGene_Start\tGene_End\tIdentity" > "$SUMMARY_TSV"

for gl_out in "$ALIGNMENT_DIR"/*.glsearch.out; do
  region_id=$(basename "$gl_out" .glsearch.out)

  contig_start=$(echo "$region_id" | awk -F'_' '{print $(NF-1)}')
  contig_end=$(echo "$region_id" | awk -F'_' '{print $NF}')

  best_line=$(awk '/^The best scores are:/ {getline; print}' "$gl_out")
  if [[ -n "$best_line" && "$best_line" =~ ^[^[:space:]]+ ]]; then
    gene_id=$(echo "$best_line" | awk '{print $1}')
    if [[ "$gene_id" =~ ([^:]+):([0-9]+)-([0-9]+) ]]; then
      gene_start="${BASH_REMATCH[2]}"
      gene_end="${BASH_REMATCH[3]}"
    else
      gene_start="NA"
      gene_end="NA"
    fi

    identity=$(awk '/global\/local score:/ {
      for (i=1; i<=NF; i++) {
        if ($i ~ /%/) {
          gsub("%","",$i); print $i; exit
        }
      }
    }' "$gl_out")
  else
    gene_id="NA"; gene_start="NA"; gene_end="NA"; identity="NA"
  fi

  echo -e "${region_id}\t${contig_start}\t${contig_end}\t${gene_id}\t${gene_start}\t${gene_end}\t${identity}" >> "$SUMMARY_TSV"
done

echo ""
echo "[INFO] Summary written to $SUMMARY_TSV"

# === Cleanup ===
rm -rf "$TMP_FASTA_DIR"

echo ""
echo "=== DONE ==="
echo "  Unaligned regions: $OUTPUT_UNALIGNED_FASTA"
echo "  Alignment outputs: $ALIGNMENT_DIR"
echo "  Summary: $SUMMARY_TSV"
