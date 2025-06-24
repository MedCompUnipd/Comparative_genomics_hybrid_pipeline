#!/bin/bash
set -euo pipefail # Exit immediately if a command exits with a non-zero status, if a variable is unset, or if a command in a pipeline fails.
shopt -s nullglob # Allows patterns that match no files to expand to a null string.

# === SAVE ORIGINAL DIRECTORY ===
ORIG_PWD="$(pwd)"
# Trap errors to return to original directory and exit cleanly on failure
trap 'cd "$ORIG_PWD" || exit; echo "[ERROR] Script interrupted or failed. Exiting."; exit 1' ERR

# --- Conda setup (not directly used by glsearch36, but good practice for full scripts) ---
# Function to check if a command exists
command_exists () {
    command -v "$1" &> /dev/null
}

if ! command_exists conda; then
  echo "[WARNING] 'conda' not found. Some dependent scripts might require it. Proceeding without Conda activation."
fi
# eval "$(conda shell.bash hook)" # Do not activate if this script doesn't explicitly need a conda env
# conda activate fasta3_env # This script's dependencies are expected to be in PATH

echo "--- Reciprocal Best Hit (RBH) Alignment Script ---"
echo "--------------------------------------------------"

# === FUNCTION: Write Temporary FASTA File ===
# Used to create single gene or contig FASTA files
write_temp_fasta() {
    local header="$1"
    local seq="$2"
    local file="$3"
    echo -e ">$header\n$seq" > "$file" || { echo "[ERROR] Failed to write temporary FASTA file: $file"; return 1; }
}

# === FUNCTION: Sanitize Filename ===
# Makes FASTA headers safe for use as filenames
sanitize_filename() {
    local name="$1"
    # Replaces problematic characters with underscores and removes non-alphanumeric except _. -
    echo "$name" | tr ' /|[]()=' '______' | tr -cd '[:alnum:]_.-'
}

# === USER INPUT ===
read -r -p "Enter path to reference genes FASTA file: " REFERENCE_GENES
read -r -p "Enter path to de novo genome FASTA file (your contig): " DENOVO_GENOME
read -r -p "Enter path to output directory for glsearch results: " OUTPUT_DIR

# === REQUIREMENTS CHECK ===
# Check if glsearch36 is installed (expected to be in PATH from setup script)
if ! command_exists glsearch36; then
    echo "[ERROR] 'glsearch36' not found. Please ensure the FASTA package is installed and its binaries are in your PATH."
    echo "Refer to the setup_environment.sh script for installation instructions."
    exit 1
fi

# Validate input files
if [[ ! -f "$REFERENCE_GENES" ]]; then
    echo "[ERROR] Reference genes FASTA file not found: $REFERENCE_GENES. Exiting."
    exit 1
fi
if [[ ! -f "$DENOVO_GENOME" ]]; then
    echo "[ERROR] De novo genome FASTA file not found: $DENOVO_GENOME. Exiting."
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR/contig_vs_genes" "$OUTPUT_DIR/genes_vs_contig" || { echo "[ERROR] Failed to create output directories."; exit 1; }

# --- FUNCTION: Align Contig vs. Genes ---
align_contig_to_genes() {
    local contig_file="$1"
    local genes_file="$2"
    local out_dir="$3"

    echo "[INFO] Running glsearch36: Contig ($contig_file) vs. each reference gene in ($genes_file)"
    local TEMP_GENE_FILE="${out_dir}/tmp_gene.fasta"
    local CURRENT_HEADER=""
    local CURRENT_SEQ=""
    local FIRST_GENE=true

    # Read genes file line by line
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" =~ ^> ]]; then
            if [[ -n "$CURRENT_HEADER" && -n "$CURRENT_SEQ" ]]; then
                local SAFE_HEADER=$(sanitize_filename "$CURRENT_HEADER")
                write_temp_fasta "$CURRENT_HEADER" "$CURRENT_SEQ" "$TEMP_GENE_FILE" || return 1
                echo "→ Aligning $SAFE_HEADER vs. Contig..."
                glsearch36 -m 0 "$contig_file" "$TEMP_GENE_FILE" > "$out_dir/${SAFE_HEADER}.glsearch.out" || {
                    echo "[ERROR] glsearch36 failed for $SAFE_HEADER vs contig."
                    rm -f "$TEMP_GENE_FILE"
                    return 1
                }
            fi
            CURRENT_HEADER=$(echo "$line" | sed 's/^>//')
            CURRENT_SEQ=""
        else
            CURRENT_SEQ+="$line"
        fi
    done < "$genes_file"

    # Process the last gene
    if [[ -n "$CURRENT_HEADER" && -n "$CURRENT_SEQ" ]]; then
        local SAFE_HEADER=$(sanitize_filename "$CURRENT_HEADER")
        write_temp_fasta "$CURRENT_HEADER" "$CURRENT_SEQ" "$TEMP_GENE_FILE" || return 1
        echo "→ Aligning $SAFE_HEADER vs. Contig (last gene)..."
        glsearch36 -m 0 "$contig_file" "$TEMP_GENE_FILE" > "$out_dir/${SAFE_HEADER}.glsearch.out" || {
            echo "[ERROR] glsearch36 failed for last gene ($SAFE_HEADER) vs contig."
            rm -f "$TEMP_GENE_FILE"
            return 1
        }
    fi
    rm -f "$TEMP_GENE_FILE" # Clean up temporary file
    echo "[INFO] Contig vs. genes alignment completed."
}

# --- FUNCTION: Align Genes vs. Contig ---
align_genes_to_contig() {
    local genes_file="$1"
    local contig_file="$2"
    local out_dir="$3"

    echo "[INFO] Running glsearch36: Each reference gene in ($genes_file) vs. Contig ($contig_file)"
    local TEMP_CONTIG_FILE="${out_dir}/tmp_contig.fasta"
    write_temp_fasta "denovo_contig" "$(cat "$contig_file" | grep -v '^>')" "$TEMP_CONTIG_FILE" || return 1

    local CURRENT_HEADER=""
    local CURRENT_SEQ=""

    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" =~ ^> ]]; then
            if [[ -n "$CURRENT_HEADER" && -n "$CURRENT_SEQ" ]]; then
                local SAFE_HEADER=$(sanitize_filename "$CURRENT_HEADER")
                local TEMP_GENE_FILE="${out_dir}/tmp_gene_for_${SAFE_HEADER}.fasta"
                write_temp_fasta "$CURRENT_HEADER" "$CURRENT_SEQ" "$TEMP_GENE_FILE" || {
                    echo "[ERROR] Failed to create temporary gene file for $CURRENT_HEADER"; rm -f "$TEMP_CONTIG_FILE"; return 1;
                }
                echo "→ Aligning Gene $SAFE_HEADER vs. Contig..."
                glsearch36 -m 0 "$TEMP_GENE_FILE" "$TEMP_CONTIG_FILE" > "$out_dir/${SAFE_HEADER}.glsearch.out" || {
                    echo "[ERROR] glsearch36 failed for gene ($SAFE_HEADER) vs contig."
                    rm -f "$TEMP_CONTIG_FILE" "$TEMP_GENE_FILE"
                    return 1
                }
                rm -f "$TEMP_GENE_FILE"
            fi
            CURRENT_HEADER=$(echo "$line" | sed 's/^>//')
            CURRENT_SEQ=""
        else
            CURRENT_SEQ+="$line"
        fi
    done < "$genes_file"

    # Process the last gene
    if [[ -n "$CURRENT_HEADER" && -n "$CURRENT_SEQ" ]]; then
        local SAFE_HEADER=$(sanitize_filename "$CURRENT_HEADER")
        local TEMP_GENE_FILE="${out_dir}/tmp_gene_for_${SAFE_HEADER}.fasta"
        write_temp_fasta "$CURRENT_HEADER" "$CURRENT_SEQ" "$TEMP_GENE_FILE" || {
            echo "[ERROR] Failed to create temporary gene file for the last gene $CURRENT_HEADER"; rm -f "$TEMP_CONTIG_FILE"; return 1;
        }
        echo "→ Aligning last Gene $SAFE_HEADER vs. Contig..."
        glsearch36 -m 0 "$TEMP_GENE_FILE" "$TEMP_CONTIG_FILE" > "$out_dir/${SAFE_HEADER}.glsearch.out" || {
            echo "[ERROR] glsearch36 failed for last gene ($SAFE_HEADER) vs contig."
            rm -f "$TEMP_CONTIG_FILE" "$TEMP_GENE_FILE"
            return 1
        }
        rm -f "$TEMP_GENE_FILE"
    fi

    rm -f "$TEMP_CONTIG_FILE" # Clean up temporary files
    echo "[INFO] Genes vs. contig alignment completed."
}


# === MAIN EXECUTION ===
echo "Starting FASTA alignment pipeline..."

echo "Executing alignment: contig vs each reference gene..."
align_contig_to_genes "$DENOVO_GENOME" "$REFERENCE_GENES" "$OUTPUT_DIR/contig_vs_genes" || { echo "[ERROR] Contig vs Genes alignment failed."; exit 1; }

echo "Executing alignment: each reference gene vs contig ..."
align_genes_to_contig "$REFERENCE_GENES" "$DENOVO_GENOME" "$OUTPUT_DIR/genes_vs_contig" || { echo "[ERROR] Genes vs Contig alignment failed."; exit 1; }

echo ""
echo "[COMPLETED] Alignments finished successfully!"
echo "Detailed results are in: $OUTPUT_DIR/contig_vs_genes/ and $OUTPUT_DIR/genes_vs_contig/"
