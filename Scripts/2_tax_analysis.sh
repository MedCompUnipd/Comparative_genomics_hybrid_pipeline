#!/bin/bash
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

echo "--- Taxonomic Selection and QC Script ---"
echo "-----------------------------------------"

# === UTILITY FUNCTION ===
# Checks if an output file already exists for a given step.
# If it exists, it asks the user whether to skip the step or overwrite it.
should_run_step() {
  local output_check="$1"
  local step_name="$2"
  local step_dir
  step_dir=$(dirname "$output_check")

  if [[ -e "$output_check" ]]; then
    echo "[INFO] Output for '$step_name' already exists at: $output_check"
    read -r -p "Do you want to skip this step? (y to skip / n to overwrite): " choice
    if [[ "$choice" =~ ^[Yy]$ ]]; then
      echo "[INFO] Skipping $step_name..."
      return 1 # Skip the step
    else
      echo "[INFO] Overwriting $step_name: removing $step_dir..."
      rm -rf "$step_dir"
      mkdir -p "$step_dir"
      return 0
    fi
  else
    mkdir -p "$step_dir"
    return 0
  fi
}

# === USER INPUT ===
read -p "Directory containing FASTQ files (from Dorado/trimmed): " FASTQ_DIR
read -p "Kraken2 database path: " DB_PATH
read -p "Kraken2 output directory: " KRAKEN2_OUT_DIR
read -p "Filtered output directory: " FILTERED_OUT_DIR
read -p "TaxID to filter for (e.g., 10359 for HCMV): " TAXID
read -p "NanoPlot output directory: " NANOPLOT_DIR
read -p "Threads to use: " THREADS

MERGED_FASTQ="${FILTERED_OUT_DIR}/reads_taxid${TAXID}_merged.fastq"

# === SETUP ===
echo "[INFO] Creating output directories..."
mkdir -p "$KRAKEN2_OUT_DIR" "$FILTERED_OUT_DIR" "$NANOPLOT_DIR"

# === [1/4] Kraken2 Classification ===
if should_run_step "$KRAKEN2_OUT_DIR/all_reads_classified.out" "Kraken2 Classification"; then
    echo "[1/4] Running Kraken2 classification..."
    conda activate kraken_env || { echo "[ERROR] Failed to activate 'kraken_env'. Please check your Conda environments."; exit 1; }

    if [ ! -d "$DB_PATH" ]; then
        echo "[ERROR] Kraken2 database path not found: $DB_PATH. Exiting."
        conda deactivate
        exit 1
    fi

    shopt -s nullglob
    INPUT_FASTQ_FILES=("$FASTQ_DIR"/*.fastq)
    shopt -u nullglob

    if [ ${#INPUT_FASTQ_FILES[@]} -eq 0 ]; then
        echo "[ERROR] No FASTQ files found in $FASTQ_DIR. Aborting Kraken2."
        conda deactivate
        exit 1
    fi

    for fastq_file in "${INPUT_FASTQ_FILES[@]}"; do
        base_name=$(basename "$fastq_file" .fastq)
        k2_output="${KRAKEN2_OUT_DIR}/${base_name}_classified.out"
        k2_report="${KRAKEN2_OUT_DIR}/${base_name}_report.txt"
        echo "→ Processing $fastq_file with Kraken2..."
        kraken2 --db "$DB_PATH" --threads "$THREADS" --output "$k2_output" --report "$k2_report" "$fastq_file" || { echo "[ERROR] Kraken2 failed for $fastq_file."; conda deactivate; exit 1; }
    done
    echo "[INFO] Kraken2 classification completed. Reports in $KRAKEN2_OUT_DIR"
    conda deactivate
else
    echo "[1/4] Skipping Kraken2 Classification."
fi

# === [2/4] Filtering reads by TaxID and Merging ===
if should_run_step "$MERGED_FASTQ" "TaxID Filtering and Merging"; then
    echo "[2/4] Filtering reads for TaxID $TAXID and merging..."
    conda activate kraken_env || { echo "[ERROR] Failed to activate 'kraken_env' for filtering. Please check your Conda environments."; exit 1; }

    # Ensure previous output is cleared if overwriting
    > "$MERGED_FASTQ" # Clear the merged file if it exists, for fresh start

    shopt -s nullglob
    CLASSIFIED_FILES=("$KRAKEN2_OUT_DIR"/*_classified.out)
    shopt -u nullglob

    if [ ${#CLASSIFIED_FILES[@]} -eq 0 ]; then
        echo "[ERROR] No Kraken2 classified output files found in $KRAKEN2_OUT_DIR. Cannot filter."
        conda deactivate
        exit 1
    fi

    for CLASSIFIED_FILE in "${CLASSIFIED_FILES[@]}"; do
        BASE=$(basename "$CLASSIFIED_FILE" _classified.out)
        FASTQ_FILE_FOR_FILTERING="$FASTQ_DIR/${BASE}.fastq"
        IDS_FILE="$FILTERED_OUT_DIR/${BASE}_selected_ids.txt"
        OUT_FASTQ_FILE="$FILTERED_OUT_DIR/${BASE}_taxid${TAXID}.fastq"

        if [ ! -f "$FASTQ_FILE_FOR_FILTERING" ]; then
            echo "[WARNING] No matching FASTQ file found for $BASE ($FASTQ_FILE_FOR_FILTERING). Skipping."
            continue
        fi

        # Extract read IDs for the target TaxID
        awk -v taxid="$TAXID" '$1 == "C" && $3 == taxid {print $2}' "$CLASSIFIED_FILE" > "$IDS_FILE"

        if [ -s "$IDS_FILE" ]; then
            # Filter FASTQ based on extracted IDs
            seqtk subseq "$FASTQ_FILE_FOR_FILTERING" "$IDS_FILE" > "$OUT_FASTQ_FILE" || { echo "[ERROR] seqtk filtering failed for $FASTQ_FILE_FOR_FILTERING."; conda deactivate; exit 1; }
            cat "$OUT_FASTQ_FILE" >> "$MERGED_FASTQ"
            echo "[INFO] Filtered and added: $OUT_FASTQ_FILE"
        else
            echo "[INFO] No reads found for TaxID $TAXID in $BASE"
        fi
    done

    if [ ! -s "$MERGED_FASTQ" ]; then
        echo "[ERROR] No reads retained after filtering for TaxID $TAXID. Exiting."
        conda deactivate
        exit 1
    else
        echo "[INFO] Merged FASTQ file created: $MERGED_FASTQ"
    fi
    conda deactivate
else
    echo "[2/4] Skipping TaxID Filtering and Merging."
fi

# === [3/4] NanoPlot QC ===
if should_run_step "$NANOPLOT_DIR/NanoStats.txt" "NanoPlot QC"; then
    echo "[3/4] Running NanoPlot QC..."
    conda activate qc_env || { echo "[ERROR] Failed to activate 'qc_env'. Please check your Conda environments."; exit 1; }

    if ! command_exists NanoPlot; then
        echo "[ERROR] 'NanoPlot' is not installed. Aborting NanoPlot QC. Ensure it's in 'qc_env'."
        conda deactivate
        exit 1
    fi

    if [ ! -f "$MERGED_FASTQ" ]; then
        echo "[ERROR] Merged FASTQ file not found: $MERGED_FASTQ. Aborting NanoPlot QC."
        conda deactivate
        exit 1
    fi

    NanoPlot -t "$THREADS" --fastq "$MERGED_FASTQ" -o "$NANOPLOT_DIR" || { echo "[ERROR] NanoPlot failed."; conda deactivate; exit 1; }
    echo "[INFO] NanoPlot QC completed. Reports in $NANOPLOT_DIR"
    conda deactivate
else
    echo "[3/4] Skipping NanoPlot QC."
fi

# === [4/4] Qualimap QC (requires BAM, usually after alignment) ===
# This part of 2_tax_sel_qc.sh seems to imply post-assembly mapping.
# For standalone QC on filtered reads, a direct mapping is needed.
# Since the original script had NanoPlot as the primary QC for FASTQ,
# and Qualimap is more commonly used on BAM files from assembly/polishing,
# I will keep this step as a placeholder or remove it if not directly applicable here.
# Assuming this Qualimap step is intended for the *filtered reads* against something,
# or it's a remnant. The script is structured to run after NanoPlot.
# If you mean Qualimap for the filtered FASTQ files, it needs an alignment.
echo "[4/4] Running Qualimap QC (requires alignment if on FASTQ, otherwise for BAMs)..."
read -p "Do you want to run Qualimap QC on the filtered reads? (Requires an assembly to map to. y/n): " RUN_QUALIMAP_QC
if [[ "$RUN_QUALIMAP_QC" =~ ^[Yy]$ ]]; then
    read -p "Path to the reference FASTA for Qualimap mapping: " QUALIMAP_REF_FASTA
    if [ ! -f "$QUALIMAP_REF_FASTA" ]; then
        echo "[ERROR] Reference FASTA not found for Qualimap mapping: $QUALIMAP_REF_FASTA. Skipping Qualimap QC."
    else
        QUALIMAP_OUT_DIR="${NANOPLOT_DIR}/qualimap_out" # Subdirectory under NanoPlot for consistency
        mkdir -p "$QUALIMAP_OUT_DIR"
        echo "[INFO] Running Minimap2 mapping for Qualimap..."
        conda activate preprocessing_env || { echo "[ERROR] Failed to activate 'preprocessing_env' for minimap2. Please check your Conda environments."; exit 1; }
        MINIMAP_PAF="${QUALIMAP_OUT_DIR}/reads_vs_ref.paf"
        MINIMAP_SAM="${QUALIMAP_OUT_DIR}/reads_vs_ref.sam"
        MINIMAP_BAM="${QUALIMAP_OUT_DIR}/reads_vs_ref.bam"

        minimap2 -ax map-ont "$QUALIMAP_REF_FASTA" "$MERGED_FASTQ" > "$MINIMAP_PAF" 2> "${QUALIMAP_OUT_DIR}/minimap2.log" || { echo "[ERROR] Minimap2 failed."; conda deactivate; exit 1; }
        echo "[INFO] Converting PAF to SAM/BAM and sorting..."
        # Convert PAF to SAM, then to BAM, then sort and index
        # minimap2 -a output is already SAM format when using -a
        minimap2 -ax map-ont "$QUALIMAP_REF_FASTA" "$MERGED_FASTQ" | samtools view -bS - | samtools sort -o "$MINIMAP_BAM" - || { echo "[ERROR] Samtools conversion/sort failed."; conda deactivate; exit 1; }
        samtools index "$MINIMAP_BAM" || { echo "[ERROR] Samtools index failed."; conda deactivate; exit 1; }

        echo "[INFO] Running Qualimap multi-sample bamqc..."
        conda deactivate # Deactivate preprocessing_env
        conda activate qc_env || { echo "[ERROR] Failed to activate 'qc_env' for qualimap. Please check your Conda environments."; exit 1; }
        qualimap bamqc -bam "$MINIMAP_BAM" -outdir "$QUALIMAP_OUT_DIR/report" --java-mem-size=4G || { echo "[ERROR] Qualimap failed."; conda deactivate; exit 1; }
        echo "[INFO] Qualimap QC completed. Reports in $QUALIMAP_OUT_DIR/report"
        conda deactivate
    fi
else
    echo "[4/4] Skipping Qualimap QC."
fi

echo ""
echo "-----------------------------------------"
echo "Taxonomic selection and QC pipeline completed successfully!"
echo "Kraken2 outputs: $KRAKEN2_OUT_DIR"
echo "Filtered FASTQ (merged): $MERGED_FASTQ"
echo "NanoPlot reports: $NANOPLOT_DIR"
if [[ "$RUN_QUALIMAP_QC" =~ ^[Yy]$ ]]; then
    echo "Qualimap reports: ${NANOPLOT_DIR}/qualimap_out/report"
fi
