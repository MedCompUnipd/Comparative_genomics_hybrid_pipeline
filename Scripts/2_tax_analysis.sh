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
# Checks if an output file or directory already exists for a given step.
# If it exists, it asks the user whether to skip the step or overwrite it.
# If overwriting, it removes the relevant directory for a clean start.
should_run_step() {
    local output_check="$1" # Path to a key output file/directory to check for existence
    local step_name="$2"    # Name of the step for user feedback
    local step_dir=""       # Directory associated with the step's output

    # Determine the directory to check/clear. If output_check is a directory, use it.
    # Otherwise, use its parent directory.
    if [[ -d "$output_check" ]]; then
        step_dir="$output_check"
    else
        step_dir=$(dirname "$output_check")
    fi

    if [[ -e "$output_check" || -d "$output_check" ]]; then
        echo "[INFO] Output for '$step_name' already exists at: $output_check"
        read -r -p "Do you want to skip this step? (y to skip / n to overwrite): " choice
        if [[ "$choice" =~ ^[Yy]$ ]]; then
            echo "[INFO] Skipping $step_name..."
            return 1 # Skip the step
        else
            # Ensure step_dir is not empty to prevent 'rm -rf /'
            if [[ -n "$step_dir" && "$step_dir" != "." && "$step_dir" != "/" ]]; then
                echo "[INFO] Overwriting $step_name: removing $step_dir..."
                rm -rf "$step_dir"
            fi
            # Recreate the base directory if it was removed and is needed by the step
            mkdir -p "$step_dir"
            return 0 # Overwrite, proceed with the step
        fi
    else
        # If output doesn't exist, ensure the output directory for the step is created
        mkdir -p "$step_dir"
        return 0 # Run the step
    fi
}


# === USER INPUT ===
read -p "Directory containing FASTQ files (from Dorado/trimmed): " FASTQ_DIR
read -p "Kraken2 database path (e.g., /path/to/kraken2_db): " DB_PATH
read -p "Kraken2 output directory: " KRAKEN2_OUT_DIR
read -p "Filtered output directory: " FILTERED_OUT_DIR
read -p "TaxID to filter for (e.g., 10359 for HCMV): " TAXID
read -p "NanoPlot output directory: " NANOPLOT_DIR
read -p "Threads to use: " THREADS

MERGED_FASTQ="${FILTERED_OUT_DIR}/reads_taxid${TAXID}_merged.fastq"

# === SETUP ===
echo "[INFO] Creating output directories..."
# Directories are also created by should_run_step, but this ensures initial existence
mkdir -p "$KRAKEN2_OUT_DIR" "$FILTERED_OUT_DIR" "$NANOPLOT_DIR"

# === [1/4] Kraken2 Database Build ===
# Check if the Kraken2 database needs to be built or if user wants to overwrite
if should_run_step "${DB_PATH}/hash.k2d" "Kraken2 Database Build"; then
    echo "[1/4] Setting up Kraken2 database..."
    conda activate kraken2_env || { echo "[ERROR] Failed to activate 'kraken2_env'. Please check your Conda environments."; exit 1; }

    # Check if DB_PATH is provided
    if [ -z "$DB_PATH" ]; then
        echo "[ERROR] Kraken2 database path is empty. Please provide a valid path."
        conda deactivate
        exit 1
    fi

    # Ensure the DB_PATH directory exists before building
    mkdir -p "$DB_PATH"

    echo "[INFO] Attempting to download and build Kraken2 standard database in $DB_PATH..."
    # kraken2-build --standard automatically downloads and builds the database
    kraken2-build --standard --db "$DB_PATH" \
        || { echo "[ERROR] Failed to download or build Kraken2 standard database at $DB_PATH."; conda deactivate; exit 1; }

    if [[ ! -f "$DB_PATH/hash.k2d" ]]; then
        echo "[ERROR] Kraken2 database was not properly built (missing hash.k2d file in $DB_PATH). Exiting."
        conda deactivate
        exit 1
    fi
    echo "[INFO] Kraken2 database built successfully."
    conda deactivate
else
    echo "[1/4] Skipping Kraken2 Database Build."
fi


# === [2/4] Kraken2 Classification ===
# Check if Kraken2 classification output exists and if user wants to skip
if should_run_step "$KRAKEN2_OUT_DIR" "Kraken2 Classification"; then
    echo "[2/4] Running Kraken2 classification..."
    conda activate kraken2_env || { echo "[ERROR] Failed to activate 'kraken2_env'. Please check your Conda environments."; exit 1; }

    # Verify Kraken2 database path exists and contains necessary files
    if [ ! -d "$DB_PATH" ] || [ ! -f "$DB_PATH/hash.k2d" ] || [ ! -f "$DB_PATH/taxo.k2d" ]; then
        echo "[ERROR] Kraken2 database is incomplete or not found at $DB_PATH. Please ensure it is built."
        conda deactivate
        exit 1
    fi

    shopt -s nullglob # Enable nullglob to prevent pattern from expanding to itself if no matches
    INPUT_FASTQ_FILES=("$FASTQ_DIR"/*.fastq "$FASTQ_DIR"/*.fastq.gz) # Now includes .fastq.gz files
    shopt -u nullglob # Disable nullglob

    if [ ${#INPUT_FASTQ_FILES[@]} -eq 0 ]; then
        echo "[ERROR] No FASTQ or FASTQ.GZ files found in $FASTQ_DIR. Aborting Kraken2 classification."
        conda deactivate
        exit 1
    fi

    for fastq_file in "${INPUT_FASTQ_FILES[@]}"; do
        if [ ! -s "$fastq_file" ]; then # Check if file is empty
            echo "[WARNING] Input FASTQ file is empty or missing: $fastq_file. Skipping this file."
            continue
        fi

        # Check if the file might be gzipped but misnamed (e.g., .fastq instead of .fastq.gz)
        # First two bytes of a gzip file are 1f 8b (hex)
        first_bytes=$(head -c 2 "$fastq_file" | xxd -p)
        if [[ "$first_bytes" == "1f8b" ]] && [[ "$fastq_file" != *.gz ]]; then
            echo "[WARNING] File '$fastq_file' appears to be gzipped but lacks a .gz extension. Kraken2 usually handles this, but if errors persist, consider renaming it to .fastq.gz or decompressing it."
        fi

        # Robust base_name extraction to handle both .fastq and .fastq.gz extensions
        file_name_only=$(basename "$fastq_file") # Removed 'local'
        base_name=""                             # Removed 'local'
        if [[ "$file_name_only" == *.fastq.gz ]]; then
            base_name="${file_name_only%.fastq.gz}"
        elif [[ "$file_name_only" == *.fastq ]]; then
            base_name="${file_name_only%.fastq}"
        else
            # Fallback for unexpected extensions, though globbing should reduce this.
            # In production, you might want to error out here or skip.
            echo "[WARNING] Unexpected file extension for $file_name_only. Using full name as base."
            base_name="$file_name_only"
        fi

        k2_output="${KRAKEN2_OUT_DIR}/${base_name}_classified.out"
        k2_report="${KRAKEN2_OUT_DIR}/${base_name}_report.txt"
        echo "? Processing $fastq_file with Kraken2..."
        kraken2 --db "$DB_PATH" --output "$k2_output" --report "$k2_report" "$fastq_file" \
            || { echo "[ERROR] Kraken2 failed for $fastq_file. Please check its format and integrity. Exiting."; conda deactivate; exit 1; }
    done
    echo "[INFO] Kraken2 classification completed. Reports in $KRAKEN2_OUT_DIR"
    conda deactivate
else
    echo "[2/4] Skipping Kraken2 Classification."
fi

# === [3/4] Filtering reads by TaxID and Merging ===
# Check if the merged FASTQ file exists and if user wants to skip
if should_run_step "$MERGED_FASTQ" "TaxID Filtering and Merging"; then
    echo "[3/4] Filtering reads for TaxID $TAXID and merging..."
    conda activate kraken2_env || { echo "[ERROR] Failed to activate 'kraken2_env' for filtering. Please check your Conda environments."; exit 1; }

    # Ensure seqtk is available in the environment
    if ! command_exists seqtk; then
        echo "[ERROR] 'seqtk' is not installed or not in PATH. Aborting filtering. Ensure it's in 'kraken2_env'."
        conda deactivate
        exit 1
    fi

    # Clear the merged file if it exists, for fresh start (should be handled by should_run_step via rm -rf on its dir)
    # Re-touch the file to ensure it exists for appending
    touch "$MERGED_FASTQ"

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
        # Use a more robust approach to find the original FASTQ file, considering both .fastq and .fastq.gz
        FASTQ_FILE_FOR_FILTERING=""
        if [ -f "$FASTQ_DIR/${BASE}.fastq" ]; then
            FASTQ_FILE_FOR_FILTERING="$FASTQ_DIR/${BASE}.fastq"
        elif [ -f "$FASTQ_DIR/${BASE}.fastq.gz" ]; then
            FASTQ_FILE_FOR_FILTERING="$FASTQ_DIR/${BASE}.fastq.gz"
        fi

        IDS_FILE="$FILTERED_OUT_DIR/${BASE}_selected_ids.txt"
        OUT_FASTQ_FILE="$FILTERED_OUT_DIR/${BASE}_taxid${TAXID}.fastq"

        if [ -z "$FASTQ_FILE_FOR_FILTERING" ]; then
            echo "[WARNING] No matching original FASTQ or FASTQ.GZ file found for base '$BASE' in '$FASTQ_DIR'. Skipping filtering for this classified file."
            continue
        fi

        # Extract read IDs for the target TaxID
        awk -v taxid="$TAXID" '$1 == "C" && $3 == taxid {print $2}' "$CLASSIFIED_FILE" > "$IDS_FILE"

        if [ -s "$IDS_FILE" ]; then
            # Filter FASTQ based on extracted IDs
            seqtk subseq "$FASTQ_FILE_FOR_FILTERING" "$IDS_FILE" > "$OUT_FASTQ_FILE" \
                || { echo "[ERROR] seqtk filtering failed for $FASTQ_FILE_FOR_FILTERING."; conda deactivate; exit 1; }
            cat "$OUT_FASTQ_FILE" >> "$MERGED_FASTQ"
            echo "[INFO] Filtered and added: $OUT_FASTQ_FILE"
        else
            echo "[INFO] No reads found for TaxID $TAXID in '$BASE'"
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
    echo "[3/4] Skipping TaxID Filtering and Merging."
fi

# === [4/4] NanoPlot QC ===
# Check if NanoPlot output exists and if user wants to skip
if should_run_step "$NANOPLOT_DIR" "NanoPlot QC"; then
    echo "[4/4] Running NanoPlot QC..."
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

    NanoPlot \
        -t "$THREADS" \
        --fastq "$MERGED_FASTQ" \
        --loglength \
        --plots dot kde \
        --title "QC Report - Merged TaxID ${TAXID}" \
        --prefix "taxid${TAXID}_report" \
        -o "$NANOPLOT_DIR" \
        || { echo "[ERROR] NanoPlot failed."; conda deactivate; exit 1; }
    echo "[INFO] NanoPlot QC completed. Reports in $NANOPLOT_DIR"
    conda deactivate
else
    echo "[4/4] Skipping NanoPlot QC."
fi

echo ""
echo "-----------------------------------------"
echo "Taxonomic selection and QC pipeline completed successfully!"
echo "Kraken2 DB path: $DB_PATH"
echo "Kraken2 outputs: $KRAKEN2_OUT_DIR"
echo "Filtered FASTQ (merged): $MERGED_FASTQ"
echo "NanoPlot reports: $NANOPLOT_DIR"
