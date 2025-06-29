#!/bin/bash
set -euo pipefail # Exit immediately if a command exits with a non-zero status, if a variable is unset, or if a command in a pipeline fails.

# --- Conda Setup ---
# Function to check if a command exists
command_exists () {
    command -v "$1" &> /dev/null
}

if ! command_exists conda; then
  echo "[ERROR] Conda not found. Please ensure conda is installed and initialized."
  exit 1
fi
eval "$(conda shell.bash hook)" # Initialize Conda for the current shell session

# --- User Input for Logfile ---
read -p "Enter path to save the LOGfile (e.g., data/project_name/pipeline_ont.log): " LOGFILE_PATH

# Ensure the log file directory exists
LOGFILE_DIR=$(dirname "$LOGFILE_PATH")
mkdir -p "$LOGFILE_DIR" || { echo "[ERROR] Failed to create log directory: $LOGFILE_DIR. Exiting."; exit 1; }

echo "=== RBH & Masking & Alignment Pipeline Started ===" | tee "$LOGFILE_PATH"
echo "" | tee -a "$LOGFILE_PATH"

# --- Preliminary Suggestions ---
echo "Preliminary suggestions:" | tee -a "$LOGFILE_PATH"
echo " 1) Download reference gene FASTA files (nucleotide and protein) and GFF3 annotations from NCBI." | tee -a "$LOGFILE_PATH"
echo "    For example:" | tee -a "$LOGFILE_PATH"
echo "      - Visit https://www.ncbi.nlm.nih.gov/genome/ to download genomes of reference strains." | tee -a "$LOGFILE_PATH"
echo "      - Use wget for FTP downloads, e.g.:" | tee -a "$LOGFILE_PATH"
echo "        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/XXXX/YYYYY/GCF_XXXXX_genomic.fna.gz" | tee -a "$LOGFILE_PATH"
echo "      - Use Entrez Direct for batch downloads:" | tee -a "$LOGFILE_PATH"
echo "        https://www.ncbi.nlm.nih.gov/books/NBK179288/" | tee -a "$LOGFILE_PATH"
echo "" | tee -a "$LOGFILE_PATH"
echo " 2) Make sure all your individual pipeline scripts (1_*.sh, 7_*.py etc.) are in the same directory as this script" | tee -a "$LOGFILE_PATH"
echo "    or update the paths in this script accordingly." | tee -a "$LOGFILE_PATH"
echo "" | tee -a "$LOGFILE_PATH"

read -p "Press ENTER to start the pipeline..."

# --- Function to run a pipeline step ---
run_step() {
    local STEP_NAME="$1"
    local COMMAND="$2"
    local CONDA_ENV="$3" # Optional: Conda environment to activate for this step

    echo "" | tee -a "$LOGFILE_PATH"
    echo "--- $STEP_NAME ---" | tee -a "$LOGFILE_PATH"
    echo "" | tee -a "$LOGFILE_PATH"

    read -p "Do you want to run this step? (y/n): " answer
    if [[ "$answer" == "y" || "$answer" == "Y" ]]; then
        if [[ -n "$CONDA_ENV" ]]; then
            echo "[INFO] Activating Conda environment: $CONDA_ENV" | tee -a "$LOGFILE_PATH"
            # Use `|| { ...; exit 1; }` for robust error handling during activation
            conda activate "$CONDA_ENV" 2>&1 | tee -a "$LOGFILE_PATH" || { echo "[ERROR] Failed to activate '$CONDA_ENV'. Exiting."; exit 1; }
        fi

        echo "[INFO] Executing command: $COMMAND" | tee -a "$LOGFILE_PATH"
        eval "$COMMAND" 2>&1 | tee -a "$LOGFILE_PATH" || {
            echo "[ERROR] Step '$STEP_NAME' failed. Check the log file: $LOGFILE_PATH" | tee -a "$LOGFILE_PATH"
            exit 1
        }
        echo "[INFO] Step '$STEP_NAME' completed successfully." | tee -a "$LOGFILE_PATH"
        
        if [[ -n "$CONDA_ENV" ]]; then
            echo "[INFO] Deactivating Conda environment: $CONDA_ENV" | tee -a "$LOGFILE_PATH"
            conda deactivate 2>&1 | tee -a "$LOGFILE_PATH"
        fi
    else
        echo "Skipped: $STEP_NAME" | tee -a "$LOGFILE_PATH"
    fi
}

# --- Pipeline Steps ---
# Each call to run_step now includes the optional Conda environment to activate.
# The environment will be activated before the command and deactivated afterwards.

run_step "Step 1: Dorado basecalling and trimming adapters" "bash 1_basecalling.sh" "preprocessing_env"
run_step "Step 2: Taxonomic analysis and reads selection" "bash 2_tax_analysis.sh" "kraken2_env" # kraken_env needed for initial part, NanoPlot/Qualimap use qc_env which is activated internally by 2_tax_sel_qc.sh
run_step "Step 3: Assembly and polishing" "bash 3_asm_pol.sh" "" # Or leave empty if handled internally, but 3_asm_pol.sh itself activates many
run_step "Step 4: Using Illumina reads for hybrid assembly and polishing" "bash 4_hybrid.sh" "" # 4_hybrid.sh handles its own env activations
run_step "Step 5: Comparing de novo assembly to reference and performing variant calling" "bash 5_variant_calling.sh" "qc_env" # or leave empty if handled internally
run_step "Step 6: Performing Reciprocal Best Hit (RBH) alignment" "bash 6_rbh.sh" "" # glsearch36 is expected in PATH, no specific conda env needed by this script
run_step "Step 7: Parse RBH glsearch output and create lists per identity" "python3 7_rbh_lists.py" "base" # Assuming pandas/Biopython in base, or specify another
run_step "Step 8: Extracting sequences based on TSV coordinates" "python3 8_extract_genes.py" "fasta3_env"
run_step "Step 9: Masking sequences with identity thresholds" "python3 9_masking.py" "" # No specific env required for this Python script beyond Python itself
run_step "Step 10: Extracting unaligned regions and aligning to whole NCBI dataset" "bash 10_extract_align.sh" "fasta3_env"
run_step "Step 11: Translating nucleotidic sequences to aminoacidic sequences" "bash 11_translation.sh" "fasta3_env"
run_step "Step 12: Exploratory analysis of transcriptomic data" "bash 12_diff_analysis" "rnaseq_env"

# --- Final Suggestions and Pipeline End ---
echo "" | tee -a "$LOGFILE_PATH"
echo "=== Pipeline completed successfully ===" | tee -a "$LOGFILE_PATH"
echo "Next steps and tips:" | tee -a "$LOGFILE_PATH"
echo " - Use the Excel and TSV files generated for detailed analysis of identities and regions." | tee -a "$LOGFILE_PATH"
echo " - Download additional genomes or annotation files from:" | tee -a "$LOGFILE_PATH"
echo "   - https://www.ncbi.nlm.nih.gov/genome/" | tee -a "$LOGFILE_PATH"
echo "   - https://www.ncbi.nlm.nih.gov/assembly/" | tee -a "$LOGFILE_PATH"
echo " - Visualize your genomic data with tools like IGV." | tee -a "$LOGFILE_PATH"
echo "" | tee -a "$LOGFILE_PATH"
echo "Log file saved as: $LOGFILE_PATH" | tee -a "$LOGFILE_PATH"
echo "Happy analyzing!" | tee -a "$LOGFILE_PATH"
