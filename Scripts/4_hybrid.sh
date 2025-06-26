#!/bin/bash
set -euo pipefail
shopt -s nullglob

# === CONFIG: NON-INTERACTIVE MODE ===
AUTO_MODE=false  # Set to true to skip interactive prompts for skipping steps

# === SAVE ORIGINAL DIRECTORY ===
ORIG_PWD="$(pwd)"
# Trap errors to return to original directory and exit cleanly on failure
trap 'cd "$ORIG_PWD" || exit; echo "[ERROR] Script interrupted or failed. Exiting."; exit 1' ERR



# --- Conda setup ---
# Function to check if a command exists
command_exists () {
    command -v "$1" &> /dev/null
}

if ! command_exists conda; then
  echo "[ERROR] 'conda' not found. Make sure Conda is installed and initialized."
  exit 1
fi
eval "$(conda shell.bash hook)"

# === SAFE CONDA ACTIVATION FUNCTION ===
safe_activate() {
  local env_name="$1"
  echo "[INFO] Activating Conda environment: $env_name"
  if ! conda activate "$env_name" &>/dev/null; then
    echo "[ERROR] Conda environment '$env_name' not found or activation failed. Check your Conda environments."
    exit 1
  fi
}

# === STEP CONTROL FUNCTION ===
# This function checks if an output file or directory already exists for a given step.
# If it exists, it asks the user whether to skip the step or overwrite it.
should_run_step() {
  local output_check="$1"
  local step_name="$2"
  local step_dir # This will be the directory where output_check resides, or output_check itself if it's a directory

  if [[ -d "$output_check" ]]; then
    step_dir="$output_check"
  else
    step_dir=$(dirname "$output_check")
  fi

  if [[ -e "$output_check" ]]; then
    echo "[INFO] Output for '$step_name' already exists at: $output_check"
    if [[ "$AUTO_MODE" == true ]]; then
      echo "[INFO] AUTO_MODE enabled: skipping $step_name"
      return 1 # Skip the step
    else
      read -r -p "Do you want to skip this step? (s to skip / n to overwrite): " choice
      if [[ "$choice" =~ ^[Ss]$ ]]; then
        echo "[INFO] Skipping $step_name..."
        return 1 # Skip the step
      else
        echo "[INFO] Overwriting $step_name: removing $step_dir..."
        rm -rf "$step_dir"
        mkdir -p "$step_dir"
        return 0
      fi
    fi
  else
    mkdir -p "$step_dir"
    return 0
  fi
}

echo "--- Integration with Illumina reads ---"
echo "---------------------------------------"

# === USER INPUT ===
read -p "Enter the path to the folder containing Illumina reads (FASTQ.gz): " ILLUMINA_READS_DIR
read -p "Enter the path to the folder containing ONT reads (FASTQ or FASTQ.gz): " ONT_READS_DIR
read -p "Enter the path to the FASTA file of the de novo assembly (e.g. from Script 3): " DENOVO_ASSEMBLY_FASTA
read -p "Enter the project name (e.g. hybrid_assembly_project): " PROJECT_NAME
read -p "Enter the desired root output directory for this script: " OUTPUT_ROOT_DIR
read -p "Enter the path to create/use the Kraken2 database: " DB_PATH
read -p "Enter the TaxID to filter (e.g. 10359 for HCMV): " TAXID
read -p "Enter the Medaka model (e.g. r1041_e82_400bps_sup_variant_v4.2.0): " MEDAKA_MODEL
read -p "Enter the number of threads to use: " THREADS
read -p "Enter the full path to pilon.jar (e.g. /path/to/pilon/pilon.jar): " PILON_JAR
read -p "Enter the amount of RAM to allocate for Java (e.g. 32G): " JAVA_RAM

# --- Validate Inputs ---
if [[ ! -d "$ILLUMINA_READS_DIR" ]]; then
    echo "[ERROR] Illumina reads directory not found: $ILLUMINA_READS_DIR. Exiting."
    exit 1
fi
if [[ ! -d "$ONT_READS_DIR" ]]; then
    echo "[ERROR] ONT reads directory not found: $ONT_READS_DIR. Exiting."
    exit 1
fi
if [[ ! -f "$DENOVO_ASSEMBLY_FASTA" ]]; then
    echo "[ERROR] De novo assembly FASTA file not found: $DENOVO_ASSEMBLY_FASTA. Exiting."
    exit 1
fi
if [[ ! -f "$PILON_JAR" ]]; then
    echo "[ERROR] pilon.jar file not found: $PILON_JAR. Exiting."
    exit 1
fi
mkdir -p "$OUTPUT_ROOT_DIR" || { echo "[ERROR] Could not create root output directory: $OUTPUT_ROOT_DIR"; exit 1; }

# === CONDA ENVIRONMENT NAMES ===
KRAKEN_ENV="kraken2_env"
QC_ENV="qc_env"
ILLUMINA_ENV="illuminareads_env"
POLISHING_ENV="polishing_env"
PREPROCESSING_ENV="preprocessing_env"

# === MAIN WORKFLOW ===
LOCAL_OUTPUT_DIR="${OUTPUT_ROOT_DIR}/${PROJECT_NAME}_hybrid_polishing"
mkdir -p "$LOCAL_OUTPUT_DIR" || { echo "[ERROR] Could not create local output directory: $LOCAL_OUTPUT_DIR"; exit 1; }

echo ""
echo "--- Starting Hybrid Polishing Pipeline ---"

# === STEP 1: KRAKEN2 DB BUILD ===
if should_run_step "$DB_PATH/hash.k2d" "Kraken2 DB Build"; then
  echo "[STEP 1] Building Kraken2 database..."
  safe_activate "$KRAKEN_ENV"
  kraken2-build --standard --db "$DB_PATH" --threads "$THREADS" || { echo "[ERROR] Kraken2 standard DB build failed."; exit 1; }
  kraken2-build --build --db "$DB_PATH" --threads "$THREADS" || { echo "[ERROR] Kraken2 final DB build failed."; exit 1; }
  conda deactivate
  echo "[INFO] Kraken2 database built at: $DB_PATH"
else
  echo "[STEP 1] Skipping Kraken2 Database Build."
fi
# === STEP 2: CONTAMINANT FILTERING ===
FILTERED_DIR="${LOCAL_OUTPUT_DIR}/filtered_taxid_${TAXID}"
if should_run_step "${FILTERED_DIR}/filtered_illumina_R1.fastq.gz" "Contaminant Filtering"; then
  echo "[STEP 2] Filtering reads by taxonomic ID with Kraken2..."
  safe_activate "$QC_ENV"
  
  # Illumina
  kraken2 --db "$DB_PATH" --threads "$THREADS" --report "${FILTERED_DIR}/illumina_kraken_report.txt" \
    --paired "$ILLUMINA_READS_DIR"/*_R1*.fastq.gz "$ILLUMINA_READS_DIR"/*_R2*.fastq.gz \
    --output "${FILTERED_DIR}/illumina_kraken_output.txt" --use-names --gzip-compressed \
    --classified-out "${FILTERED_DIR}/classified_illumina#.fq"
  
  mv "${FILTERED_DIR}/classified_illumina_1.fq" "${FILTERED_DIR}/filtered_illumina_R1.fastq"
  mv "${FILTERED_DIR}/classified_illumina_2.fq" "${FILTERED_DIR}/filtered_illumina_R2.fastq"
  pigz "${FILTERED_DIR}/filtered_illumina_R1.fastq"
  pigz "${FILTERED_DIR}/filtered_illumina_R2.fastq"

  # ONT
  kraken2 --db "$DB_PATH" --threads "$THREADS" --report "${FILTERED_DIR}/ont_kraken_report.txt" \
    "$ONT_READS_DIR"/*.fastq* \
    --output "${FILTERED_DIR}/ont_kraken_output.txt" --use-names --gzip-compressed \
    --classified-out "${FILTERED_DIR}/filtered_ont.fastq"
  
  conda deactivate
  echo "[INFO] Filtering completed. Filtered files in $FILTERED_DIR"
else
  echo "[STEP 2] Skipping contaminant filtering."
fi

# === STEP 3: POLISHING ===
POLISHING_DIR="${LOCAL_OUTPUT_DIR}/polishing"
if should_run_step "${POLISHING_DIR}/final_polished.fasta" "Polishing Pipeline"; then
  echo "[STEP 3] Running hybrid polishing pipeline..."
  safe_activate "$POLISHING_ENV"
  mkdir -p "$POLISHING_DIR"

  # --- Medaka ---
  echo "[INFO] Running Medaka..."
  medaka_consensus -i "${FILTERED_DIR}/filtered_ont.fastq" -d "$DENOVO_ASSEMBLY_FASTA" \
    -o "${POLISHING_DIR}/medaka_output" -t "$THREADS" -m "$MEDAKA_MODEL"

  # --- Align Illumina reads to Medaka-polished assembly ---
  echo "[INFO] Aligning Illumina reads to Medaka output..."
  bwa index "${POLISHING_DIR}/medaka_output/consensus.fasta"
  bwa mem -t "$THREADS" "${POLISHING_DIR}/medaka_output/consensus.fasta" \
    "${FILTERED_DIR}/filtered_illumina_R1.fastq.gz" "${FILTERED_DIR}/filtered_illumina_R2.fastq.gz" |
    samtools sort -@ "$THREADS" -o "${POLISHING_DIR}/illumina_sorted.bam"
  samtools index "${POLISHING_DIR}/illumina_sorted.bam"

  # --- Pilon ---
  echo "[INFO] Running Pilon..."
  java -Xmx"$JAVA_RAM" -jar "$PILON_JAR" --genome "${POLISHING_DIR}/medaka_output/consensus.fasta" \
    --frags "${POLISHING_DIR}/illumina_sorted.bam" --output final_polished --outdir "$POLISHING_DIR" \
    --threads "$THREADS" --changes --vcf || { echo "[ERROR] Pilon failed."; exit 1; }

  conda deactivate
  echo "[INFO] Polishing completed. Final assembly in: ${POLISHING_DIR}/final_polished.fasta"
else
  echo "[STEP 3] Skipping polishing step."
fi

# === FINAL ===
cd "$ORIG_PWD" || exit
echo "--------------------------------------------------------"
echo "--- HYBRID POLISHING PIPELINE COMPLETED SUCCESSFULLY ---"
echo "Output directory: $LOCAL_OUTPUT_DIR"
