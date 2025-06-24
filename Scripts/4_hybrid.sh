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
  echo "[ERROR] 'conda' not found. Please ensure Conda is installed and initialized."
  exit 1
fi
eval "$(conda shell.bash hook)"

# === SAFE CONDA ACTIVATION FUNCTION ===
safe_activate() {
  local env_name="$1"
  echo "[INFO] Activating Conda environment: $env_name"
  if ! conda activate "$env_name" &>/dev/null; then
    echo "[ERROR] Conda environment '$env_name' not found or failed to activate. Please check your Conda environments."
    exit 1
  fi
}

# === STEP CONTROL FUNCTION ===
# This function checks if an output file already exists for a given step.
# If it exists, it asks the user whether to skip the step or overwrite it.
should_run_step() {
  local output_check="$1"
  local step_name="$2"
  local step_dir
  step_dir=$(dirname "$output_check")

  if [[ -e "$output_check" ]]; then
    echo "[INFO] Output for '$step_name' already exists: $output_check"
    if [[ "$AUTO_MODE" == true ]]; then
      echo "[INFO] AUTO_MODE active: skipping $step_name"
      return 1
    else
      read -p "Do you want to skip this step? (y to skip / n to overwrite): " choice
      [[ "$choice" =~ ^[Yy]$ ]] && return 1
    fi
  fi
  # If not skipping, ensure directory exists (or recreate if overwriting implicitly)
  mkdir -p "$step_dir" || { echo "[ERROR] Failed to create directory: $step_dir"; exit 1; }
  return 0
}

# === USER INPUT ===
read -p "Enter the path to the folder containing Illumina reads (FASTQ.gz or FASTQ): " ILLUMINA_READS_DIR
read -p "Enter the path to the folder containing ONT reads (FASTQ or FASTQ.gz): " ONT_READS_DIR
read -p "Enter the path to the de novo assembly FASTA file (e.g., from Script 3): " DENOVO_ASSEMBLY_FASTA
read -p "Enter the project name (e.g., hybrid_assembly_project): " PROJECT_NAME
read -p "Enter desired root output directory for this script: " OUTPUT_ROOT_DIR
read -p "Enter Medaka model (e.g., r1041_e82_400bps_sup_variant_v4.2.0): " MEDAKA_MODEL
read -p "Enter threads to use: " THREADS

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
mkdir -p "$OUTPUT_ROOT_DIR" || { echo "[ERROR] Failed to create output root directory: $OUTPUT_ROOT_DIR"; exit 1; }

# === CONDA ENVIRONMENT NAMES ===
POLISHING_ENV="polishing_env"
PREPROCESSING_ENV="preprocessing_env"
QC_ENV="qc_env"

# === MAIN WORKFLOW ===
LOCAL_OUTPUT_DIR="${OUTPUT_ROOT_DIR}/${PROJECT_NAME}_hybrid_polishing"
mkdir -p "$LOCAL_OUTPUT_DIR" || { echo "[ERROR] Failed to create local output directory: $LOCAL_OUTPUT_DIR"; exit 1; }

echo ""
echo "--- Starting Hybrid Polishing Pipeline ---"

# === STEP 1: Merge Illumina Reads ===
ILLUMINA_MERGED_FASTQ_GZ="${LOCAL_OUTPUT_DIR}/1_illumina_merge/${PROJECT_NAME}_illumina_merged.fastq.gz"
if should_run_step "$ILLUMINA_MERGED_FASTQ_GZ" "Illumina Read Merging"; then
    echo "[STEP 1] Merging Illumina reads..."
    mkdir -p "$(dirname "$ILLUMINA_MERGED_FASTQ_GZ")"

    # Find all .fastq or .fastq.gz files in the Illumina directory
    shopt -s nullglob
    ILLUMINA_FILES=("$ILLUMINA_READS_DIR"/*.fastq "$ILLUMINA_READS_DIR"/*.fastq.gz)
    shopt -u nullglob

    if [ ${#ILLUMINA_FILES[@]} -eq 0 ]; then
        echo "[ERROR] No .fastq or .fastq.gz files found in $ILLUMINA_READS_DIR. Aborting Illumina merge."
        exit 1
    fi

    cat "${ILLUMINA_FILES[@]}" > "$ILLUMINA_MERGED_FASTQ_GZ" || { echo "[ERROR] Failed to merge Illumina reads."; exit 1; }
    echo "[INFO] Illumina reads merged to: $ILLUMINA_MERGED_FASTQ_GZ"
fi

# === STEP 2: Pilon Polishing ===
PILON_POLISHED_FASTA="${LOCAL_OUTPUT_DIR}/2_pilon_pol/pilon_polished.fasta"
if should_run_step "$PILON_POLISHED_FASTA" "Pilon Polishing"; then
    echo "[STEP 2] Running Pilon polishing with Illumina reads..."
    mkdir -p "$(dirname "$PILON_POLISHED_FASTA")"

    safe_activate "$POLISHING_ENV"

    # 1. Map Illumina reads to the de novo assembly using BWA-MEM
    echo "[INFO] Mapping Illumina reads to de novo assembly with BWA-MEM..."
    local BWA_BAM="${LOCAL_OUTPUT_DIR}/2_pilon_pol/illumina_mapped.bam"
    bwa index "$DENOVO_ASSEMBLY_FASTA" || { echo "[ERROR] BWA indexing failed for de novo assembly."; exit 1; }
    bwa mem -t "$THREADS" "$DENOVO_ASSEMBLY_FASTA" "$ILLUMINA_MERGED_FASTQ_GZ" | samtools view -Sb - | samtools sort -o "$BWA_BAM" - || { echo "[ERROR] BWA mapping/Samtools sort failed."; exit 1; }
    samtools index "$BWA_BAM" || { echo "[ERROR] Samtools indexing failed for Illumina BAM."; exit 1; }
    echo "[INFO] Illumina reads mapped to: $BWA_BAM"

    # 2. Run Pilon
    echo "[INFO] Running Pilon..."
    # Pilon often requires Java. openjdk=8 should be in polishing_env
    java -Xmx$(echo "$THREADS * 2" | bc)G -jar "$CONDA_PREFIX/share/pilon/pilon.jar" --genome "$DENOVO_ASSEMBLY_FASTA" --frags "$BWA_BAM" --output pilon_polished --outdir "$(dirname "$PILON_POLISHED_FASTA")" --threads "$THREADS" || { echo "[ERROR] Pilon polishing failed."; exit 1; }
    echo "[INFO] Pilon polishing completed. Output: $PILON_POLISHED_FASTA"
    conda deactivate
fi

# === STEP 3: Medaka Polishing ===
MEDAKA_POLISHED_FASTA="${LOCAL_OUTPUT_DIR}/3_medaka_pol/medaka_polished.fasta"
if should_run_step "$MEDAKA_POLISHED_FASTA" "Medaka Polishing"; then
    echo "[STEP 3] Running Medaka polishing with ONT reads..."
    mkdir -p "$(dirname "$MEDAKA_POLISHED_FASTA")"

    safe_activate "$POLISHING_ENV"

    # Merge ONT reads if not a single file (similar to 2_tax_sel_qc.sh)
    ONT_MERGED_FASTQ="${LOCAL_OUTPUT_DIR}/3_medaka_pol/${PROJECT_NAME}_ont_merged.fastq"
    find "$ONT_READS_DIR" -type f -name "*.fastq" -print0 | xargs -0 cat > "$ONT_MERGED_FASTQ" || \
    find "$ONT_READS_DIR" -type f -name "*.fastq.gz" -print0 | xargs -0 gzip -cd >> "$ONT_MERGED_FASTQ" || { echo "[ERROR] Failed to merge ONT reads."; exit 1; }

    if [ ! -s "$ONT_MERGED_FASTQ" ]; then
        echo "[ERROR] No ONT FASTQ files found or merged file is empty. Aborting Medaka."
        exit 1
    fi

    medaka_consensus -i "$ONT_MERGED_FASTQ" -d "$PILON_POLISHED_FASTA" -o "$(dirname "$MEDAKA_POLISHED_FASTA")" -t "$THREADS" -m "$MEDAKA_MODEL" || { echo "[ERROR] Medaka polishing failed."; exit 1; }
    # Medaka output is consensus.fasta within the output directory
    mv "$(dirname "$MEDAKA_POLISHED_FASTA")/consensus.fasta" "$MEDAKA_POLISHED_FASTA" || { echo "[ERROR] Failed to move Medaka output consensus.fasta."; exit 1; }
    echo "[INFO] Medaka polishing completed. Output: $MEDAKA_POLISHED_FASTA"
    conda deactivate
fi

# === STEP 4: QUAST QC for Hybrid Assembly ===
HYBRID_QUAST_REPORT_DIR="${LOCAL_OUTPUT_DIR}/4_quast"
if should_run_step "$HYBRID_QUAST_REPORT_DIR/report.html" "Hybrid Assembly QUAST QC"; then
  echo "[STEP 4] Running QUAST for hybrid assembly..."
  safe_activate "$QC_ENV"
  quast.py -t "$THREADS" -o "$HYBRID_QUAST_REPORT_DIR" "$MEDAKA_POLISHED_FASTA" || { echo "[ERROR] QUAST failed for hybrid assembly."; exit 1; }
  echo "[INFO] QUAST QC for hybrid assembly completed. Reports in $HYBRID_QUAST_REPORT_DIR"
  conda deactivate
fi

# === STEP 5: Mapping with minimap2 and Qualimap (for the final hybrid assembly) ===
HYBRID_MAPPING_DIR="${LOCAL_OUTPUT_DIR}/5_hybrid_mapping"
if should_run_step "$HYBRID_MAPPING_DIR/qualimap_report/genome_results.txt" "Hybrid Mapping and Qualimap"; then
  echo "[STEP 5] Mapping ONT reads to hybrid assembly with minimap2 and running Qualimap..."
  mkdir -p "$HYBRID_MAPPING_DIR"

  safe_activate "$PREPROCESSING_ENV"
  ONT_MERGED_FASTQ="${LOCAL_OUTPUT_DIR}/3_medaka_pol/${PROJECT_NAME}_ont_merged.fastq" # Re-use merged ONT reads

  if [ ! -f "$ONT_MERGED_FASTQ" ]; then
      echo "[ERROR] Merged ONT FASTQ not found: $ONT_MERGED_FASTQ. Aborting mapping/Qualimap."
      exit 1
  fi

  local HYBRID_BAM="${HYBRID_MAPPING_DIR}/${PROJECT_NAME}_hybrid_mapped.bam"
  minimap2 -ax map-ont "$MEDAKA_POLISHED_FASTA" "$ONT_MERGED_FASTQ" | samtools view -bS - | samtools sort -o "$HYBRID_BAM" - || { echo "[ERROR] Minimap2/Samtools mapping failed for hybrid assembly."; exit 1; }
  samtools index "$HYBRID_BAM" || { echo "[ERROR] Samtools indexing failed for hybrid assembly BAM."; exit 1; }
  echo "[INFO] ONT reads mapped to hybrid assembly. BAM file: $HYBRID_BAM"
  conda deactivate

  safe_activate "$QC_ENV"
  qualimap bamqc -bam "$HYBRID_BAM" -outdir "${HYBRID_MAPPING_DIR}/qualimap_report" --java-mem-size=4G || { echo "[ERROR] Qualimap failed for hybrid assembly mapping."; exit 1; }
  echo "[INFO] Qualimap QC for hybrid assembly completed. Reports in ${HYBRID_MAPPING_DIR}/qualimap_report"
  conda deactivate
fi

echo ""
echo "-------------------------------------"
echo "Hybrid polishing pipeline completed successfully!"
echo "Final polished assembly: $MEDAKA_POLISHED_FASTA"
echo "QUAST report: $HYBRID_QUAST_REPORT_DIR"
echo "Qualimap report: ${HYBRID_MAPPING_DIR}/qualimap_report"
