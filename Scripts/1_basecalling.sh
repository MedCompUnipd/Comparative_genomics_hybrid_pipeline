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
conda activate preprocessing_env || { echo "[ERROR] Failed to activate 'preprocessing_env'. Please check your Conda environments."; exit 1; }

echo "--- Dorado Basecalling and Trimming Script ---"
echo "----------------------------------------------"

# === USER INPUT ===
read -p "Do you want to convert .fast5 to .pod5? (y/n): " CONVERT_CHOICE

FAST5_DIR=""
if [[ "$CONVERT_CHOICE" =~ ^[Yy]$ ]]; then
    read -p "Directory containing .fast5 files: " FAST5_DIR
    if [[ ! -d "$FAST5_DIR" ]]; then
        echo "[ERROR] FAST5 directory not found: $FAST5_DIR. Exiting."
        exit 1
    fi
fi

read -p "Directory for .pod5 files (output or existing): " POD5_DIR
read -p "Output directory for FASTQ files: " OUT_DIR

THREADS=8 # Default threads, can be made interactive if needed
LOG_DIR="${OUT_DIR}/logs"

# === SETUP ===
echo "[INFO] Creating output directories..."
mkdir -p "$POD5_DIR" "$OUT_DIR" "$LOG_DIR"

# === [1/3] OPTIONAL: CONVERT FAST5 → POD5 ===
if [[ "$CONVERT_CHOICE" =~ ^[Yy]$ ]]; then
    echo "[1/3] Converting .fast5 to .pod5..."

    shopt -s nullglob
    FAST5_FILES=("$FAST5_DIR"/*.fast5)
    shopt -u nullglob

    if [ ${#FAST5_FILES[@]} -eq 0 ]; then
        echo "[ERROR] No .fast5 files found in $FAST5_DIR. Aborting conversion."
        exit 1
    fi

    for file in "${FAST5_FILES[@]}"; do
        base_name=$(basename "$file" .fast5)
        output_file="${POD5_DIR}/${base_name}.pod5"
        echo "→ Converting $file to $output_file..."
        pod5 convert from_fast5 "$file" -o "$output_file" || { echo "[ERROR] pod5 conversion failed for $file."; exit 1; }
    done

    echo "Conversion complete: ${#FAST5_FILES[@]} files processed."
else
    echo "[1/3] Skipping conversion step. Using existing .pod5 files."
fi

# === [2/3] DORADO BASECALLING ===
echo "[2/3] Running Dorado basecalling..."
echo "To perform basecalling, Dorado needs a basecalling model."
read -p "Please insert the path of the directory for Dorado basecalling models download: " DORADO_MODELS_DIR
mkdir -p "$DORADO_MODELS_DIR"

pushd "$DORADO_MODELS_DIR" > /dev/null
# Check if dorado command exists before trying to download models
if ! command_exists dorado; then
    echo "[ERROR] 'dorado' command not found. Please ensure Dorado is installed and in your PATH."
    popd > /dev/null
    exit 1
fi
dorado download --model all || { echo "[ERROR] Dorado model download failed."; popd > /dev/null; exit 1; }
echo "Available Dorado basecalling models:"
ls -1
read -p "Insert path to the Dorado basecalling model (e.g. /home/usr/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0): " MODEL_PATH
popd > /dev/null

shopt -s nullglob
POD5_FILES=("$POD5_DIR"/*.pod5)
shopt -u nullglob

if [ ${#POD5_FILES[@]} -eq 0 ]; then
    echo "[ERROR] No .pod5 files found in $POD5_DIR. Aborting basecalling."
    exit 1
fi

for pod5_file in "${POD5_FILES[@]}"; do
    base_name=$(basename "$pod5_file" .pod5)
    dorado_out_fastq="${OUT_DIR}/${base_name}.fastq"
    dorado_log="${LOG_DIR}/${base_name}_dorado.log"
    echo "→ Basecalling $pod5_file. Output to $dorado_out_fastq"
    dorado basecaller "$MODEL_PATH" "$pod5_file" --threads "$THREADS" > "$dorado_out_fastq" 2> "$dorado_log" || { echo "[ERROR] Dorado basecalling failed for $pod5_file. Check $dorado_log"; exit 1; }
done
echo "Dorado basecalling completed."

# === [3/3] OPTIONAL: ADAPTER TRIMMING with Nanoq (if fastq files exist) ===
echo "[3/3] Performing adapter trimming with Nanoq..."

shopt -s nullglob
FASTQ_FILES=("$OUT_DIR"/*.fastq)
shopt -u nullglob

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "[WARNING] No FASTQ files found in $OUT_DIR for trimming. Skipping trimming step."
else
    TRIMMED_DIR="${OUT_DIR}/trimmed"
    mkdir -p "$TRIMMED_DIR"
    for fastq_file in "${FASTQ_FILES[@]}"; do
        base_name=$(basename "$fastq_file" .fastq)
        trimmed_fastq="${TRIMMED_DIR}/${base_name}_trimmed.fastq"
        echo "→ Trimming adapters from $fastq_file to $trimmed_fastq"
        nanoq trim -i "$fastq_file" -o "$trimmed_fastq" || { echo "[ERROR] Nanoq trimming failed for $fastq_file."; exit 1; }
    done
    echo "Nanoq adapter trimming completed."
    echo "[INFO] Trimmed FASTQ files are in $TRIMMED_DIR"
fi

echo ""
echo "----------------------------------------------"
echo "Dorado basecalling and trimming pipeline completed successfully!"
echo "Raw FASTQ files (if generated): $OUT_DIR"
echo "Trimmed FASTQ files (if generated): $TRIMMED_DIR"
