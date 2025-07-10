#!/bin/bash
set -euo pipefail
shopt -s nullglob

# === LOAD CONDA (portable) ===
# Verify Conda installation and initialize it.
if ! command -v conda &> /dev/null; then
  echo "[ERROR] Conda not found. Make sure conda is installed and initialized."
  exit 1
fi
eval "$(conda shell.bash hook)"

# ---

# === FUNCTION: Check if step should run ===
# This function checks if an output file already exists for a given step.
# If it exists, it asks the user whether to skip the step or overwrite it.
should_run_step() {
  local output_check="$1"
  local step_name="$2"
  local step_dir
  step_dir=$(dirname "$output_check")

  if [[ -e "$output_check" ]]; then
    echo "[INFO] Output for '$step_name' already exists at: $output_check"
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
  else
    mkdir -p "$step_dir"
    return 0
  fi
}

# ---
echo "--- Variant calling script ---"
echo "------------------------------"
# === USER INPUT ===
# Prompt the user for necessary file paths and parameters.
read -r -p "Enter the path to the reference genome FASTA file: " REFERENCE
read -r -p "Enter the path to the de novo genome FASTA file (your polished assembly): " DENOVO_GENOME
read -r -p "Enter the desired output directory path: " OUTDIR
read -r -p "Enter the number of threads to use: " THREADS

# ---

# === CHECK FILES ===
# Verify that all input files exist before proceeding.


for FILE in "${REQUIRED_FILES[@]}"; do
  [[ -f "$FILE" ]] || { echo "[ERROR] Missing input: $FILE"; exit 1; }
done
# Create the main output directory.
mkdir -p "$OUTDIR"
echo "[INFO] Starting the comparison pipeline between your polished de novo assembly ($DENOVO_GENOME) and the reference genome ($REFERENCE)."

# ---

# === STEP 0: MUMmer Comparison ===
# This step compares the large-scale structure of your polished de novo genome with the reference genome.
# Useful for identifying large rearrangements, inversions, or translocations.
MUMMER_OUTPUT_CHECK="$OUTDIR/aln_mummer/ref_vs_denovo.coords"
if should_run_step "$MUMMER_OUTPUT_CHECK" "MUMmer Comparison"; then
  echo "[0] Running MUMmer comparison (reference vs your de novo assembly)..."
  conda activate mummer_env || { echo "[ERROR] Could not activate 'mummer_env'. Check your Conda environments."; exit 1; }
  nucmer --prefix="$OUTDIR/aln_mummer/ref_vs_denovo" "$DENOVO_GENOME" "$REFERENCE" || { echo "[ERROR] nucmer failed."; exit 1; }
  delta-filter -1 "$OUTDIR/aln_mummer/ref_vs_denovo.delta" > "$OUTDIR/aln_mummer/ref_vs_denovo.filtered.delta" || { echo "[ERROR] delta-filter failed."; exit 1; }
  show-coords -rclT "$OUTDIR/aln_mummer/ref_vs_denovo.filtered.delta" > "$OUTDIR/aln_mummer/ref_vs_denovo.tsv" || { echo "[ERROR] show-coords failed."; exit 1; }
  
  # Check if gnuplot is available before trying mummerplot
  if ! command -v gnuplot &> /dev/null; then
      echo "[WARNING] 'gnuplot' not found. MUMmer plots cannot be generated. Install gnuplot manually (e.g., sudo apt install gnuplot)."
  else
      mummerplot --png --layout --prefix "$OUTDIR/aln_mummer/ref_vs_denovo" "$OUTDIR/aln_mummer/ref_vs_denovo.filtered.delta" || { echo "[ERROR] mummerplot failed."; }
  fi
  
  conda deactivate
  echo "[INFO] MUMmer comparison results saved in $OUTDIR/aln_mummer/"
else
  echo "[0] Skipping MUMmer comparison."
fi

# ---

# === STEP 1: Variant Calling (De novo assembly vs Reference) using BCFtools ===
# This step identifies single nucleotide variants (SNPs) and small indels between your polished de novo genome
# and the reference genome using BCFtools.
BCFTOOLS_VCF_OUTPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref.vcf.gz"
if should_run_step "$BCFTOOLS_VCF_OUTPUT" "BCFtools-based Variant Calling (De novo vs Ref)"; then
  echo "[1] Calling variants from your de novo assembly vs reference using BCFtools..."
  conda activate preprocessing_env || { echo "[ERROR] Could not activate 'preprocessing_env'. Check your Conda environments."; exit 1; }
  minimap2 -ax asm5 "$REFERENCE" "$DENOVO_GENOME" > "$OUTDIR/variant_calling_bcftools/aln.sam" || { echo "[ERROR] minimap2 failed."; exit 1; }
  conda deactivate

  conda activate qc_env || { echo "[ERROR] Could not activate 'qc_env'. Check your Conda environments."; exit 1; }
  samtools sort -@ "$THREADS" -o "$OUTDIR/variant_calling_bcftools/aln.sorted.bam" "$OUTDIR/variant_calling_bcftools/aln.sam" || { echo "[ERROR] samtools sort failed."; exit 1; }
  samtools index "$OUTDIR/variant_calling_bcftools/aln.sorted.bam" || { echo "[ERROR] samtools index failed."; exit 1; }
  rm "$OUTDIR/variant_calling_bcftools/aln.sam"

  bcftools mpileup -Ou -f "$REFERENCE" "$OUTDIR/variant_calling_bcftools/aln.sorted.bam" | \
    bcftools call -mv -Oz -o "$BCFTOOLS_VCF_OUTPUT" || { echo "[ERROR] bcftools mpileup/call failed."; exit 1; }
  bcftools index "$BCFTOOLS_VCF_OUTPUT" || { echo "[ERROR] bcftools index failed."; exit 1; }
  conda deactivate
  echo "[INFO] BCFtools-based variants de novo vs Reference saved in $BCFTOOLS_VCF_OUTPUT"
else
  echo "[1] Skipping BCFtools-based variant calling (De novo vs Ref)."
fi

# ---

# === STEP 2: Filter BCFtools Variants for High Confidence (QUAL >= 30) ===
# Filter the variants called by BCFtools (de novo vs reference) to keep only high-confidence calls.
FILTERED_VCF_OUTPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref_filtered.vcf.gz"
if should_run_step "$FILTERED_VCF_OUTPUT" "Filter BCFtools-based Variants"; then
  echo "[2] Filtering BCFtools-based variants de novo vs reference (QUAL >= 30) for high confidence..."
  VCF_INPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref.vcf.gz"
  VCF_TEMP_OUTPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref_filtered.vcf"
  conda activate qc_env || { echo "[ERROR] Could not activate 'qc_env'. Check your Conda environments."; exit 1; }
  bcftools filter -e 'QUAL<30' "$VCF_INPUT" -o "$VCF_TEMP_OUTPUT" -Ov || { echo "[ERROR] bcftools filter failed."; exit 1; }
  bgzip -c "$VCF_TEMP_OUTPUT" > "$FILTERED_VCF_OUTPUT" || { echo "[ERROR] bgzip failed."; exit 1; }
  bcftools index "$FILTERED_VCF_OUTPUT" || { echo "[ERROR] bcftools index failed."; exit 1; }
  rm "$VCF_TEMP_OUTPUT"
  conda deactivate
  echo "[INFO] High-confidence filtered variants de novo vs Reference: $FILTERED_VCF_OUTPUT"
else
  echo "[2] Skipping filtering of BCFtools-based variants."
fi

# ---

echo "-----------------------------------------------------------------"
echo " De novo vs reference comparison pipeline completed successfully."
echo "--- Key comparison outputs ---"
echo "-> Structural differences (MUMmer plots and coordinates): $OUTDIR/aln_mummer/"

echo ""
echo "This script provides a robust comparison of your polished de novo assembly against the reference genome."
