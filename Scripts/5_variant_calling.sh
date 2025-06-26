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

# === USER INPUT ===
# Prompt the user for necessary file paths and parameters.
read -r -p "Enter the path to the reference genome FASTA file: " REFERENCE
read -r -p "Enter the path to the de novo genome FASTA file (your polished assembly): " DENOVO_GENOME
read -r -p "Enter the desired output directory path: " OUTDIR
read -r -p "Enter the number of threads to use: " THREADS

# ---
# Conditional annotation input
PERFORM_ANNOTATION="n" # Default to no annotation
read -r -p "Do you want to perform variant annotation? (y/n): " ANNOTATION_CHOICE
if [[ "$ANNOTATION_CHOICE" =~ ^[SsYy]$ ]]; then
  PERFORM_ANNOTATION="y"
  read -r -p "Enter the path to the GFF3 annotation file (used by bcftools csq): " GFF
  echo "[TIP] You can find GFF3 files from the EMBL-EBI ENA database (https://www.ebi.ac.uk/ena/). Search using the accession number of your reference genome and download the GFF3 format."
else
  echo "[INFO] Variant annotation will be skipped."
  # Set dummy value for GFF if annotation is skipped
  GFF=""
fi

# ---

# === CHECK FILES ===
# Verify that all input files exist before proceeding.
REQUIRED_FILES=("$REFERENCE" "$DENOVO_GENOME")
if [[ "$PERFORM_ANNOTATION" == "y" ]]; then
  REQUIRED_FILES+=("$GFF")
fi

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
  nucmer --prefix="$OUTDIR/aln_mummer/ref_vs_denovo" "$REFERENCE" "$DENOVO_GENOME" || { echo "[ERROR] nucmer failed."; exit 1; }
  delta-filter -1 "$OUTDIR/aln_mummer/ref_vs_denovo.delta" > "$OUTDIR/aln_mummer/ref_vs_denovo.filtered.delta" || { echo "[ERROR] delta-filter failed."; exit 1; }
  show-coords -rcl "$OUTDIR/aln_mummer/ref_vs_denovo.filtered.delta" > "$OUTDIR/aln_mummer/ref_vs_denovo.coords" || { echo "[ERROR] show-coords failed."; exit 1; }
  
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

# === STEP 3: Annotate Filtered Variants using BCFtools csq (Conditional) ===
ANNOTATED_VCF_OUTPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref_filtered_annotated_bcftools.vcf"
if [[ "$PERFORM_ANNOTATION" == "y" ]]; then
  if should_run_step "$ANNOTATED_VCF_OUTPUT" "Annotation of Filtered De novo vs Reference Variants (BCFtools csq)"; then
    echo "[3] Annotating filtered de novo vs reference variants using BCFtools csq..."
    conda activate qc_env || { echo "[ERROR] Could not activate 'qc_env'. Check your Conda environments."; exit 1; }

    # Ensure the reference FASTA is indexed for bcftools csq
    if [[ ! -f "$REFERENCE.fai" ]]; then
      echo "[INFO] Creating FASTA index for the reference genome ($REFERENCE.fai)..."
      samtools faidx "$REFERENCE" || { echo "[ERROR] samtools faidx failed."; exit 1; }
    fi

    bcftools csq -f "$REFERENCE" -g "$GFF" -o "$ANNOTATED_VCF_OUTPUT" "$FILTERED_VCF_OUTPUT" || { echo "[ERROR] bcftools csq annotation failed."; exit 1; }
    conda deactivate
    echo "[INFO] Annotated filtered variants de novo vs Reference (BCFtools csq): $ANNOTATED_VCF_OUTPUT"
  else
    echo "[3] Skipping annotation of filtered variants (BCFtools csq)."
  fi
else
  echo "[INFO] Variant annotation was skipped as requested."
fi

# ---

echo ""
echo "[DONE] De novo vs reference comparison pipeline completed successfully."
echo "--- Key comparison outputs ---"
echo "-> Structural differences (MUMmer plots and coordinates): $OUTDIR/aln_mummer/"
if [[ "$PERFORM_ANNOTATION" == "y" ]]; then
  echo "-> High-confidence point mutations (de novo vs Reference, BCFtools-based, annotated with BCFtools csq): $ANNOTATED_VCF_OUTPUT"
else
  echo "-> High-confidence point mutations (de novo vs Reference, BCFtools-based, NOT ANNOTATED): $FILTERED_VCF_OUTPUT"
fi
echo ""
echo "This script provides a robust comparison of your polished de novo assembly against the reference genome."
