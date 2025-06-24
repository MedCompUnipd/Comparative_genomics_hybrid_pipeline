#!/bin/bash
set -euo pipefail
shopt -s nullglob

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
  echo "[ERROR] Conda not found. Please ensure conda is installed and initialized."
  exit 1
fi
eval "$(conda shell.bash hook)"

echo "--- De Novo Assembly to Reference Comparison Script ---"
echo "-----------------------------------------------------"

# === FUNCTION: Check if step should run ===
# This function checks if an output file or directory already exists for a given step.
# If it exists, it asks the user whether to skip the step or overwrite it.
should_run_step() {
  local output_check="$1"
  local step_name="$2"
  local choice=""

  if [[ -e "$output_check" ]]; then
    echo "[INFO] Output for '$step_name' already exists at: $output_check"
    read -r -p "Do you want to skip this step? (y to skip / n to overwrite): " choice
    if [[ "$choice" =~ ^[Yy]$ ]]; then
      echo "[INFO] Skipping $step_name..."
      return 1 # Skip the step
    else
      echo "[INFO] Overwriting $step_name: removing $output_check..."
      rm -rf "$output_check" # Remove the conflicting output
      return 0
    fi
  else
    return 0
  fi
}

# === USER INPUT ===
read -r -p "Enter path to the de novo assembly FASTA file (e.g., polished.fasta): " DENOVO_ASSEMBLY
read -r -p "Enter path to the reference genome FASTA file: " REFERENCE_GENOME
read -r -p "Enter path to the reference genome GFF/GTF file (for annotation, optional): " REF_GFF
read -r -p "Enter desired output directory for comparison results: " OUTPUT_DIR
read -r -p "Enter threads to use: " THREADS

# --- Validate Inputs ---
if [[ ! -f "$DENOVO_ASSEMBLY" ]]; then
    echo "[ERROR] De novo assembly FASTA not found: $DENOVO_ASSEMBLY. Exiting."
    exit 1
fi
if [[ ! -f "$REFERENCE_GENOME" ]]; then
    echo "[ERROR] Reference genome FASTA not found: $REFERENCE_GENOME. Exiting."
    exit 1
fi
mkdir -p "$OUTPUT_DIR" || { echo "[ERROR] Failed to create output directory: $OUTPUT_DIR"; exit 1; }

# === CONDA ENVIRONMENT NAMES ===
QC_ENV="qc_env"

# --- MAIN WORKFLOW ---
echo ""
echo "--- Starting Comparison Pipeline ---"

# === STEP 1: QUAST Comparison ===
QUAST_COMP_DIR="${OUTPUT_DIR}/quast_comparison"
if should_run_step "$QUAST_COMP_DIR/report.html" "QUAST Comparison"; then
  echo "[1/3] Running QUAST comparison..."
  mkdir -p "$QUAST_COMP_DIR"
  conda activate "$QC_ENV" || { echo "[ERROR] Failed to activate 'qc_env'. Please check your Conda environments."; exit 1; }
  quast.py -r "$REFERENCE_GENOME" -t "$THREADS" -o "$QUAST_COMP_DIR" "$DENOVO_ASSEMBLY" || { echo "[ERROR] QUAST comparison failed."; exit 1; }
  echo "[INFO] QUAST comparison completed. Reports in $QUAST_COMP_DIR"
  conda deactivate
else
  echo "[1/3] Skipping QUAST Comparison."
fi

# === STEP 2: MUMmer Alignment and Plotting ===
MUMMER_ALIGN_DIR="${OUTPUT_DIR}/aln_mummer"
if should_run_step "$MUMMER_ALIGN_DIR" "MUMmer Alignment and Plotting"; then
  echo "[2/3] Running MUMmer alignment and plotting..."
  mkdir -p "$MUMMER_ALIGN_DIR"
  conda activate "$QC_ENV" || { echo "[ERROR] Failed to activate 'qc_env'. Please check your Conda environments."; exit 1; }

  # 1. Run nucmer
  echo "[INFO] Running nucmer..."
  nucmer --prefix="${MUMMER_ALIGN_DIR}/denovo_vs_ref" "$REFERENCE_GENOME" "$DENOVO_ASSEMBLY" || { echo "[ERROR] nucmer failed."; exit 1; }

  # 2. Filter alignments
  echo "[INFO] Filtering delta file with delta-filter..."
  delta-filter -q "${MUMMER_ALIGN_DIR}/denovo_vs_ref.delta" > "${MUMMER_ALIGN_DIR}/denovo_vs_ref.filtered.delta" || { echo "[ERROR] delta-filter failed."; exit 1; }

  # 3. Generate alignment coordinates
  echo "[INFO] Generating alignment coordinates with show-coords..."
  show-coords -r -c -l -T "${MUMMER_ALIGN_DIR}/denovo_vs_ref.filtered.delta" > "${MUMMER_ALIGN_DIR}/denovo_vs_ref.coords" || { echo "[ERROR] show-coords failed."; exit 1; }

  # 4. Generate plots (requires gnuplot, which is usually a system dependency)
  echo "[INFO] Generating alignment plots with mummerplot..."
  # Check if gnuplot is available
  if ! command_exists gnuplot; then
      echo "[WARNING] 'gnuplot' not found. MUMmer plots cannot be generated. Please install gnuplot manually (e.g., sudo apt install gnuplot)."
  else
      # Using -fatb for larger dots, -x and -y for axis ranges from assembly sizes if known
      mummerplot --postscript --prefix="${MUMMER_ALIGN_DIR}/denovo_vs_ref" -p "${MUMMER_ALIGN_DIR}/denovo_vs_ref" "${MUMMER_ALIGN_DIR}/denovo_vs_ref.filtered.delta" || { echo "[ERROR] mummerplot failed."; } # mummerplot often fails gracefully
      # Convert PostScript to PDF for easier viewing
      if [[ -f "${MUMMER_ALIGN_DIR}/denovo_vs_ref.ps" ]]; then
          echo "[INFO] Converting PostScript plot to PDF..."
          ps2pdf "${MUMMER_ALIGN_DIR}/denovo_vs_ref.ps" "${MUMMER_ALIGN_DIR}/denovo_vs_ref.pdf" || echo "[WARNING] ps2pdf conversion failed. PDF plot might not be generated."
      fi
  fi
  echo "[INFO] MUMmer alignment completed. Outputs in $MUMMER_ALIGN_DIR"
  conda deactivate
else
  echo "[2/3] Skipping MUMmer Alignment and Plotting."
fi

# === STEP 3: Variant Calling with BCFtools and Annotation ===
BCFTOOLS_VARIANT_DIR="${OUTPUT_DIR}/variant_calling_bcftools"
if should_run_step "$BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered_annotated_bcftools.vcf" "BCFtools Variant Calling and Annotation"; then
  echo "[3/3] Running BCFtools variant calling and optional annotation..."
  mkdir -p "$BCFTOOLS_VARIANT_DIR"
  conda activate "$QC_ENV" || { echo "[ERROR] Failed to activate 'qc_env'. Please check your Conda environments."; exit 1; }

  # 1. Index reference genome for BWA if not already done
  if [[ ! -f "$REFERENCE_GENOME.fai" ]]; then
      echo "[INFO] Indexing reference genome for BWA (required for BCFtools mpileup)..."
      bwa index "$REFERENCE_GENOME" || { echo "[ERROR] BWA indexing of reference failed."; exit 1; }
  fi

  # 2. Map de novo assembly to reference using BWA-MEM
  echo "[INFO] Mapping de novo assembly to reference with BWA-MEM..."
  local BWA_BAM="${BCFTOOLS_VARIANT_DIR}/denovo_vs_ref.bam"
  bwa mem -t "$THREADS" "$REFERENCE_GENOME" "$DENOVO_ASSEMBLY" | samtools view -Sb - | samtools sort -o "$BWA_BAM" - || { echo "[ERROR] BWA mapping/Samtools sort failed."; exit 1; }
  samtools index "$BWA_BAM" || { echo "[ERROR] Samtools indexing failed."; exit 1; }

  # 3. Call variants with bcftools mpileup and call
  echo "[INFO] Calling variants with BCFtools mpileup and call..."
  bcftools mpileup -f "$REFERENCE_GENOME" "$BWA_BAM" | bcftools call -mv -o "$BCFTOOLS_VARIANT_DIR/denovo_vs_ref.vcf" || { echo "[ERROR] BCFtools mpileup/call failed."; exit 1; }
  echo "[INFO] Variants called: $BCFTOOLS_VARIANT_DIR/denovo_vs_ref.vcf"

  # 4. Filter variants (e.g., remove low quality)
  echo "[INFO] Filtering variants..."
  bcftools view -f PASS "$BCFTOOLS_VARIANT_DIR/denovo_vs_ref.vcf" | bcftools norm -d all -f "$REFERENCE_GENOME" -o "$BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered.vcf" || { echo "[ERROR] BCFtools filtering failed."; exit 1; }
  bgzip -f "$BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered.vcf" || { echo "[ERROR] bgzip failed."; exit 1; }
  tabix -p vcf "$BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered.vcf.gz" || { echo "[ERROR] tabix failed."; exit 1; }
  echo "[INFO] Filtered variants: $BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered.vcf.gz"

  # 5. Optional: Annotate variants with BCFtools csq
  local PERFORM_ANNOTATION="n"
  if [[ -f "$REF_GFF" ]]; then
      read -r -p "Reference GFF/GTF file provided. Do you want to perform variant annotation with BCFtools csq? (y/n): " PERFORM_ANNOTATION
  fi

  if [[ "$PERFORM_ANNOTATION" =~ ^[Yy]$ ]]; then
    if [[ ! -f "$REF_GFF" ]]; then
      echo "[WARNING] GFF/GTF file not found: $REF_GFF. Skipping BCFtools csq annotation."
    else
      echo "[INFO] Annotating variants using BCFtools csq..."
      # Ensure the reference FASTA is indexed for bcftools csq
      if [[ ! -f "$REFERENCE_GENOME.fai" ]]; then
        echo "[INFO] Creating FASTA index for reference genome ($REFERENCE_GENOME.fai)..."
        samtools faidx "$REFERENCE_GENOME" || { echo "[ERROR] Samtools faidx failed."; exit 1; }
      fi
      # Annotate with bcftools csq
      bcftools csq -f "$REFERENCE_GENOME" -g "$REF_GFF" -o "$BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered_annotated_bcftools.vcf" "$BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered.vcf.gz" || { echo "[ERROR] BCFtools csq annotation failed."; exit 1; }
      echo "[INFO] Annotated filtered de novo vs Reference variants (BCFtools csq): $BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered_annotated_bcftools.vcf"
    fi
  else
    echo "[INFO] Skipping variant annotation as requested."
  fi
  conda deactivate
else
  echo "[3/3] Skipping BCFtools Variant Calling and Annotation."
fi

echo ""
echo "[DONE] De novo assembly to Reference comparison pipeline completed successfully."
echo "--- Key Comparison Outputs ---"
echo "-> QUAST reports: $QUAST_COMP_DIR/"
echo "-> Structural differences (MUMmer plots and coords): $MUMMER_ALIGN_DIR/"
echo "-> High-confidence point mutations (BCFtools-based): ${BCFTOOLS_VARIANT_DIR}/denovo_vs_ref_filtered.vcf.gz"
if [[ "$PERFORM_ANNOTATION" =~ ^[Yy]$ && -f "$BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered_annotated_bcftools.vcf" ]]; then
  echo "-> Annotated variants (BCFtools csq): $BCFTOOLS_VARIANT_DIR/denovo_vs_ref_filtered_annotated_bcftools.vcf"
fi
