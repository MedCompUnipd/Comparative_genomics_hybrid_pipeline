#!/bin/bash
set -euo pipefail
shopt -s nullglob

# === CONFIG: NON-INTERACTIVE MODE ===
AUTO_MODE=false  # Set to true to skip interactive prompts

# === SAVE ORIGINAL DIRECTORY ===
ORIG_PWD="$(pwd)"
trap 'cd "$ORIG_PWD"; echo "[ERROR] Unexpected interruption or failure."; exit 1' ERR

# === REQUIRED TOOLS CHECK ===
required_tools=(conda kraken2 bwa samtools fastqc fastp minimap2 qualimap)
for tool in "${required_tools[@]}"; do
  if ! command -v "$tool" &>/dev/null; then
    echo "[ERROR] Required tool '$tool' not found in PATH."
    exit 1
  fi
done

# === LOAD CONDA ===
if ! command -v conda &> /dev/null; then
  echo "[ERROR] 'conda' not found. Please ensure Conda is installed and initialized."
  exit 1
fi
eval "$(conda shell.bash hook)"

safe_activate() {
  local env="$1"
  if ! conda activate "$env" &>/dev/null; then
    echo "[ERROR] Conda environment '$env' not found or failed to activate."
    exit 1
  fi
}

# === STEP CONTROL FUNCTION ===
should_run_step() {
  local output_check="$1"
  local step_name="$2"

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
  return 0
}
echo "--- Integration with Illumina reads ---"
echo "---------------------------------------"
# === USER INPUT ===
read -p "Enter the path to the folder containing Illumina reads (FASTQ.gz): " REMOTE_DIR
read -p "Enter the path to save results: " LOCAL_OUTPUT_DIR
read -p "Enter the path to create the Kraken2 database: " DB_PATH
read -p "Enter the TaxID to filter for (e.g., 10359): " TAXID
read -p "Enter the path to the base folder containing polished assemblies: " BASE_GENOME_PATH
read -p "Enter the assembler used (e.g., flye, canu): " ASSEMBLER
read -p "Enter number of threads to use: " THREADS

GENOME="${BASE_GENOME_PATH}/${ASSEMBLER}/medaka_output/consensus.fasta"
if [[ ! -f "$GENOME" ]]; then
  echo "[ERROR] Polished genome not found at: $GENOME"
  exit 1
fi

mkdir -p "$LOCAL_OUTPUT_DIR"
echo "[INFO] Using polished genome: $GENOME"

# === STEP 1: KRAKEN2 DB BUILD ===
if should_run_step "$DB_PATH/hash.k2d" "Kraken2 DB Build"; then
  echo "[STEP 1] Building Kraken2 database..."
  safe_activate kraken2_env
  kraken2-build --standard --db "$DB_PATH" --threads "$THREADS"
  kraken2-build --build --db "$DB_PATH" --threads "$THREADS"
  conda deactivate
fi

# === STEP 2: KRAKEN2 CLASSIFICATION ===
if should_run_step "$LOCAL_OUTPUT_DIR/1_kraken" "Kraken2 Classification"; then
  echo "[STEP 2] Classifying reads with Kraken2..."
  safe_activate kraken2_env
  mkdir -p "$LOCAL_OUTPUT_DIR/1_kraken"
  for R1 in "$REMOTE_DIR"/*_R1.fastq.gz; do
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    kraken2 --db "$DB_PATH" --threads "$THREADS" --paired \
      --report "$LOCAL_OUTPUT_DIR/1_kraken/${SAMPLE}.report" \
      --output "$LOCAL_OUTPUT_DIR/1_kraken/${SAMPLE}.kraken" \
      "$R1" "$R2"
  done
  conda deactivate
fi

# === STEP 3: FILTERING READS ===
if should_run_step "$LOCAL_OUTPUT_DIR/2_seqtk_filt" "Read Filtering"; then
  echo "[STEP 3] Filtering reads for TaxID $TAXID..."
  safe_activate kraken2_env
  mkdir -p "$LOCAL_OUTPUT_DIR/2_seqtk_filt"

  USE_SEQTK=false
  if command -v seqtk &> /dev/null; then
    USE_SEQTK=true
    echo "[INFO] Using seqtk for filtering."
  else
    echo "[WARNING] seqtk not found. Falling back to awk filtering."
  fi

  for KRAKEN_FILE in "$LOCAL_OUTPUT_DIR/1_kraken"/*.kraken; do
    SAMPLE=$(basename "$KRAKEN_FILE" .kraken)
    R1="$REMOTE_DIR/${SAMPLE}_R1.fastq.gz"
    R2="$REMOTE_DIR/${SAMPLE}_R2.fastq.gz"
    ID_BASE="$LOCAL_OUTPUT_DIR/2_seqtk_filt/${SAMPLE}_tax${TAXID}_IDs"

    awk -v taxid="$TAXID" '$1 == "C" && $3 == taxid { split($2,a,"/"); split(a[1],b," "); print b[1] }' "$KRAKEN_FILE" | sort -u > "${ID_BASE}_all.txt"
    zcat "$R1" | awk 'NR%4==1 {gsub(/^@/, "", $1); split($1,a," "); print a[1]}' | sort -u > "${ID_BASE}_r1.txt"
    zcat "$R2" | awk 'NR%4==1 {gsub(/^@/, "", $1); split($1,a," "); print a[1]}' | sort -u > "${ID_BASE}_r2.txt"
    comm -12 "${ID_BASE}_all.txt" "${ID_BASE}_r1.txt" | comm -12 - "${ID_BASE}_r2.txt" > "${ID_BASE}_paired.txt"

    if [[ "$USE_SEQTK" == true ]]; then
      for END in 1 2; do
        zcat "$REMOTE_DIR/${SAMPLE}_R${END}.fastq.gz" | \
          seqtk subseq - "${ID_BASE}_paired.txt" | gzip > "$LOCAL_OUTPUT_DIR/2_seqtk_filt/${SAMPLE}_tax${TAXID}_R${END}.fastq.gz"
      done
    else
      awk '{ print "@" $1 }' "${ID_BASE}_paired.txt" > "${ID_BASE}_paired.grep"
      for END in 1 2; do
        zcat "$REMOTE_DIR/${SAMPLE}_R${END}.fastq.gz" | \
          paste - - - - | \
          awk -F'\t' '{ split($1, h, " "); print h[1] "\t" $0 }' | \
          grep -F -f "${ID_BASE}_paired.grep" | \
          cut -f2- | \
          tr '\t' '\n' | \
          gzip > "$LOCAL_OUTPUT_DIR/2_seqtk_filt/${SAMPLE}_tax${TAXID}_R${END}.fastq.gz"
      done
    fi
  done
  conda deactivate
fi

# === STEP 4: FASTQC AND FASTP ===
if should_run_step "$LOCAL_OUTPUT_DIR/3_fastqc" "FastQC and fastp"; then
  echo "[STEP 4] Running FastQC and fastp..."
  safe_activate qc_env
  mkdir -p "$LOCAL_OUTPUT_DIR/3_fastqc" "$LOCAL_OUTPUT_DIR/4_fastp"
  for R1 in "$LOCAL_OUTPUT_DIR"/2_seqtk_filt/*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    fastqc -t "$THREADS" -o "$LOCAL_OUTPUT_DIR/3_fastqc" "$R1" "$R2"
    fastp \
      -i "$R1" -I "$R2" \
      -o "$LOCAL_OUTPUT_DIR/4_fastp/${SAMPLE}_R1.fastq.gz" \
      -O "$LOCAL_OUTPUT_DIR/4_fastp/${SAMPLE}_R2.fastq.gz" \
      --detect_adapter_for_pe \
      --thread "$THREADS" \
      -j "$LOCAL_OUTPUT_DIR/4_fastp/${SAMPLE}_fastp.json" \
      -h "$LOCAL_OUTPUT_DIR/4_fastp/${SAMPLE}_fastp.html"
  done
  conda deactivate
fi

# === STEP 5: BWA + PILON ===
if should_run_step "$LOCAL_OUTPUT_DIR/5_bwamem_pol/pilon_polished.fasta" "BWA and Pilon"; then
  echo "[STEP 5] Aligning with BWA and polishing with Pilon..."
  safe_activate illuminareads_env
  mkdir -p "$LOCAL_OUTPUT_DIR/5_bwamem_pol"
  cd "$LOCAL_OUTPUT_DIR/5_bwamem_pol"

  [[ ! -f "${GENOME}.bwt" ]] && bwa index "$GENOME"

  BAM_LIST=()
  for R1 in "$LOCAL_OUTPUT_DIR"/4_fastp/*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    BAM="${SAMPLE}.bam"
    bwa mem -t "$THREADS" "$GENOME" "$R1" "$R2" | samtools sort -@ "$THREADS" -o "$BAM"
    samtools index "$BAM"
    BAM_LIST+=("$BAM")
  done

  samtools merge -@ "$THREADS" merged.bam "${BAM_LIST[@]}"
  samtools index merged.bam

  read -p "Enter the path to pilon.jar: " PILON_JAR
  [[ ! -f "$PILON_JAR" ]] && { echo "[ERROR] pilon.jar not found."; exit 1; }
  read -p "Enter amount of RAM for Java (e.g., 32G): " JAVA_RAM

  java -Xmx${JAVA_RAM} -jar "$PILON_JAR" \
    --genome "$GENOME" \
    --frags merged.bam \
    --output pilon_polished \
    --changes --vcf --mindepth 5 --fix all \
    > pilon.log 2>&1

  conda deactivate
fi

# === STEP 6: QUAST ===
if should_run_step "$LOCAL_OUTPUT_DIR/6_quast" "QUAST Analysis"; then
  echo "[STEP 6] Running QUAST..."
  mkdir -p "$LOCAL_OUTPUT_DIR/6_quast"
  quast.py -t "$THREADS" -o "$LOCAL_OUTPUT_DIR/6_quast" "$LOCAL_OUTPUT_DIR/5_bwamem_pol/pilon_polished.fasta"
fi

# === STEP 7: QUALIMAP ===
if should_run_step "$LOCAL_OUTPUT_DIR/7_mapping/qualimap_report" "Qualimap Analysis"; then
  echo "[STEP 7] Mapping and QC with Qualimap..."
  mkdir -p "$LOCAL_OUTPUT_DIR/7_mapping"
  read -p "Enter path to filtered FASTQ files (wildcard, e.g. /path/*.fastq.gz): " READ_GLOB
  READ_FILES=( $READ_GLOB )
  [[ ${#READ_FILES[@]} -eq 0 ]] && { echo "[ERROR] No FASTQ files match: $READ_GLOB"; exit 1; }

  MERGED_READS="$LOCAL_OUTPUT_DIR/7_mapping/merged_reads.fastq"
  echo "[INFO] Merging ${#READ_FILES[@]} FASTQ files..."
  zcat "${READ_FILES[@]}" > "$MERGED_READS"

  cd "$LOCAL_OUTPUT_DIR/7_mapping"
  safe_activate preprocessing_env
  minimap2 -ax map-ont "$LOCAL_OUTPUT_DIR/5_bwamem_pol/pilon_polished.fasta" "$MERGED_READS" -t "$THREADS" > mapped.sam
  conda deactivate

  safe_activate qc_env
  samtools view -bS mapped.sam | samtools sort -o mapped.sorted.bam
  samtools index mapped.sorted.bam
  qualimap bamqc \
    -bam mapped.sorted.bam \
    -outdir qualimap_report \
    -outformat PDF:HTML \
    --java-mem-size=16G
  conda deactivate
  rm mapped.sam
fi

cd "$ORIG_PWD"
echo "[DONE] Full Illumina read processing pipeline completed successfully."
