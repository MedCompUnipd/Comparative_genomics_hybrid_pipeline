#!/bin/bash
set -euo pipefail
shopt -s nullglob

# === CONFIG: NON-INTERACTIVE MODE ===
AUTO_MODE=false  # Set to true to skip interactive prompts for skipping steps

# === SAVE ORIGINAL DIRECTORY ===
ORIG_PWD="$(pwd)"
# Trap errors to return to original directory and exit cleanly on failure
trap 'cd "$ORIG_PWD" || exit; echo "[ERROR] Script interrotto o fallito. Uscita."; exit 1' ERR

# --- Conda setup ---
# Function to check if a command exists
command_exists () {
    command -v "$1" &> /dev/null
}

if ! command_exists conda; then
  echo "[ERROR] 'conda' non trovato. Assicurati che Conda sia installato e inizializzato."
  exit 1
fi
eval "$(conda shell.bash hook)"

# === SAFE CONDA ACTIVATION FUNCTION ===
safe_activate() {
  local env_name="$1"
  echo "[INFO] Attivazione ambiente Conda: $env_name"
  if ! conda activate "$env_name" &>/dev/null; then
    echo "[ERROR] Ambiente Conda '$env_name' non trovato o attivazione fallita. Controlla i tuoi ambienti Conda."
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

  # Determine step_dir: if output_check is a directory, use it directly. Otherwise, use its dirname.
  if [[ -d "$output_check" ]]; then
    step_dir="$output_check"
  else
    step_dir=$(dirname "$output_check")
  fi

  if [[ -e "$output_check" ]]; then
    echo "[INFO] L'output per '$step_name' esiste già in: $output_check"
    if [[ "$AUTO_MODE" == true ]]; then
      echo "[INFO] AUTO_MODE attivo: salto di $step_name"
      return 1 # Skip the step
    else
      read -r -p "Vuoi saltare questo passaggio? (s per saltare / n per sovrascrivere): " choice
      if [[ "$choice" =~ ^[Ss]$ ]]; then
        echo "[INFO] Salto di $step_name..."
        return 1 # Skip the step
      else
        echo "[INFO] Sovrascrittura di $step_name: rimozione di $step_dir..."
        rm -rf "$step_dir" # Remove the conflicting output directory
        mkdir -p "$step_dir" # Recreate it
        return 0
      fi
    fi
  else
    mkdir -p "$step_dir" # Ensure parent directory or output directory exists if output_check doesn't exist
    return 0
  fi
}


# === USER INPUT ===
read -p "Inserisci il percorso della cartella contenente le letture Illumina (FASTQ.gz): " ILLUMINA_READS_DIR # E.g., /path/to/illumina_raw
read -p "Inserisci il percorso della cartella contenente le letture ONT (FASTQ o FASTQ.gz): " ONT_READS_DIR # E.g., /path/to/ont_raw
read -p "Inserisci il percorso del file FASTA dell'assemblaggio de novo (es. dall'Script 3): " DENOVO_ASSEMBLY_FASTA
read -p "Inserisci il nome del progetto (es. hybrid_assembly_project): " PROJECT_NAME
read -p "Inserisci la directory di output radice desiderata per questo script: " OUTPUT_ROOT_DIR
read -p "Inserisci il percorso per creare/usare il database Kraken2: " DB_PATH # E.g., /path/to/kraken_db
read -p "Inserisci il TaxID da filtrare (es. 10359 per HCMV): " TAXID
read -p "Inserisci il modello Medaka (es. r1041_e82_400bps_sup_variant_v4.2.0): " MEDAKA_MODEL
read -p "Inserisci il numero di thread da utilizzare: " THREADS
read -p "Inserisci il percorso completo a pilon.jar (es. /path/to/pilon/pilon.jar): " PILON_JAR
read -p "Inserisci la quantità di RAM da allocare per Java (es. 32G): " JAVA_RAM

# --- Validate Inputs ---
if [[ ! -d "$ILLUMINA_READS_DIR" ]]; then
    echo "[ERROR] Directory letture Illumina non trovata: $ILLUMINA_READS_DIR. Uscita."
    exit 1
fi
if [[ ! -d "$ONT_READS_DIR" ]]; then
    echo "[ERROR] Directory letture ONT non trovata: $ONT_READS_DIR. Uscita."
    exit 1
fi
if [[ ! -f "$DENOVO_ASSEMBLY_FASTA" ]]; then
    echo "[ERROR] File FASTA assemblaggio de novo non trovato: $DENOVO_ASSEMBLY_FASTA. Uscita."
    exit 1
fi
if [[ ! -f "$PILON_JAR" ]]; then
    echo "[ERROR] File pilon.jar non trovato: $PILON_JAR. Uscita."
    exit 1
fi
mkdir -p "$OUTPUT_ROOT_DIR" || { echo "[ERROR] Impossibile creare directory di output radice: $OUTPUT_ROOT_DIR"; exit 1; }

# === CONDA ENVIRONMENT NAMES ===
KRAKEN_ENV="kraken2_env"
QC_ENV="qc_env"
ILLUMINA_ENV="illuminareads_env" # For BWA and Pilon
POLISHING_ENV="polishing_env" # For Medaka
PREPROCESSING_ENV="preprocessing_env" # For minimap2 for Qualimap


# === MAIN WORKFLOW ===
LOCAL_OUTPUT_DIR="${OUTPUT_ROOT_DIR}/${PROJECT_NAME}_hybrid_polishing"
mkdir -p "$LOCAL_OUTPUT_DIR" || { echo "[ERROR] Impossibile creare directory di output locale: $LOCAL_OUTPUT_DIR"; exit 1; }

echo ""
echo "--- Avvio Pipeline di Polishing Ibrido ---"

# === STEP 1: KRAKEN2 DB BUILD ===
if should_run_step "$DB_PATH/hash.k2d" "Kraken2 DB Build"; then
  echo "[STEP 1] Costruzione del database Kraken2..."
  safe_activate "$KRAKEN_ENV"
  kraken2-build --standard --db "$DB_PATH" --threads "$THREADS" || { echo "[ERROR] Kraken2 DB standard build fallito."; exit 1; }
  kraken2-build --build --db "$DB_PATH" --threads "$THREADS" || { echo "[ERROR] Kraken2 DB final build fallito."; exit 1; }
  conda deactivate
  echo "[INFO] Database Kraken2 costruito in: $DB_PATH"
else
  echo "[STEP 1] Salto della Costruzione del Database Kraken2."
fi

# === STEP 2: KRAKEN2 CLASSIFICATION (Illumina Reads) ===
KRAKEN_ILLUMINA_OUT_DIR="${LOCAL_OUTPUT_DIR}/1_kraken_illumina"
if should_run_step "$KRAKEN_ILLUMINA_OUT_DIR" "Kraken2 Classification (Illumina)"; then
  echo "[STEP 2] Classificazione delle letture Illumina con Kraken2..."
  safe_activate "$KRAKEN_ENV"
  # mkdir -p "$KRAKEN_ILLUMINA_OUT_DIR" # Handled by should_run_step

  shopt -s nullglob
  ILLUMINA_R1_FILES=("$ILLUMINA_READS_DIR"/*_R1.fastq.gz)
  shopt -u nullglob

  if [ ${#ILLUMINA_R1_FILES[@]} -eq 0 ]; then
      echo "[ERROR] Nessun file R1 FastQ.gz trovato in $ILLUMINA_READS_DIR. Interruzione della classificazione Kraken2 per Illumina."
      conda deactivate; exit 1
  fi

  for R1 in "${ILLUMINA_R1_FILES[@]}"; do
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    KRAKEN_OUTPUT_FILE="$KRAKEN_ILLUMINA_OUT_DIR/${SAMPLE}.kraken"
    KRAKEN_REPORT_FILE="$KRAKEN_ILLUMINA_OUT_DIR/${SAMPLE}.report"

    echo "→ Classificazione ${SAMPLE} con Kraken2..."
    kraken2 --db "$DB_PATH" --threads "$THREADS" --paired \
      --report "$KRAKEN_REPORT_FILE" \
      --output "$KRAKEN_OUTPUT_FILE" \
      "$R1" "$R2" || { echo "[ERROR] Kraken2 classification fallita per ${SAMPLE}."; conda deactivate; exit 1; }
  done
  echo "[INFO] Classificazione Kraken2 completata. Report in $KRAKEN_ILLUMINA_OUT_DIR"
  conda deactivate
else
  echo "[STEP 2] Salto della Classificazione Kraken2 (Illumina)."
fi

# === STEP 3: FILTERING READS by TaxID (Illumina) ===
FILTERED_ILLUMINA_READS_DIR="${LOCAL_OUTPUT_DIR}/2_illumina_filtered"
if should_run_step "$FILTERED_ILLUMINA_READS_DIR" "Read Filtering (Illumina by TaxID)"; then
  echo "[STEP 3] Filtraggio delle letture Illumina per TaxID $TAXID..."
  safe_activate "$KRAKEN_ENV"
  # mkdir -p "$FILTERED_ILLUMINA_READS_DIR" # Handled by should_run_step

  USE_SEQTK=false
  if command_exists seqtk; then
    USE_SEQTK=true
    echo "[INFO] Utilizzo di seqtk per il filtraggio."
  else
    echo "[WARNING] seqtk non trovato. Ripiego sul filtraggio con awk."
  fi

  shopt -s nullglob
  KRAKEN_FILES=("$KRAKEN_ILLUMINA_OUT_DIR"/*.kraken)
  shopt -u nullglob

  if [ ${#KRAKEN_FILES[@]} -eq 0 ]; then
      echo "[ERROR] Nessun file Kraken2 .kraken trovato in $KRAKEN_ILLUMINA_OUT_DIR. Impossibile filtrare."
      conda deactivate; exit 1
  fi

  for KRAKEN_FILE in "${KRAKEN_FILES[@]}"; do
    SAMPLE=$(basename "$KRAKEN_FILE" .kraken)
    R1_RAW="$ILLUMINA_READS_DIR/${SAMPLE}_R1.fastq.gz"
    R2_RAW="$ILLUMINA_READS_DIR/${SAMPLE}_R2.fastq.gz"
    ID_BASE="$FILTERED_ILLUMINA_READS_DIR/${SAMPLE}_tax${TAXID}_IDs"

    # Extract read IDs for the target TaxID
    awk -v taxid="$TAXID" '$1 == "C" && $3 == taxid { split($2,a,"/"); split(a[1],b," "); print b[1] }' "$KRAKEN_FILE" | sort -u > "${ID_BASE}_all.txt"
    zcat "$R1_RAW" | awk 'NR%4==1 {gsub(/^@/, "", $1); split($1,a," "); print a[1]}' | sort -u > "${ID_BASE}_r1.txt"
    zcat "$R2_RAW" | awk 'NR%4==1 {gsub(/^@/, "", $1); split($1,a," "); print a[1]}' | sort -u > "${ID_BASE}_r2.txt"
    # Find common IDs present in both R1, R2, and classified as target TaxID
    comm -12 "${ID_BASE}_all.txt" "${ID_BASE}_r1.txt" | comm -12 - "${ID_BASE}_r2.txt" > "${ID_BASE}_paired.txt"

    if [[ "$USE_SEQTK" == true ]]; then
      for END in 1 2; do
        zcat "$ILLUMINA_READS_DIR/${SAMPLE}_R${END}.fastq.gz" | \
          seqtk subseq - "${ID_BASE}_paired.txt" | gzip > "$FILTERED_ILLUMINA_READS_DIR/${SAMPLE}_tax${TAXID}_R${END}.fastq.gz" || { echo "[ERROR] seqtk subseq fallito per ${SAMPLE}_R${END}."; conda deactivate; exit 1; }
      done
    else
      awk '{ print "@" $1 }' "${ID_BASE}_paired.txt" > "${ID_BASE}_paired.grep"
      for END in 1 2; do
        zcat "$ILLUMINA_READS_DIR/${SAMPLE}_R${END}.fastq.gz" | \
          paste - - - - | \
          awk -F'\t' '{ split($1, h, " "); print h[1] "\t" $0 }' | \
          grep -F -f "${ID_BASE}_paired.grep" | \
          cut -f2- | \
          tr '\t' '\n' | \
          gzip > "$FILTERED_ILLUMINA_READS_DIR/${SAMPLE}_tax${TAXID}_R${END}.fastq.gz" || { echo "[ERROR] awk filtering fallito per ${SAMPLE}_R${END}."; conda deactivate; exit 1; }
      done
    fi
    # Clean up temporary ID files
    rm -f "${ID_BASE}_all.txt" "${ID_BASE}_r1.txt" "${ID_BASE}_r2.txt" "${ID_BASE}_paired.txt" "${ID_BASE}_paired.grep"
  done
  echo "[INFO] Filtraggio per TaxID completato. Letture filtrate in $FILTERED_ILLUMINA_READS_DIR"
  conda deactivate
else
  echo "[STEP 3] Salto del Filtraggio delle Letture (Illumina by TaxID)."
fi

# === STEP 4: FASTQC AND FASTP (on filtered Illumina reads) ===
FASTQC_OUT_DIR="${LOCAL_OUTPUT_DIR}/3_fastqc"
FASTP_OUT_DIR="${LOCAL_OUTPUT_DIR}/4_fastp"
if should_run_step "$FASTQC_OUT_DIR" "FastQC and fastp"; then # Using fastqc dir as check for the pair of tools
  echo "[STEP 4] Esecuzione di FastQC e fastp sulle letture Illumina filtrate..."
  safe_activate "$QC_ENV" # fastqc and fastp are in qc_env
  # mkdir -p "$FASTQC_OUT_DIR" "$FASTP_OUT_DIR" # Handled by should_run_step for FASTQC_OUT_DIR

  shopt -s nullglob
  FILTERED_R1_FILES=("$FILTERED_ILLUMINA_READS_DIR"/*_R1.fastq.gz)
  shopt -u nullglob

  if [ ${#FILTERED_R1_FILES[@]} -eq 0 ]; then
      echo "[ERROR] Nessun file FastQ.gz filtrato R1 trovato in $FILTERED_ILLUMINA_READS_DIR. Interruzione di FastQC/fastp."
      conda deactivate; exit 1
  fi

  for R1 in "${FILTERED_R1_FILES[@]}"; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    
    echo "→ Running FastQC for ${SAMPLE}..."
    fastqc -t "$THREADS" -o "$FASTQC_OUT_DIR" "$R1" "$R2" || { echo "[ERROR] FastQC fallito per ${SAMPLE}."; conda deactivate; exit 1; }

    echo "→ Running fastp for ${SAMPLE}..."
    fastp \
      -i "$R1" -I "$R2" \
      -o "$FASTP_OUT_DIR/${SAMPLE}_fastp_R1.fastq.gz" \
      -O "$FASTP_OUT_DIR/${SAMPLE}_fastp_R2.fastq.gz" \
      --detect_adapter_for_pe \
      --thread "$THREADS" \
      -j "$FASTP_OUT_DIR/${SAMPLE}_fastp.json" \
      -h "$FASTP_OUT_DIR/${SAMPLE}_fastp.html" || { echo "[ERROR] fastp fallito per ${SAMPLE}."; conda deactivate; exit 1; }
  done
  echo "[INFO] FastQC e fastp completati. Report in $FASTQC_OUT_DIR e letture trimmed in $FASTP_OUT_DIR"
  conda deactivate
else
  echo "[STEP 4] Salto di FastQC e fastp."
fi


# === STEP 5: PILON POLISHING (with Fastp-processed Illumina reads) ===
PILON_POLISHED_FASTA="${LOCAL_OUTPUT_DIR}/5_pilon_pol/pilon_polished.fasta"
if should_run_step "$PILON_POLISHED_FASTA" "Pilon Polishing"; then
    echo "[STEP 5] Esecuzione del polishing con Pilon utilizzando le letture Illumina processate..."
    # mkdir -p "$(dirname "$PILON_POLISHED_FASTA")" # Handled by should_run_step

    safe_activate "$ILLUMINA_ENV" # BWA and Pilon are in illuminareads_env

    # 1. Map Fastp-processed Illumina reads to the de novo assembly using BWA-MEM
    echo "[INFO] Mappatura delle letture Illumina processate sull'assemblaggio de novo con BWA-MEM..."
    
    # Ensure genome is indexed for BWA
    if [[ ! -f "${DENOVO_ASSEMBLY_FASTA}.bwt" ]]; then
        echo "[INFO] Indicizzazione del genoma de novo per BWA..."
        bwa index "$DENOVO_ASSEMBLY_FASTA" || { echo "[ERROR] Indicizzazione BWA fallita per l'assemblaggio de novo."; exit 1; }
    fi

    # Process each fastp-processed sample and merge BAMs
    BAM_LIST=()
    shopt -s nullglob
    FASTP_R1_FILES=("$FASTP_OUT_DIR"/*_fastp_R1.fastq.gz)
    shopt -u nullglob

    if [ ${#FASTP_R1_FILES[@]} -eq 0 ]; then
        echo "[ERROR] Nessun file FastQ.gz processato da fastp R1 trovato in $FASTP_OUT_DIR. Interruzione del polishing Pilon."
        conda deactivate; exit 1
    fi

    for R1 in "${FASTP_R1_FILES[@]}"; do
        SAMPLE=$(basename "$R1" _fastp_R1.fastq.gz)
        R2="${R1/_fastp_R1.fastq.gz/_fastp_R2.fastq.gz}"
        BAM="${LOCAL_OUTPUT_DIR}/5_pilon_pol/${SAMPLE}.bam"

        echo "→ Allineamento ${SAMPLE} con BWA-MEM..."
        bwa mem -t "$THREADS" "$DENOVO_ASSEMBLY_FASTA" "$R1" "$R2" | samtools sort -@ "$THREADS" -o "$BAM" - || { echo "[ERROR] Mappatura BWA/Samtools sort fallita per ${SAMPLE}."; conda deactivate; exit 1; }
        samtools index "$BAM" || { echo "[ERROR] Indicizzazione Samtools fallita per ${SAMPLE}."; conda deactivate; exit 1; }
        BAM_LIST+=("$BAM")
    done

    local MERGED_BAM="${LOCAL_OUTPUT_DIR}/5_pilon_pol/illumina_mapped_merged.bam"
    echo "[INFO] Unione dei file BAM individuali in $MERGED_BAM..."
    samtools merge -@ "$THREADS" "$MERGED_BAM" "${BAM_LIST[@]}" || { echo "[ERROR] Unione dei file BAM fallita."; conda deactivate; exit 1; }
    samtools index "$MERGED_BAM" || { echo "[ERROR] Indicizzazione del BAM unito fallita."; conda deactivate; exit 1; }
    
    echo "[INFO] Letture Illumina allineate. File BAM: $MERGED_BAM"

    # 2. Run Pilon
    echo "[INFO] Esecuzione di Pilon..."
    # Pilon often requires Java. Assumed openjdk=8 is in illuminareads_env
    java -Xmx${JAVA_RAM} -jar "$PILON_JAR" --genome "$DENOVO_ASSEMBLY_FASTA" --frags "$MERGED_BAM" --output pilon_polished --outdir "$(dirname "$PILON_POLISHED_FASTA")" --threads "$THREADS" --changes --vcf --mindepth 5 --fix all > "${LOCAL_OUTPUT_DIR}/5_pilon_pol/pilon.log" 2>&1 || { echo "[ERROR] Polishing Pilon fallito. Controlla il log: ${LOCAL_OUTPUT_DIR}/5_pilon_pol/pilon.log"; exit 1; }
    
    echo "[INFO] Polishing Pilon completato. Output: $PILON_POLISHED_FASTA"
    conda deactivate
else
    echo "[STEP 5] Salto del Polishing Pilon."
fi


# === STEP 6: MEDAKA POLISHING (with ONT reads) ===
MEDAKA_POLISHED_FASTA="${LOCAL_OUTPUT_DIR}/6_medaka_pol/medaka_polished.fasta"
if should_run_step "$MEDAKA_POLISHED_FASTA" "Medaka Polishing"; then
    echo "[STEP 6] Esecuzione del polishing Medaka con letture ONT..."
    # mkdir -p "$(dirname "$MEDAKA_POLISHED_FASTA")" # Handled by should_run_step

    safe_activate "$POLISHING_ENV"

    # Merge ONT reads if not a single file (similar to 2_tax_analysis.sh)
    ONT_MERGED_FASTQ="${LOCAL_OUTPUT_DIR}/6_medaka_pol/${PROJECT_NAME}_ont_merged.fastq"
    
    # Check for existing merged ONT file first
    if [ ! -f "$ONT_MERGED_FASTQ" ] || [ ! -s "$ONT_MERGED_FASTQ" ]; then
        echo "[INFO] Nessun file ONT unito trovato. Unione di tutti i file .fastq e .fastq.gz in: $ONT_READS_DIR"
        find "$ONT_READS_DIR" -type f -name "*.fastq" -print0 | xargs -0 cat > "$ONT_MERGED_FASTQ" || \
        find "$ONT_READS_DIR" -type f -name "*.fastq.gz" -print0 | xargs -0 gzip -cd >> "$ONT_MERGED_FASTQ" || { echo "[ERROR] Impossibile unire le letture ONT."; conda deactivate; exit 1; }
    fi

    if [ ! -s "$ONT_MERGED_FASTQ" ]; then
        echo "[ERROR] Nessun file FastQ ONT trovato o il file unito è vuoto. Interruzione di Medaka."
        conda deactivate; exit 1
    fi

    medaka_consensus -i "$ONT_MERGED_FASTQ" -d "$PILON_POLISHED_FASTA" -o "$(dirname "$MEDAKA_POLISHED_FASTA")" -t "$THREADS" -m "$MEDAKA_MODEL" || { echo "[ERROR] Polishing Medaka fallito."; conda deactivate; exit 1; }
    # Medaka output is consensus.fasta within the output directory, rename it
    mv "$(dirname "$MEDAKA_POLISHED_FASTA")/consensus.fasta" "$MEDAKA_POLISHED_FASTA" || { echo "[ERROR] Impossibile spostare il file consensus.fasta di Medaka."; conda deactivate; exit 1; }
    
    echo "[INFO] Polishing Medaka completato. Output: $MEDAKA_POLISHED_FASTA"
    conda deactivate
else
    echo "[STEP 6] Salto del Polishing Medaka."
fi

# === STEP 7: QUAST QC for Hybrid Assembly ===
HYBRID_QUAST_REPORT_DIR="${LOCAL_OUTPUT_DIR}/7_quast"
if should_run_step "$HYBRID_QUAST_REPORT_DIR/report.html" "Hybrid Assembly QUAST QC"; then
  echo "[STEP 7] Esecuzione di QUAST per l'assemblaggio ibrido..."
  safe_activate "$QC_ENV"
  quast.py -t "$THREADS" -o "$HYBRID_QUAST_REPORT_DIR" "$MEDAKA_POLISHED_FASTA" || { echo "[ERROR] QUAST fallito per l'assemblaggio ibrido."; conda deactivate; exit 1; }
  echo "[INFO] QC QUAST per l'assemblaggio ibrido completato. Report in $HYBRID_QUAST_REPORT_DIR"
  conda deactivate
else
  echo "[STEP 7] Salto del QC QUAST per l'assemblaggio ibrido."
fi

# === STEP 8: Mapping with minimap2 and Qualimap (for the final hybrid assembly) ===
HYBRID_MAPPING_DIR="${LOCAL_OUTPUT_DIR}/8_hybrid_mapping"
if should_run_step "$HYBRID_MAPPING_DIR/qualimap_report/genome_results.txt" "Hybrid Mapping and Qualimap"; then
  echo "[STEP 8] Mappatura delle letture ONT sull'assemblaggio ibrido con minimap2 ed esecuzione di Qualimap..."
  # mkdir -p "$HYBRID_MAPPING_DIR" # Handled by should_run_step

  safe_activate "$PREPROCESSING_ENV" # minimap2 and samtools are in preprocessing_env
  ONT_MERGED_FASTQ="${LOCAL_OUTPUT_DIR}/6_medaka_pol/${PROJECT_NAME}_ont_merged.fastq" # Re-use merged ONT reads from Medaka step

  if [ ! -f "$ONT_MERGED_FASTQ" ]; then
      echo "[ERROR] File FASTQ ONT unito non trovato: $ONT_MERGED_FASTQ. Interruzione di mappatura/Qualimap."
      conda deactivate; exit 1
  fi

  local HYBRID_BAM="${HYBRID_MAPPING_DIR}/${PROJECT_NAME}_hybrid_mapped.bam"
  echo "→ Mappatura con minimap2..."
  minimap2 -ax map-ont "$MEDAKA_POLISHED_FASTA" "$ONT_MERGED_FASTQ" | samtools view -bS - | samtools sort -o "$HYBRID_BAM" - || { echo "[ERROR] Mappatura minimap2/Samtools fallita per l'assemblaggio ibrido."; conda deactivate; exit 1; }
  samtools index "$HYBRID_BAM" || { echo "[ERROR] Indicizzazione Samtools fallita per il BAM ibrido."; conda deactivate; exit 1; }
  echo "[INFO] Letture ONT mappate sull'assemblaggio ibrido. File BAM: $HYBRID_BAM"
  conda deactivate

  safe_activate "$QC_ENV" # Qualimap is in qc_env
  echo "→ Esecuzione di Qualimap..."
  qualimap bamqc -bam "$HYBRID_BAM" -outdir "${HYBRID_MAPPING_DIR}/qualimap_report" --java-mem-size=${JAVA_RAM} || { echo "[ERROR] Qualimap fallito per la mappatura dell'assemblaggio ibrido."; conda deactivate; exit 1; }
  echo "[INFO] QC Qualimap per l'assemblaggio ibrido completato. Report in ${HYBRID_MAPPING_DIR}/qualimap_report"
  conda deactivate
else
  echo "[STEP 8] Salto della Mappatura Ibrida e Qualimap."
fi

echo ""
echo "-------------------------------------"
echo "Pipeline di polishing ibrido completata con successo!"
echo "Assemblaggio finale lucidato: $MEDAKA_POLISHED_FASTA"
echo "Report QUAST: $HYBRID_QUAST_REPORT_DIR"
echo "Report Qualimap: ${HYBRID_MAPPING_DIR}/qualimap_report"

