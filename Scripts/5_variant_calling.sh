#!/bin/bash
set -euo pipefail
shopt -s nullglob

# === LOAD CONDA (portable) ===
# Verify Conda installation and initialize it.
if ! command -v conda &> /dev/null; then
  echo "[ERROR] Conda non trovato. Assicurati che conda sia installato e inizializzato."
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
    echo "[INFO] L'output per '$step_name' esiste già in: $output_check"
    read -r -p "Vuoi saltare questo passaggio? (s per saltare / n per sovrascrivere): " choice
    if [[ "$choice" =~ ^[Ss]$ ]]; then
      echo "[INFO] Salto di $step_name..."
      return 1 # Skip the step
    else
      echo "[INFO] Sovrascrittura di $step_name: rimozione di $step_dir..."
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
read -r -p "Inserisci il percorso del file FASTA del genoma di riferimento: " REFERENCE
read -r -p "Inserisci il percorso del file FASTA del genoma de novo (il tuo assemblaggio pulito): " DENOVO_GENOME
read -r -p "Inserisci il percorso della directory di output desiderata: " OUTDIR
read -r -p "Inserisci il numero di thread da utilizzare: " THREADS

# ---
# Conditional annotation input
PERFORM_ANNOTATION="n" # Default to no annotation
read -r -p "Vuoi eseguire l'annotazione delle varianti? (s/n): " ANNOTATION_CHOICE
if [[ "$ANNOTATION_CHOICE" =~ ^[Ss]$ ]]; then
  PERFORM_ANNOTATION="y"
  read -r -p "Inserisci il percorso del file GFF3 di annotazione (verrà usato da bcftools csq): " GFF
  echo "[SUGGERIMENTO] Per il file GFF3, puoi cercarlo nel database EMBL-EBI tramite ENA (https://www.ebi.ac.uk/ena/). Cerca il numero di accesso del tuo genoma di riferimento e scarica il file in formato GFF3."
else
  echo "[INFO] L'annotazione delle varianti verrà saltata."
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
  [[ -f "$FILE" ]] || { echo "[ERROR] Input mancante: $FILE"; exit 1; }
done
# Create the main output directory.
mkdir -p "$OUTDIR"
echo "[INFO] Avvio della pipeline di confronto tra il tuo assemblaggio de novo pulito ($DENOVO_GENOME) e il riferimento ($REFERENCE)."

# ---

# === STEP 0: MUMmer Comparison ===
# This step compares the large-scale structure of your polished de novo genome with the reference genome.
# Useful for identifying large rearrangements, inversions, or translocations.
MUMMER_OUTPUT_CHECK="$OUTDIR/aln_mummer/ref_vs_denovo.coords"
if should_run_step "$MUMMER_OUTPUT_CHECK" "MUMmer Comparison"; then
  echo "[0] Esecuzione del confronto MUMmer (riferimento vs tuo assemblaggio de novo)..."
  # mkdir -p "$OUTDIR/aln_mummer" # Handled by should_run_step
  conda activate mummer_env || { echo "[ERROR] Impossibile attivare 'mummer_env'. Controlla i tuoi ambienti Conda."; exit 1; }
  nucmer --prefix="$OUTDIR/aln_mummer/ref_vs_denovo" "$REFERENCE" "$DENOVO_GENOME" || { echo "[ERROR] nucmer fallito."; exit 1; }
  delta-filter -1 "$OUTDIR/aln_mummer/ref_vs_denovo.delta" > "$OUTDIR/aln_mummer/ref_vs_denovo.filtered.delta" || { echo "[ERROR] delta-filter fallito."; exit 1; }
  show-coords -rcl "$OUTDIR/aln_mummer/ref_vs_denovo.filtered.delta" > "$OUTDIR/aln_mummer/ref_vs_denovo.coords" || { echo "[ERROR] show-coords fallito."; exit 1; }
  
  # Check if gnuplot is available before trying mummerplot
  if ! command -v gnuplot &> /dev/null; then
      echo "[WARNING] 'gnuplot' non trovato. I grafici MUMmer non possono essere generati. Installa gnuplot manualmente (es. sudo apt install gnuplot)."
  else
      mummerplot --png --layout --prefix "$OUTDIR/aln_mummer/ref_vs_denovo" "$OUTDIR/aln_mummer/ref_vs_denovo.filtered.delta" || { echo "[ERROR] mummerplot fallito."; } # mummerplot often fails gracefully
  fi
  
  conda deactivate
  echo "[INFO] I risultati del confronto MUMmer si trovano in $OUTDIR/aln_mummer/"
else
  echo "[0] Salto del confronto MUMmer."
fi

# ---

# === STEP 1: Variant Calling (De novo assembly vs Reference) using BCFtools ===
# This step identifies single nucleotide variants (SNPs) and small Indels between your polished de novo genome
# and the reference genome using BCFtools.
BCFTOOLS_VCF_OUTPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref.vcf.gz"
if should_run_step "$BCFTOOLS_VCF_OUTPUT" "BCFtools-based Variant Calling (De novo vs Ref)"; then
  echo "[1] Chiamata delle varianti dal tuo assemblaggio de novo vs riferimento usando BCFtools..."
  # mkdir -p "$OUTDIR/variant_calling_bcftools" # Handled by should_run_step

  conda activate preprocessing_env || { echo "[ERROR] Impossibile attivare 'preprocessing_env'. Controlla i tuoi ambienti Conda."; exit 1; }
  minimap2 -ax asm5 "$REFERENCE" "$DENOVO_GENOME" > "$OUTDIR/variant_calling_bcftools/aln.sam" || { echo "[ERROR] minimap2 fallito."; exit 1; }
  conda deactivate

  conda activate qc_env || { echo "[ERROR] Impossibile attivare 'qc_env'. Controlla i tuoi ambienti Conda."; exit 1; }
  samtools sort -@ "$THREADS" -o "$OUTDIR/variant_calling_bcftools/aln.sorted.bam" "$OUTDIR/variant_calling_bcftools/aln.sam" || { echo "[ERROR] samtools sort fallito."; exit 1; }
  samtools index "$OUTDIR/variant_calling_bcftools/aln.sorted.bam" || { echo "[ERROR] samtools index fallito."; exit 1; }
  rm "$OUTDIR/variant_calling_bcftools/aln.sam"

  bcftools mpileup -Ou -f "$REFERENCE" "$OUTDIR/variant_calling_bcftools/aln.sorted.bam" | \
    bcftools call -mv -Oz -o "$BCFTOOLS_VCF_OUTPUT" || { echo "[ERROR] bcftools mpileup/call fallito."; exit 1; }
  bcftools index "$BCFTOOLS_VCF_OUTPUT" || { echo "[ERROR] bcftools index fallito."; exit 1; }
  conda deactivate
  echo "[INFO] Le varianti BCFtools-based de novo vs Reference si trovano in $BCFTOOLS_VCF_OUTPUT"
else
  echo "[1] Salto della chiamata delle varianti BCFtools-based (De novo vs Ref)."
fi

# ---

# === STEP 2: Filter BCFtools Variants for High Confidence (QUAL >= 30) ===
# Filter the variants called by BCFtools (de novo vs reference) to keep only high-confidence calls.
FILTERED_VCF_OUTPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref_filtered.vcf.gz"
if should_run_step "$FILTERED_VCF_OUTPUT" "Filter BCFtools-based Variants"; then
  echo "[2] Filtro delle varianti BCFtools-based de novo vs riferimento (QUAL >= 30) per alta confidenza..."
  # mkdir -p "$OUTDIR/variant_calling_bcftools" # Handled by should_run_step
  VCF_INPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref.vcf.gz"
  VCF_TEMP_OUTPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref_filtered.vcf" # Non-gzipped temporary file
  conda activate qc_env || { echo "[ERROR] Impossibile attivare 'qc_env'. Controlla i tuoi ambienti Conda."; exit 1; } # Assuming bcftools is here
  bcftools filter -e 'QUAL<30' "$VCF_INPUT" -o "$VCF_TEMP_OUTPUT" -Ov || { echo "[ERROR] bcftools filter fallito."; exit 1; }
  bgzip -c "$VCF_TEMP_OUTPUT" > "$FILTERED_VCF_OUTPUT" || { echo "[ERROR] bgzip fallito."; exit 1; }
  bcftools index "$FILTERED_VCF_OUTPUT" || { echo "[ERROR] bcftools index fallito."; exit 1; }
  rm "$VCF_TEMP_OUTPUT" # Clean up temporary non-gzipped file
  conda deactivate
  echo "[INFO] Varianti filtrate ad alta confidenza de novo vs Reference: $FILTERED_VCF_OUTPUT"
else
  echo "[2] Salto del filtro delle varianti BCFtools-based."
fi

# ---

# === STEP 3: Annotate Filtered Variants using BCFtools csq (Conditional) ===
ANNOTATED_VCF_OUTPUT="$OUTDIR/variant_calling_bcftools/denovo_vs_ref_filtered_annotated_bcftools.vcf"
if [[ "$PERFORM_ANNOTATION" == "y" ]]; then
  if should_run_step "$ANNOTATED_VCF_OUTPUT" "Annotation of Filtered De novo vs Reference Variants (BCFtools csq)"; then
    echo "[3] Annotazione delle varianti filtrate de novo vs riferimento usando BCFtools csq..."
    conda activate qc_env || { echo "[ERROR] Impossibile attivare 'qc_env'. Controlla i tuoi ambienti Conda."; exit 1; }

    # Ensure the reference FASTA is indexed for bcftools csq
    if [[ ! -f "$REFERENCE.fai" ]]; then
      echo "[INFO] Creazione dell'indice FASTA per il genoma di riferimento ($REFERENCE.fai)..."
      samtools faidx "$REFERENCE" || { echo "[ERROR] samtools faidx fallito."; exit 1; }
    fi

    # Annotate with bcftools csq
    bcftools csq -f "$REFERENCE" -g "$GFF" -o "$ANNOTATED_VCF_OUTPUT" "$FILTERED_VCF_OUTPUT" || { echo "[ERROR] bcftools csq annotation fallito."; exit 1; }
    conda deactivate
    echo "[INFO] Varianti annotate filtrate de novo vs Reference (BCFtools csq): $ANNOTATED_VCF_OUTPUT"
  else
    echo "[3] Salto dell'annotazione delle varianti filtrate (BCFtools csq)."
  fi
else
  echo "[INFO] L'annotazione delle varianti è stata saltata come richiesto."
fi

# ---

echo ""
echo "[FATTO] Pipeline di confronto assemblaggio de novo vs riferimento completata con successo."
echo "--- Output chiave del confronto ---"
echo "-> Differenze strutturali (grafici MUMmer e coordinate): $OUTDIR/aln_mummer/"
if [[ "$PERFORM_ANNOTATION" == "y" ]]; then
  echo "-> Mutazioni puntiformi ad alta confidenza (de novo vs Reference, basate su BCFtools, annotate con BCFtools csq): $ANNOTATED_VCF_OUTPUT"
else
  echo "-> Mutazioni puntiformi ad alta confidenza (de novo vs Reference, basate su BCFtools, NON ANNOTATE): $FILTERED_VCF_OUTPUT"
fi
echo ""
echo "Questo script fornisce un confronto robusto del tuo assemblaggio de novo pulito rispetto al genoma di riferimento."
