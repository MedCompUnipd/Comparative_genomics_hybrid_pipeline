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

echo "--- Assembly and Polishing Script ---"
echo "-------------------------------------"

# === USER INPUT ===
read -p "==> Input reads (FASTQ file or directory): " READS
read -p "==> Root output directory for assembly steps: " ROOT
read -p "==> Estimated genome size (e.g. 235k): " GENOME_SIZE
read -p "==> Medaka model (e.g. r1041_e82_400bps_sup_variant_v4.2.0): " MODEL
read -p "==> Threads to use: " THREADS
read -p "==> Project name (e.g. HCMV): " PROJECT_NAME

# === ENVIRONMENT NAMES ===
PREPROCESSING_ENV="preprocessing_env"
ASSEMBLY_ENV="assembly_env"
POLISHING_ENV="polishing_env"
QC_ENV="qc_env"
NDN_ENV="ndn_env"


declare -A assemblies
mkdir -p "$ROOT"

# Control input reads
if [[ ! -e "$READS" ]]; then
  echo "[ERROR] Input reads not found: $READS"
  exit 1
fi

if [[ -d "$READS" ]]; then
  MERGED_FILE=$(find "$READS" -type f \( -iname "*merged.fastq.gz" -o -iname "*merged.fastq" -o -iname "*merged.fq.gz" -o -iname "*merged.fq" \) | head -n 1)
elif [[ -f "$READS" ]]; then
  MERGED_FILE="$READS"
else
  echo "[ERROR] Reads input is neither a file nor a directory: $READS"
  exit 1
fi

if [[ ! -f "$MERGED_FILE" ]]; then
  echo "[ERROR] No merged FASTQ file found in input reads path."
  exit 1
fi

echo "[INFO] Using merged FASTQ: $MERGED_FILE"


# ---------NEXTDENOVO--------------
read -p "==> Run NextDenovo? [y/n] " run_nextdenovo
if [[ "$run_nextdenovo" == "y" ]]; then
  echo "[INFO] Running NextDenovo..."
  conda deactivate || true
  conda activate "$NDN_ENV"

  ND_DIR="$ROOT/nextdenovo"
  mkdir -p "$ND_DIR"

  ND_CFG="$ND_DIR/run.cfg"
  ND_OUTPUT=$(find "$ND_DIR" -name "*nd.asm*.fasta" | head -n 1)

  if [[ ! -s "$ND_OUTPUT" ]]; then
    if [[ -d "$READS" ]]; then
      find "$READS" -type f \( -iname "*.fastq*" -o -iname "*.fq*" \) > "$ND_DIR/input.fofn"
    else
      echo "$READS" > "$ND_DIR/input.fofn"
    fi

    cat > "$ND_CFG" <<EOF
[General]
job_type = local
job_prefix = $PROJECT_NAME
task = all
rewrite = yes
dataroot = $ROOT
workdir = $ND_DIR
input_fofn = input.fofn
read_type = ont

[correct_option]
read_cutoff = 1k
genome_size = $GENOME_SIZE
seed_cutoff = 10k
sort_options = -m 16g
correction_options = -p $THREADS
minimap2_options = -t $THREADS

[assemble_option]
minimap2_options = -t $THREADS
EOF

    nextDenovo "$ND_CFG" | tee "$ND_DIR/nextdenovo.log"
    ND_OUTPUT=$(find "$ND_DIR" -name "*nd.asm*.fasta" | head -n 1)
  fi

  conda deactivate

  if [[ -s "$ND_OUTPUT" ]]; then
    assemblies[nextdenovo]="$ND_OUTPUT"
  else
    echo "[WARNING] NextDenovo output not found or empty."
  fi
else
  echo "[INFO] Skipping NextDenovo."
fi




# === CANU ===
read -p "==> Do you want to run Canu? [y/n] " run_canu
Canu_DIR="$ROOT/canu"
Canu_OUT="$Canu_DIR/${PROJECT_NAME}.contigs.fasta"
conda deactivate || true
conda activate assembly_env
if [[ "$run_canu" == "y" ]]; then
  if [[ ! -s "$Canu_OUT" ]]; then
    echo "[INFO] Running Canu..."
    mkdir -p "$Canu_DIR"
    canu -p "$PROJECT_NAME" -d "$Canu_DIR" genomeSize="$GENOME_SIZE" \
         -nanopore-raw "$MERGED_FILE" useGrid=false minReadLength=500
  fi
  assemblies[canu]="$Canu_OUT"
else
  echo "[INFO] Skipping Canu."
fi
conda deactivate

# === FLYE ===
read -p "==> Do you want to run Flye? [y/n] " run_flye
conda deactivate || true
conda activate "$ASSEMBLY_ENV"
if [[ "$run_flye" == "y" ]]; then
  Flye_DIR="$ROOT/flye"
  Flye_OUT="$Flye_DIR/assembly.fasta"
  if [[ ! -s "$Flye_OUT" ]]; then
    echo "[INFO] Running Flye..."
    mkdir -p "$Flye_DIR"
    flye --nano-raw "$MERGED_FILE" --out-dir "$Flye_DIR" \
         --threads "$THREADS" --genome-size "$GENOME_SIZE" \
         --min-overlap 1000 --asm-coverage 50 | tee "$Flye_DIR/flye.log"
  fi
  assemblies[flye]="$Flye_OUT"
else
  echo "[INFO] Skipping Flye."
fi

# === WTDBG2 ===
read -p "==> Do you want to run Wtdbg2? [y/n] " run_wtdbg2
if [[ "$run_wtdbg2" == "y" ]]; then
  if [[ ! -f "$Canu_DIR/${PROJECT_NAME}.trimmedReads.fasta.gz" ]]; then
    echo "[ERROR] trimmedReads file not found in $Canu_DIR"
    exit 1
  fi
  Wtdbg2_DIR="$ROOT/wtdbg2"
  Wtdbg2_OUT="$Wtdbg2_DIR/wtdbg_assembly.fasta"
  if [[ ! -s "$Wtdbg2_OUT" ]]; then
    echo "[INFO] Running Wtdbg2..."
    mkdir -p "$Wtdbg2_DIR"
    gunzip -c "$Canu_DIR/${PROJECT_NAME}.trimmedReads.fasta.gz" > "$Wtdbg2_DIR/corrected_reads.fasta"
    wtdbg2 -i "$Wtdbg2_DIR/corrected_reads.fasta" -o "$Wtdbg2_DIR/wtdbg_assembly" \
           -t "$THREADS" -g "$GENOME_SIZE" -p 0 -k 15 -AS 2 -s 0.05 -L 5000
    wtpoa-cns -t "$THREADS" -i "$Wtdbg2_DIR/wtdbg_assembly.ctg.lay.gz" -fo "$Wtdbg2_OUT"
  fi
  assemblies[wtdbg2]="$Wtdbg2_OUT"
else
  echo "[INFO] Skipping Wtdbg2."
fi

# === HIFIASM ===
read -p "==> Do you want to run Hifiasm? [y/n] " run_hifiasm
if [[ "$run_hifiasm" == "y" ]]; then
  Hifiasm_DIR="$ROOT/hifiasm"
  Hifiasm_OUT="$Hifiasm_DIR/${PROJECT_NAME}_from_gfa.fasta"
  if [[ ! -s "$Hifiasm_OUT" ]]; then
    echo "[INFO] Running Hifiasm..."
    mkdir -p "$Hifiasm_DIR"
    hifiasm -o "$Hifiasm_DIR/$PROJECT_NAME" -t "$THREADS" --ont "$MERGED_FILE" \
        2> "$Hifiasm_DIR/hifiasm_error.log"
    awk '/^S/{print ">"$2"\n"$3}' "$Hifiasm_DIR/${PROJECT_NAME}.bp.p_ctg.gfa" > "$Hifiasm_OUT"
  fi
  assemblies[hifiasm]="$Hifiasm_OUT"
else
  echo "[INFO] Skipping Hifiasm."
fi

conda deactivate

# NextDenovo
if [[ "$run_nextdenovo" != "y" && -z "${assemblies[nextdenovo]+x}" ]]; then
  ND_PREV=$(find "$ROOT/nextdenovo" -name "*nd.asm*.fasta" | head -n 1)
  if [[ -s "$ND_PREV" ]]; then
    assemblies[nextdenovo]="$ND_PREV"
    echo "[INFO] NextDenovo assembly trovato: $ND_PREV"
  fi
fi

# Canu
if [[ "$run_canu" != "y" && -z "${assemblies[canu]+x}" ]]; then
  if [[ -s "$Canu_OUT" ]]; then
    assemblies[canu]="$Canu_OUT"
    echo "[INFO] Canu assembly trovato: $Canu_OUT"
  fi
fi

# Flye
if [[ "$run_flye" != "y" && -z "${assemblies[flye]+x}" ]]; then
  Flye_OUT="$ROOT/flye/assembly.fasta"
  if [[ -s "$Flye_OUT" ]]; then
    assemblies[flye]="$Flye_OUT"
    echo "[INFO] Flye assembly trovato: $Flye_OUT"
  fi
fi

# Wtdbg2
if [[ "$run_wtdbg2" != "y" && -z "${assemblies[wtdbg2]+x}" ]]; then
  Wtdbg2_OUT="$ROOT/wtdbg2/wtdbg_assembly.fasta"
  if [[ -s "$Wtdbg2_OUT" ]]; then
    assemblies[wtdbg2]="$Wtdbg2_OUT"
    echo "[INFO] Wtdbg2 assembly trovato: $Wtdbg2_OUT"
  fi
fi

# Hifiasm
if [[ "$run_hifiasm" != "y" && -z "${assemblies[hifiasm]+x}" ]]; then
  HIFI_OUT="$ROOT/hifiasm/${PROJECT_NAME}_from_gfa.fasta"
  if [[ -s "$HIFI_OUT" ]]; then
    assemblies[hifiasm]="$HIFI_OUT"
    echo "[INFO] Hifiasm assembly trovato: $HIFI_OUT"
  fi
fi


echo "[INFO] Assemblies generated/found:"
for asm in "${!assemblies[@]}"; do
  echo " - $asm: ${assemblies[$asm]}"
done

# === POLISHING, QUAST, MAPPING ===

for asm in "${!assemblies[@]}"; do
  echo "=== PROCESSING $asm assembly ==="

  ASM_FASTA="${assemblies[$asm]}"
  ASM_BASE=$(basename "$ASM_FASTA" .fasta)
  POLISH_DIR="$ROOT/polished_assemblies/$asm"
  QUAST_DIR="$ROOT/quast_results/$asm"
  MAPPING_DIR="$ROOT/mapping/$asm"
  mkdir -p "$POLISH_DIR" "$QUAST_DIR" "$MAPPING_DIR"

  # 1) QUAST pre-polishing
  conda activate "$QC_ENV"
  echo "[INFO][$asm] Running QUAST pre-polishing..."
  quast.py -t "$THREADS" -o "$QUAST_DIR/quast_pre" "$ASM_FASTA"
  conda deactivate

  # 2) Polishing Racon 3 rounds
  # Input polished files
  CURRENT_FASTA="$ASM_FASTA"
  READS_FOR_POLISH="$MERGED_FILE"

  conda activate "$POLISHING_ENV"

  for i in {1..3}; do
    echo "[INFO][$asm] Racon polishing round $i..."
    MINIMAP_OUT="$POLISH_DIR/minimap_racon_${i}.paf"
    RACON_OUT="$POLISH_DIR/racon_round_${i}.fasta"

    minimap2 -t "$THREADS" -x map-ont "$CURRENT_FASTA" "$READS_FOR_POLISH" > "$MINIMAP_OUT"
    racon -t "$THREADS" "$READS_FOR_POLISH" "$MINIMAP_OUT" "$CURRENT_FASTA" > "$RACON_OUT"

    CURRENT_FASTA="$RACON_OUT"
  done

  # 3) Polishing Medaka
  echo "[INFO][$asm] Running Medaka polishing..."
  MEDAKA_OUT="$POLISH_DIR/medaka_output"
  mkdir -p "$MEDAKA_OUT"

  medaka_consensus -i "$READS_FOR_POLISH" -d "$CURRENT_FASTA" -o "$MEDAKA_OUT" -t "$THREADS" -m "$MODEL"
  MEDAKA_FASTA="$MEDAKA_OUT/consensus.fasta"

  conda deactivate

  # 4) QUAST post-polishing
  conda activate "$QC_ENV"
  echo "[INFO][$asm] Running QUAST post-polishing..."
  quast.py -t "$THREADS" -o "$QUAST_DIR/quast_post" "$MEDAKA_FASTA"
  conda deactivate

  # 5) Mapping with minimap2 + Qualimap
  conda activate "$POLISHING_ENV"
  echo "[INFO][$asm] Mapping polished assembly reads with minimap2..."
  BAM_OUT="$MAPPING_DIR/${asm}_mapped.bam"
  minimap2 -t "$THREADS" -ax map-ont "$MEDAKA_FASTA" "$READS_FOR_POLISH" | samtools sort -@ "$THREADS" -o "$BAM_OUT"
  samtools index "$BAM_OUT"
  conda deactivate
  echo "[INFO][$asm] Running Qualimap..."
  conda activate "$QC_ENV"
  qualimap bamqc -bam "$BAM_OUT" --java-mem-size=4G \
  -outdir "$MAPPING_DIR/qualimap_report" -nt "$THREADS" -outformat PDF:HTML

  conda deactivate

done
echo ""
echo "-------------------------------------"
echo "Assembly and Polishing pipeline completed successfully!"
echo "Check QUAST and Qualimap report to decide which de novo genome has the best coverage and quality!"


