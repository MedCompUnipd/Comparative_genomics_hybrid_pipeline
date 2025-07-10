#!/bin/bash
set -euo pipefail

# --- Conda setup ---
command_exists () {
    command -v "$1" &> /dev/null
}

if ! command_exists conda; then
    echo "[ERROR] Conda not found. Please ensure conda is installed and initialized."
    exit 1
fi
eval "$(conda shell.bash hook)"

echo "--- Assembly and Polishing Pipeline ---"
echo "---------------------------------------"

# === USER INPUT ===
read -p "==> Input reads (FASTQ file or directory): " READS
read -p "==> Root output directory for assembly steps: " ROOT
read -p "==> Estimated genome size (e.g. 235k): " GENOME_SIZE
read -p "==> Medaka model (e.g. r1041_e82_400bps_sup_variant_v4.2.0): " MODEL
read -p "==> Number of threads to use: " THREADS
read -p "==> Project name (e.g. HCMV): " PROJECT_NAME

# === ENVIRONMENT NAMES ===
PREPROCESSING_ENV="preprocessing_env"
ASSEMBLY_ENV="assembly_env"
POLISHING_ENV="polishing_env"
QC_ENV="qc_env"
NDN_ENV="ndn_env"

declare -A assemblies
mkdir -p "$ROOT"

# Check input reads
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
  echo "[ERROR] No merged FASTQ file found in input path."
  exit 1
fi

echo "[INFO] Using merged FASTQ: $MERGED_FILE"

# --------- NEXTDENOVO -------------
read -p "==> Run NextDenovo? [y/n] " run_nextdenovo
ND_DIR="$ROOT/nextdenovo"
ND_CFG="$ND_DIR/run.cfg"
mkdir -p "$ND_DIR"

if [[ "$run_nextdenovo" == "y" ]]; then
  echo "[INFO] Running NextDenovo..."
  conda deactivate || true
  conda activate "$NDN_ENV"

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
  conda deactivate
fi

# Check NextDenovo output
ND_OUTPUT=$(find "$ND_DIR" -name "*nd.asm*.fasta" 2>/dev/null | head -n 1)
if [[ -s "$ND_OUTPUT" ]]; then
  assemblies[nextdenovo]="$ND_OUTPUT"
  echo "[INFO] Found NextDenovo assembly: $ND_OUTPUT"
else
  echo "[INFO] No NextDenovo assembly found or used."
fi

# --------- CANU --------------------
read -p "==> Run Canu? [y/n] " run_canu
Canu_DIR="$ROOT/canu"
Canu_OUT="$Canu_DIR/${PROJECT_NAME}.contigs.fasta"
if [[ "$run_canu" == "y" ]]; then
  echo "[INFO] Running Canu..."
  conda deactivate || true
  conda activate "$ASSEMBLY_ENV"
  mkdir -p "$Canu_DIR"
  canu -p "$PROJECT_NAME" -d "$Canu_DIR" genomeSize="$GENOME_SIZE" \
       -nanopore-raw "$MERGED_FILE" useGrid=false minReadLength=500
  conda deactivate
fi
if [[ -s "$Canu_OUT" ]]; then
  assemblies[canu]="$Canu_OUT"
  echo "[INFO] Found Canu assembly: $Canu_OUT"
fi

# --------- FLYE --------------------
read -p "==> Run Flye? [y/n] " run_flye
Flye_DIR="$ROOT/flye"
Flye_OUT="$Flye_DIR/assembly.fasta"
if [[ "$run_flye" == "y" ]]; then
  echo "[INFO] Running Flye..."
  conda deactivate || true
  conda activate "$ASSEMBLY_ENV"
  mkdir -p "$Flye_DIR"
  flye --nano-raw "$MERGED_FILE" --out-dir "$Flye_DIR" \
       --threads "$THREADS" --genome-size "$GENOME_SIZE" \
       --min-overlap 1000 --asm-coverage 50 | tee "$Flye_DIR/flye.log"
  conda deactivate
fi
if [[ -s "$Flye_OUT" ]]; then
  assemblies[flye]="$Flye_OUT"
  echo "[INFO] Found Flye assembly: $Flye_OUT"
fi

# --------- WTDBG2 ------------------
read -p "==> Run Wtdbg2? [y/n] " run_wtdbg2
Wtdbg2_DIR="$ROOT/wtdbg2"
Wtdbg2_OUT="$Wtdbg2_DIR/wtdbg_assembly.fasta"
if [[ "$run_wtdbg2" == "y" ]]; then
  echo "[INFO] Running Wtdbg2..."
  mkdir -p "$Wtdbg2_DIR"
  conda deactivate || true
  conda activate "$ASSEMBLY_ENV"
  wtdbg2 -i "$MERGED_FILE" -o "$Wtdbg2_DIR/wtdbg_assembly" -t "$THREADS" -g "$GENOME_SIZE"
  wtpoa-cns -t "$THREADS" -i "$Wtdbg2_DIR/wtdbg_assembly.ctg.lay.gz" -fo "$Wtdbg2_OUT"
  conda deactivate
fi
if [[ -s "$Wtdbg2_OUT" ]]; then
  assemblies[wtdbg2]="$Wtdbg2_OUT"
  echo "[INFO] Found Wtdbg2 assembly: $Wtdbg2_OUT"
fi

# --------- HIFIASM -----------------
read -p "==> Run Hifiasm? [y/n] " run_hifiasm
Hifiasm_DIR="$ROOT/hifiasm"
Hifiasm_OUT="$Hifiasm_DIR/${PROJECT_NAME}_from_gfa.fasta"
if [[ "$run_hifiasm" == "y" ]]; then
  echo "[INFO] Running Hifiasm..."
  mkdir -p "$Hifiasm_DIR"
  conda deactivate || true
  conda activate "$ASSEMBLY_ENV"
  hifiasm -o "$Hifiasm_DIR/$PROJECT_NAME" -t "$THREADS" --ont "$MERGED_FILE" \
      2> "$Hifiasm_DIR/hifiasm_error.log"
  awk '/^S/{print ">"$2"\n"$3}' "$Hifiasm_DIR/${PROJECT_NAME}.bp.p_ctg.gfa" > "$Hifiasm_OUT"
  conda deactivate
fi
if [[ -s "$Hifiasm_OUT" ]]; then
  assemblies[hifiasm]="$Hifiasm_OUT"
  echo "[INFO] Found Hifiasm assembly: $Hifiasm_OUT"
fi

# === ASSEMBLY SUMMARY ===
echo ""
echo "[INFO] Assemblies ready for polishing:"
for asm in "${!assemblies[@]}"; do
  echo " - $asm: ${assemblies[$asm]}"
done

# === POLISHING + QUAST + QUALIMAP ===
for asm in "${!assemblies[@]}"; do
  echo "=== POLISHING [$asm] ==="
  ASM_FASTA="${assemblies[$asm]}"
  POLISH_DIR="$ROOT/polished_assemblies/$asm"
  QUAST_DIR="$ROOT/quast_results/$asm"
  MAPPING_DIR="$ROOT/mapping/$asm"
  mkdir -p "$POLISH_DIR" "$QUAST_DIR" "$MAPPING_DIR"

  # 1. QUAST pre-polishing
  conda activate "$QC_ENV"
  quast.py -t "$THREADS" -o "$QUAST_DIR/quast_pre" "$ASM_FASTA"
  conda deactivate

  # 2. Racon polishing (3 rounds)
  conda activate "$POLISHING_ENV"
  CURRENT_FASTA="$ASM_FASTA"
  for i in {1..3}; do
    echo "[INFO][$asm] Racon round $i..."
    MINIMAP_OUT="$POLISH_DIR/minimap_racon_${i}.paf"
    RACON_OUT="$POLISH_DIR/racon_round_${i}.fasta"
    minimap2 -t "$THREADS" -x map-ont "$CURRENT_FASTA" "$MERGED_FILE" > "$MINIMAP_OUT"
    racon -t "$THREADS" "$MERGED_FILE" "$MINIMAP_OUT" "$CURRENT_FASTA" > "$RACON_OUT"
    CURRENT_FASTA="$RACON_OUT"
  done

  # 3. Medaka
  MEDAKA_OUT="$POLISH_DIR/medaka_output"
  medaka_consensus -i "$MERGED_FILE" -d "$CURRENT_FASTA" -o "$MEDAKA_OUT" -t "$THREADS" -m "$MODEL" -b 16
  MEDAKA_FASTA="$MEDAKA_OUT/consensus.fasta"
  conda deactivate

  # 4. QUAST post-polishing
  conda activate "$QC_ENV"
  quast.py -t "$THREADS" -o "$QUAST_DIR/quast_post" "$MEDAKA_FASTA"
  conda deactivate

  # 5. Mapping + Qualimap
  conda activate "$POLISHING_ENV"
  BAM_OUT="$MAPPING_DIR/${asm}_mapped.bam"
  minimap2 -t "$THREADS" -ax map-ont "$MEDAKA_FASTA" "$MERGED_FILE" | samtools sort -@ "$THREADS" -o "$BAM_OUT"
  samtools index "$BAM_OUT"
  conda deactivate

  conda activate "$QC_ENV"
  qualimap bamqc -bam "$BAM_OUT" --java-mem-size=4G -outdir "$MAPPING_DIR/qualimap_report" -nt "$THREADS" -outformat PDF:HTML
  conda deactivate
done

echo ""
echo "---------------------------------------"
echo "Pipeline completed successfully!"
echo "Check QUAST and Qualimap reports to evaluate your assemblies."
