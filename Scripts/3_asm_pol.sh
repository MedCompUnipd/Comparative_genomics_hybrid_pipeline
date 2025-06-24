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
read -p "==> Input reads (FASTQ file or directory containing FASTQ files, e.g., merged_reads.fastq or /path/to/trimmed_fastqs): " READS_INPUT
read -p "==> Root output directory for all assembly steps: " ROOT_OUT_DIR
read -p "==> Estimated genome size (e.g. 235k, 4.5m): " GENOME_SIZE
read -p "==> Medaka model (e.g. r1041_e82_400bps_sup_variant_v4.2.0): " MEDAKA_MODEL
read -p "==> Threads to use: " THREADS
read -p "==> Project name (e.g. HCMV): " PROJECT_NAME

# === ENVIRONMENT NAMES ===
PREPROCESSING_ENV="preprocessing_env" # Used for minimap2, samtools, qualimap (if applicable)
ASSEMBLY_ENV="assembly_env" # Used for Flye, wtdbg2, Canu, NextDenovo, paralleltask
POLISHING_ENV="polishing_env" # Used for Racon, Medaka, BWA, Pilon
QC_ENV="qc_env" # Used for QUAST
NDN_ENV="ndn_env" # As per original script, if NextDenovo needs a Python env
CANU_ENV="canu_env" # As per original script, if Canu needs a Python env

declare -A assemblies # Associative array to store paths to assembly FASTA files

# === SETUP ===
mkdir -p "$ROOT_OUT_DIR" || { echo "[ERROR] Failed to create root output directory: $ROOT_OUT_DIR"; exit 1; }

# Control input reads
MERGED_FILE=""
if [[ -d "$READS_INPUT" ]]; then
  echo "[INFO] Input is a directory. Searching for merged FASTQ file..."
  # Look for common merged FASTQ names. Prioritize .fastq over .fq, and .gz over non-compressed.
  MERGED_FILE=$(find "$READS_INPUT" -type f \( -iname "*merged.fastq.gz" -o -iname "*merged.fastq" -o -iname "*merged.fq.gz" -o -iname "*merged.fq" \) | head -n 1)
  if [[ -z "$MERGED_FILE" ]]; then
      echo "[WARNING] No merged FASTQ file found in $READS_INPUT. Attempting to merge all .fastq files."
      MERGED_FILE="${ROOT_OUT_DIR}/${PROJECT_NAME}_merged_reads.fastq"
      echo "[INFO] Merging all .fastq files from $READS_INPUT to $MERGED_FILE..."
      find "$READS_INPUT" -type f -name "*.fastq" -print0 | xargs -0 cat > "$MERGED_FILE" || { echo "[ERROR] Failed to merge FASTQ files."; exit 1; }
      if [ ! -s "$MERGED_FILE" ]; then
          echo "[ERROR] No FASTQ files found or merged file is empty. Aborting."
          exit 1
      fi
  else
      echo "[INFO] Found merged FASTQ file: $MERGED_FILE"
  fi
elif [[ -f "$READS_INPUT" ]]; then
  MERGED_FILE="$READS_INPUT"
else
  echo "[ERROR] Input reads not found or not a valid directory/file: $READS_INPUT"
  exit 1
fi

READS_FOR_POLISH="$MERGED_FILE" # Use this for polishing steps


# --- Function to run an assembly method ---
run_assembly() {
  local method_name="$1"
  local output_dir="$2"
  local read_file="$3"
  local genome_size="$4"
  local threads="$5"
  local project_name="$6"
  local assembly_fasta_path=""

  echo ""
  echo "=== Running $method_name Assembly ==="
  mkdir -p "$output_dir" || { echo "[ERROR] Failed to create output directory: $output_dir"; return 1; }
  pushd "$output_dir" > /dev/null

  case "$method_name" in
    "Flye")
      echo "[INFO] Activating $ASSEMBLY_ENV for Flye..."
      conda activate "$ASSEMBLY_ENV" || { echo "[ERROR] Failed to activate 'assembly_env'."; popd; return 1; }
      echo "[INFO] Running Flye..."
      flye --nano-hq "$read_file" --genome-size "$genome_size" --threads "$threads" --out-dir "./flye_output" || { echo "[ERROR] Flye assembly failed."; conda deactivate; popd; return 1; }
      assembly_fasta_path="./flye_output/assembly.fasta"
      conda deactivate
      ;;
    "wtdbg2")
      echo "[INFO] Activating $ASSEMBLY_ENV for wtdbg2..."
      conda activate "$ASSEMBLY_ENV" || { echo "[ERROR] Failed to activate 'assembly_env'."; popd; return 1; }
      echo "[INFO] Running wtdbg2..."
      wtdbg2 -i "$read_file" -o "wtdbg_output" -t "$threads" -g "$genome_size" || { echo "[ERROR] wtdbg2 assembly failed."; conda deactivate; popd; return 1; }
      wtpoa-cns -t "$threads" -i "wtdbg_output.ctg.lay.gz" -f "$read_file" -o "wtdbg_output.ctg.fasta" || { echo "[ERROR] wtpoa-cns failed."; conda deactivate; popd; return 1; }
      assembly_fasta_path="wtdbg_output.ctg.fasta"
      conda deactivate
      ;;
    "Canu")
      echo "[INFO] Activating $CANU_ENV (or assuming Canu is in PATH) for Canu..."
      conda activate "$CANU_ENV" || { echo "[ERROR] Failed to activate 'canu_env'."; popd; return 1; }
      # Canu is assumed to be in PATH due to manual installation in setup script
      if ! command_exists canu; then
          echo "[ERROR] 'canu' command not found. Please ensure Canu is installed and in your PATH."
          conda deactivate
          popd; return 1;
      fi
      echo "[INFO] Running Canu..."
      # Canu uses genomeSize= and -nanopore-raw for ONT reads
      canu -p "$project_name" -d "./canu_output" genomeSize="$genome_size" -nanopore-raw "$read_file" useGrid=false -j "$threads" || { echo "[ERROR] Canu assembly failed."; conda deactivate; popd; return 1; }
      assembly_fasta_path="./canu_output/$project_name.contigs.fasta" # or .unitigs.fasta depending on assembly type
      conda deactivate
      ;;
    "Hifiasm")
      echo "[INFO] Hifiasm (manual installation expected)."
      # Hifiasm is assumed to be in PATH due to manual installation in setup script
      if ! command_exists hifiasm; then
          echo "[ERROR] 'hifiasm' command not found. Please ensure Hifiasm is installed and in your PATH."
          popd; return 1;
      fi
      echo "[INFO] Running Hifiasm..."
      hifiasm -o "./hifiasm_output.asm" -t "$threads" "$read_file" || { echo "[ERROR] Hifiasm assembly failed."; popd; return 1; }
      assembly_fasta_path="./hifiasm_output.asm.bp.p_ctg.gfa" # This is a GFA, conversion to FASTA needed for polishing
      echo "[INFO] Hifiasm output is GFA. Converting to FASTA..."
      awk '/^S/{print ">"$2"\\n"$3}' "$assembly_fasta_path" > "./hifiasm_output.fasta" || { echo "[ERROR] GFA to FASTA conversion failed."; popd; return 1; }
      assembly_fasta_path="./hifiasm_output.fasta"
      ;;
    "NextDenovo")
      echo "[INFO] NextDenovo (manual installation expected)."
      # NextDenovo is assumed to be in PATH due to manual installation in setup script
      if ! command_exists NextDenovo; then
          echo "[ERROR] 'NextDenovo' command not found. Please ensure NextDenovo is installed and in your PATH."
          popd; return 1;
      fi
      echo "[INFO] Running NextDenovo. This requires a config file."
      echo "[INSTRUCTIONS] You need to create a NextDenovo config file (e.g., input.cfg) manually."
      echo "  Example config structure:"
      echo "  [General]"
      echo "  job_type = local"
      echo "  job_prefix = nextdenovo"
      echo "  task = all"
      echo "  read_type = ont"
      echo "  input_type = raw"
      echo "  input_fofn = input.fofn  (a file listing your FASTQ file: $read_file)"
      echo "  genome_size = $genome_size"
      echo "  threads = $threads"
      echo "  parallel_jobs = 4" # Adjust as needed
      echo "  minimap2_options = -x ont"
      echo ""
      echo "  Create a file named 'input.fofn' in the NextDenovo output directory with the path to your merged FASTQ file:"
      echo "  echo \"$read_file\" > \"nextdenovo_output/input.fofn\""
      echo ""
      read -p "Please provide the path to your NextDenovo config file (e.g., /path/to/nextdenovo_output/input.cfg): " NDN_CONFIG_FILE
      if [[ ! -f "$NDN_CONFIG_FILE" ]]; then
          echo "[ERROR] NextDenovo config file not found: $NDN_CONFIG_FILE. Aborting."
          popd; return 1;
      fi
      echo "[INFO] Creating input.fofn for NextDenovo..."
      NDN_OUT_DIR="${output_dir}/nextdenovo_output"
      mkdir -p "$NDN_OUT_DIR"
      echo "$read_file" > "${NDN_OUT_DIR}/input.fofn"

      conda activate "$NDN_ENV" || { echo "[ERROR] Failed to activate 'ndn_env'."; popd; return 1; }
      NextDenovo "$NDN_CONFIG_FILE" || { echo "[ERROR] NextDenovo assembly failed."; conda deactivate; popd; return 1; }
      assembly_fasta_path="${NDN_OUT_DIR}/03.ctg_graph/ndn.ctg.fasta"
      conda deactivate
      ;;
    *)
      echo "[ERROR] Unknown assembly method: $method_name"
      popd; return 1;
      ;;
  esac
  popd > /dev/null

  if [[ -f "$output_dir/$(basename "$assembly_fasta_path")" ]]; then # Check if the expected FASTA is there
      assemblies["$method_name"]="$output_dir/$(basename "$assembly_fasta_path")"
      echo "[SUCCESS] $method_name assembly completed. Output: ${assemblies["$method_name"]}"
  else
      echo "[ERROR] $method_name assembly failed to produce expected FASTA file."
      return 1
  fi
  return 0
}

# --- Function to run polishing steps ---
run_polishing() {
  local asm_method="$1"
  local initial_fasta="$2"
  local reads_for_polish="$3"
  local threads="$4"
  local medaka_model="$5"

  local POLISH_DIR="${ROOT_OUT_DIR}/${asm_method}_polish"
  mkdir -p "$POLISH_DIR" || { echo "[ERROR] Failed to create polishing directory: $POLISH_DIR"; return 1; }
  echo ""
  echo "=== Running Polishing for $asm_method Assembly ==="

  local CURRENT_FASTA="$initial_fasta"

  # 1) Polishing Racon (3 rounds)
  echo "[INFO][$asm_method] Running Racon polishing (3 rounds)..."
  conda activate "$POLISHING_ENV" || { echo "[ERROR] Failed to activate 'polishing_env'."; return 1; }
  for i in {1..3}; do
    echo "  → Racon round $i..."
    local MINIMAP_OUT="${POLISH_DIR}/minimap_racon_${i}.paf"
    local RACON_OUT="${POLISH_DIR}/racon_round_${i}.fasta"

    minimap2 -t "$threads" -x map-ont "$CURRENT_FASTA" "$reads_for_polish" > "$MINIMAP_OUT" || { echo "[ERROR] minimap2 for Racon round $i failed."; conda deactivate; return 1; }
    racon -t "$threads" "$reads_for_polish" "$MINIMAP_OUT" "$CURRENT_FASTA" > "$RACON_OUT" || { echo "[ERROR] racon round $i failed."; conda deactivate; return 1; }
    CURRENT_FASTA="$RACON_OUT"
  done
  echo "[INFO][$asm_method] Racon polishing completed. Output: $CURRENT_FASTA"
  conda deactivate

  # 2) Polishing Medaka
  echo "[INFO][$asm_method] Running Medaka polishing..."
  conda activate "$POLISHING_ENV" || { echo "[ERROR] Failed to activate 'polishing_env'."; return 1; }
  local MEDAKA_OUT_DIR="${POLISH_DIR}/medaka_output"
  mkdir -p "$MEDAKA_OUT_DIR" || { echo "[ERROR] Failed to create Medaka output directory."; conda deactivate; return 1; }

  medaka_consensus -i "$reads_for_polish" -d "$CURRENT_FASTA" -o "$MEDAKA_OUT_DIR" -t "$threads" -m "$medaka_model" || { echo "[ERROR] Medaka polishing failed."; conda deactivate; return 1; }
  local MEDAKA_FASTA="${MEDAKA_OUT_DIR}/consensus.fasta"
  echo "[INFO][$asm_method] Medaka polishing completed. Output: $MEDAKA_FASTA"
  conda deactivate

  # 3) QUAST post-polishing
  echo "[INFO][$asm_method] Running QUAST post-polishing..."
  conda activate "$QC_ENV" || { echo "[ERROR] Failed to activate 'qc_env'."; return 1; }
  local QUAST_DIR="${POLISH_DIR}/quast_post"
  mkdir -p "$QUAST_DIR" || { echo "[ERROR] Failed to create QUAST directory."; conda deactivate; return 1; }
  quast.py -t "$threads" -o "$QUAST_DIR" "$MEDAKA_FASTA" || { echo "[ERROR] QUAST post-polishing failed."; conda deactivate; return 1; }
  echo "[INFO][$asm_method] QUAST post-polishing completed. Reports in $QUAST_DIR"
  conda deactivate

  # 4) Mapping with minimap2 + Qualimap
  echo "[INFO][$asm_method] Mapping polished assembly reads with minimap2 and running Qualimap..."
  local MAPPING_DIR="${POLISH_DIR}/mapping"
  mkdir -p "$MAPPING_DIR" || { echo "[ERROR] Failed to create mapping directory."; return 1; }

  conda activate "$PREPROCESSING_ENV" || { echo "[ERROR] Failed to activate 'preprocessing_env' for mapping."; return 1; }
  local BAM_OUT="${MAPPING_DIR}/${asm_method}_mapped.bam"
  minimap2 -ax map-ont "$MEDAKA_FASTA" "$reads_for_polish" | samtools view -bS - | samtools sort -o "$BAM_OUT" - || { echo "[ERROR] Minimap2/Samtools mapping failed."; conda deactivate; return 1; }
  samtools index "$BAM_OUT" || { echo "[ERROR] Samtools indexing failed."; conda deactivate; return 1; }
  echo "[INFO][$asm_method] Mapping completed. BAM file: $BAM_OUT"
  conda deactivate

  conda activate "$QC_ENV" || { echo "[ERROR] Failed to activate 'qc_env' for Qualimap."; return 1; }
  qualimap bamqc -bam "$BAM_OUT" -outdir "${MAPPING_DIR}/qualimap_report" --java-mem-size=4G || { echo "[ERROR] Qualimap failed."; conda deactivate; return 1; }
  echo "[INFO][$asm_method] Qualimap QC completed. Reports in ${MAPPING_DIR}/qualimap_report"
  conda deactivate

  # Store the final polished assembly path
  assemblies["${asm_method}_polished"]="$MEDAKA_FASTA"
  return 0
}

# === Main Assembly Workflow ===
echo ""
echo "--- Assembly Methods Selection ---"
echo "Available methods: Flye, wtdbg2, Canu, Hifiasm, NextDenovo"
read -p "Enter assembly methods to run, separated by space (e.g., Flye Canu Hifiasm): " -a SELECTED_ASSEMBLIES

if [ ${#SELECTED_ASSEMBLIES[@]} -eq 0 ]; then
    echo "[ERROR] No assembly methods selected. Exiting."
    exit 1
fi

for ASM in "${SELECTED_ASSEMBLIES[@]}"; do
    ASM_LOWER=$(echo "$ASM" | tr '[:upper:]' '[:lower:]')
    ASM_OUT_DIR="${ROOT_OUT_DIR}/${ASM_LOWER}_assembly"
    run_assembly "$ASM" "$ASM_OUT_DIR" "$READS_FOR_POLISH" "$GENOME_SIZE" "$THREADS" "$PROJECT_NAME" || echo "[ERROR] $ASM assembly failed. See previous messages."
done

echo ""
echo "--- Polishing Workflow ---"
if [ ${#assemblies[@]} -eq 0 ]; then
    echo "[WARNING] No successful assemblies to polish. Skipping polishing step."
else
    for ASM_KEY in "${!assemblies[@]}"; do
        if [[ "$ASM_KEY" != *_polished ]]; then # Only polish raw assemblies, not already polished ones
            ASM_PATH="${assemblies[$ASM_KEY]}"
            run_polishing "$ASM_KEY" "$ASM_PATH" "$READS_FOR_POLISH" "$THREADS" "$MEDAKA_MODEL" || echo "[ERROR] Polishing for $ASM_KEY failed. See previous messages."
        fi
    done
fi

echo ""
echo "-------------------------------------"
echo "Assembly and Polishing pipeline completed successfully!"
echo "Final polished assemblies:"
for ASM_KEY in "${!assemblies[@]}"; do
    echo "  - $ASM_KEY: ${assemblies[$ASM_KEY]}"
done
