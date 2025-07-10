#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e
set -o pipefail

echo "#################################################################"
echo "# Hybrid Comparative Genomic Pipeline Environment Setup Script #"
echo "#################################################################"
echo ""
echo "This script will install all necessary software via Conda or provide instructions for manual installations."
echo "Please ensure Miniconda/Anaconda is already installed and initialized before proceeding."
echo "No sudo privileges are used by this script for software installation."
echo ""

# --- Function to check if a command exists ---
command_exists () {
    command -v "$1" &> /dev/null
}

# --- Determine shell profile file ---
SHELL_PROFILE=""
if [[ -n "$BASH_VERSION" ]]; then
    SHELL_PROFILE="$HOME/.bashrc"
elif [[ -n "$ZSH_VERSION" ]]; then
    SHELL_PROFILE="$HOME/.zshrc"
else
    echo "[WARNING] Unknown shell. PATH modifications might not be persistent. Assuming .bashrc."
    SHELL_PROFILE="$HOME/.bashrc"
fi
echo "[INFO] Detected shell profile: $SHELL_PROFILE"

# Load conda environment
echo "[INFO] Loading Conda hook..."
if ! command_exists conda; then
    echo "[ERROR] Conda is not installed or not in PATH. Please install Miniconda/Anaconda first."
    exit 1
fi
eval "$(conda shell.bash hook)"

echo "[INFO] Starting Conda environment setup..."
echo ""

echo "--- INSTALLING SOFTWARE AND CONDA ENVIRONMENTS ---"
echo "------------------------------------------------"

# === 1. preprocessing_env ===
echo "[1/9] Creating environment: preprocessing_env"
if ! conda env list | grep -q "preprocessing_env"; then
    conda create -y -n preprocessing_env -c bioconda -c conda-forge minimap2=2.28 || {
        echo "[ERROR] Failed to create 'preprocessing_env'. Exiting."
        exit 1
    }
else
    echo "[INFO] Environment 'preprocessing_env' already exists."
fi

echo "[INFO] Installing pod5 in preprocessing_env"
conda activate preprocessing_env || { echo "[ERROR] Failed to activate 'preprocessing_env'. Exiting."; exit 1; }
pip install pod5 || { echo "[ERROR] Failed to install pod5 in 'preprocessing_env'. Exiting."; exit 1; }
conda deactivate || { echo "[ERROR] Failed to deactivate 'preprocessing_env'. Exiting."; exit 1; }
echo "[INFO] 'preprocessing_env' setup complete."
echo ""

# === Dorado install ===
echo "[INFO] Downloading and extracting Dorado v0.9.6..."
read -p "Path to the directory where Dorado should be installed (e.g., $HOME/opt): " DORADO_DIR
DORADO_DIR="${DORADO_DIR%/}"
DORADO_VERSION="0.9.6"
DORADO_TAR_GZ="dorado-$DORADO_VERSION-linux-x64.tar.gz"
DORADO_URL="https://cdn.oxfordnanoportal.com/software/analysis/$DORADO_TAR_GZ"
DORADO_BIN="$DORADO_DIR/dorado-$DORADO_VERSION-linux-x64/bin"
DORADO_DOWNLOAD_FILENAME="dorado.tar.gz"

if [[ -x "$DORADO_BIN/dorado" ]]; then
    echo "[INFO] Dorado is already installed in $DORADO_BIN, skipping download."
else
    mkdir -p "$DORADO_DIR"
    pushd "$DORADO_DIR" > /dev/null
    echo "[INFO] Downloading Dorado..."
    curl -fL "$DORADO_URL" -o "$DORADO_DOWNLOAD_FILENAME" || { echo "[ERROR] Failed to download Dorado. Exiting."; popd > /dev/null; exit 1; }
    tar -xzf "$DORADO_DOWNLOAD_FILENAME" || { echo "[ERROR] Extraction failed. Exiting."; popd > /dev/null; exit 1; }
    rm -f "$DORADO_DOWNLOAD_FILENAME"
    popd > /dev/null
    [[ -x "$DORADO_BIN/dorado" ]] || { echo "[ERROR] Dorado binary missing after extraction. Exiting."; exit 1; }
    echo "[INFO] Dorado v$DORADO_VERSION installed successfully."
fi

if [[ ":$PATH:" != *":$DORADO_BIN:"* ]]; then
    echo "[INFO] Adding Dorado to PATH in $SHELL_PROFILE"
    echo -e "\n# Added by Dorado installer\nexport PATH=\"\$PATH:$DORADO_BIN\"" >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE"
fi
export PATH="$PATH:$DORADO_BIN"
echo ""

# === 2. kraken2_env ===
echo "[2/9] Creating environment: kraken2_env"
if ! conda env list | grep -q "kraken2_env"; then
    conda create -y -n kraken2_env -c bioconda -c conda-forge kraken2 seqtk || {
        echo "[ERROR] Failed to create 'kraken2_env'. Exiting."
        exit 1
    }
else
    echo "[INFO] Environment 'kraken2_env' already exists."
fi

# === 3. assembly_env ===
echo "[3/9] Creating environment: assembly_env"
if ! conda env list | grep -q "assembly_env"; then
    conda create -y -n assembly_env -c bioconda -c conda-forge flye wtdbg paralleltask openjdk=8 || {
        echo "[ERROR] Failed to create 'assembly_env'. Exiting."
        exit 1
    }
else
    echo "[INFO] Environment 'assembly_env' already exists."
fi
echo ""

# === Canu install ===
echo "[INFO] Installing Canu v2.2..."
read -p "Path to the directory where Canu should be installed (e.g., $HOME/opt): " CANU_DIR
CANU_DIR="${CANU_DIR%/}"
CANU_VERSION="2.2"
CANU_ARCHIVE="canu-$CANU_VERSION.Linux-amd64.tar.xz"
CANU_URL="https://github.com/marbl/canu/releases/download/v$CANU_VERSION/$CANU_ARCHIVE"
CANU_BIN="$CANU_DIR/canu-$CANU_VERSION/bin"
CANU_DOWNLOAD_PATH="/tmp/$CANU_ARCHIVE"

if [[ -x "$CANU_BIN/canu" ]]; then
    echo "[INFO] Canu is already installed in $CANU_BIN."
else
    mkdir -p "$CANU_DIR"
    curl -fL "$CANU_URL" -o "$CANU_DOWNLOAD_PATH" || { echo "[ERROR] Download failed. Exiting."; exit 1; }
    tar -xJf "$CANU_DOWNLOAD_PATH" -C "$CANU_DIR"
    rm -f "$CANU_DOWNLOAD_PATH"
    [[ -x "$CANU_BIN/canu" ]] || { echo "[ERROR] Canu binary missing. Exiting."; exit 1; }
    echo "[INFO] Canu v$CANU_VERSION installed successfully."
fi

if [[ ":$PATH:" != *":$CANU_BIN:"* ]]; then
    echo "[INFO] Adding Canu to PATH in $SHELL_PROFILE"
    echo -e "\n# Added by Canu installer\nexport PATH=\"\$PATH:$CANU_BIN\"" >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE"
fi
export PATH="$PATH:$CANU_BIN"
echo ""

# === Hifiasm install ===
echo "[INFO] Installing Hifiasm from source..."
read -p "Path to the directory where Hifiasm should be installed (e.g., $HOME/opt): " HIFIASM_DIR
HIFIASM_DIR="${HIFIASM_DIR%/}"
HIFIASM_REPO="https://github.com/chhylp123/hifiasm.git"
HIFIASM_CLONE_DIR="$HIFIASM_DIR/hifiasm"
HIFIASM_BIN="$HIFIASM_CLONE_DIR/hifiasm"

if [[ -x "$HIFIASM_BIN" ]]; then
    echo "[INFO] Hifiasm is already built."
else
    mkdir -p "$HIFIASM_DIR"
    git clone "$HIFIASM_REPO" "$HIFIASM_CLONE_DIR"
    pushd "$HIFIASM_CLONE_DIR"
    make -j$(nproc)
    popd
    [[ -x "$HIFIASM_BIN" ]] || { echo "[ERROR] Compilation failed. Exiting."; exit 1; }
    echo "[INFO] Hifiasm installed successfully."
fi

if [[ ":$PATH:" != *":$(dirname "$HIFIASM_BIN"):"* ]]; then
    echo "[INFO] Adding Hifiasm to PATH in $SHELL_PROFILE"
    echo -e "\n# Added by Hifiasm installer\nexport PATH=\"\$PATH:$(dirname "$HIFIASM_BIN")\"" >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE"
fi
export PATH="$PATH:$(dirname "$HIFIASM_BIN")"
echo ""

# === 4. ndn_env ===



# === 4. ndn_env ===
echo "[4/9] Creating environment: ndn_env"
if ! conda env list | grep -q "ndn_env"; then
    conda create -y -n ndn_env python=3.10 minimap2 samtools -c bioconda -c conda-forge || {
        echo "[ERROR] Failed to create 'ndn_env'. Exiting."
        exit 1
    }
else
    echo "[INFO] Environment 'ndn_env' already exists."
fi

echo "[INFO] Cloning and setting up NextDenovo..."
read -p "Path to the directory where NextDenovo should be installed (e.g., $HOME/opt): " NDN_DIR
NDN_DIR="${NDN_DIR%/}"
NEXTDENOVO_REPO="https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz"
NEXTDENOVO_CLONE_DIR="$NDN_DIR/NextDenovo"
NDN_BIN="$NEXTDENOVO_CLONE_DIR/nextDenovo"

if [[ -f "$NDN_BIN" ]]; then
    echo "[INFO] NextDenovo already cloned."
else
    mkdir -p "$NDN_DIR"
    NDN_TGZ="$NDN_DIR/NextDenovo.tgz"
    wget -O "$NDN_TGZ" "$NEXTDENOVO_REPO" || { echo "[ERROR] Failed to download NextDenovo. Exiting."; exit 1; }
    tar -xzf "$NDN_TGZ" -C "$NDN_DIR" || { echo "[ERROR] Failed to extract NextDenovo. Exiting."; exit 1; }
    rm "$NDN_TGZ"
fi

NDN_LIB_PATH="$NEXTDENOVO_CLONE_DIR/lib"
if [[ ":$PYTHONPATH:" != *":$NDN_LIB_PATH:"* ]]; then
    echo "[INFO] Adding NextDenovo lib to PYTHONPATH in $SHELL_PROFILE"
    echo -e "\n# Added by NextDenovo installer\nexport PYTHONPATH=\"\$PYTHONPATH:$NDN_LIB_PATH\"" >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE"
fi
export PYTHONPATH="$PYTHONPATH:$NDN_LIB_PATH"
echo "[INFO] Setup complete."
echo ""
if [[ ":$PATH:" != *":$NEXTDENOVO_CLONE_DIR:"* ]]; then
    echo "[INFO] Adding NextDenovo to PATH in $SHELL_PROFILE"
    echo -e "\n# Added by NextDenovo installer\nexport PATH=\"\$PATH:$NEXTDENOVO_CLONE_DIR\"" >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE"
fi
export PATH="$PATH:$NEXTDENOVO_CLONE_DIR"


# === 5. rnaseq_env ===

echo "[5/9] Creating environment: rnaseq_env"
if ! conda env list | grep -q "rnaseq_env"; then
    conda create -y -n rnaseq_env -c bioconda -c conda-forge hisat2 samtools subread gffread || {
        echo "[ERROR] Failed to create 'rnaseq_env'. Exiting."
        exit 1
    }
else
    echo "[INFO] Environment 'rnaseq_env' already exists."
fi
echo ""


# === 6-10. Other environments ===
declare -A envs=(
  ["polishing_env"]="medaka racon quast"
  ["qc_env"]="nanoplot fastqc qualimap bcftools samtools fastp"
  ["illuminareads_env"]="bwa pilon"
  ["mummer_env"]="mummer bedtools blast"
  ["fasta3_env"]="pandas openpyxl biopython fasta3 xlsxwriter"
)

i=5
for env in "${!envs[@]}"; do
  echo "[$i/9] Creating environment: $env"
  if ! conda env list | grep -q "$env"; then
    conda create -y -n "$env" -c bioconda -c conda-forge ${envs[$env]} || {
        echo "[ERROR] Failed to create '$env'. Exiting."
        exit 1
    }
  else
    echo "[INFO] Environment '$env' already exists."
  fi
  echo "[INFO] '$env' setup complete."
  echo ""
  ((i++))
done

echo "#####################################################"
echo "#           SETUP COMPLETE!                       #"
echo "#####################################################"
echo "All specified Conda environments have been created or verified."
echo "Manually installed tools (Dorado, Canu, Hifiasm, NextDenovo) should be in your PATH/PYTHONPATH."
echo "Please remember to restart your terminal or run 'source $SHELL_PROFILE' "
echo "to ensure all PATH and PYTHONPATH changes are active in new sessions."
echo ""
echo "Good luck with your work!"
