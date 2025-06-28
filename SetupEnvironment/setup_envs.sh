#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e
# Exit if any command in a pipeline fails.
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

# Check if conda is available and load its hook
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
    mkdir -p "$DORADO_DIR" || { echo "[ERROR] Could not create directory $DORADO_DIR. Exiting."; exit 1; }
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
    {
        echo ""
        echo "# Added by Dorado installer"
        echo "export PATH=\"\$PATH:$DORADO_BIN\""
    } >> "$SHELL_PROFILE"
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
    {
        echo ""
        echo "# Added by Canu installer"
        echo "export PATH=\"\$PATH:$CANU_BIN\""
    } >> "$SHELL_PROFILE"
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
    {
        echo ""
        echo "# Added by Hifiasm installer"
        echo "export PATH=\"\$PATH:$(dirname "$HIFIASM_BIN")\""
    } >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE"
fi
export PATH="$PATH:$(dirname "$HIFIASM_BIN")"
echo ""

# === 4. ndn_env ===
echo "[4/9] Creating environment: ndn_env"
if ! conda env list | grep -q "ndn_env"; then
    conda create -y -n ndn_env python=3.10 -c bioconda -c conda-forge minimap2 samtools || {
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
    {
        echo ""
        echo "# Added by NextDenovo installer"
        echo "export PYTHONPATH=\"\$PYTHONPATH:$NDN_LIB_PATH\""
    } >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE"
fi
export PYTHONPATH="$PYTHONPATH:$NDN_LIB_PATH"
echo "[INFO] Setup complete."

# === 5. polishing_env ===
echo "[5/9] Creating environment: polishing_env"
if ! conda env list | grep -q "polishing_env"; then
    conda create -y -n polishing_env -c bioconda -c nanoporetech -c conda-forge \
        medaka racon quast || \
        { echo "[ERROR] Failed to create 'polishing_env' Conda environment. Exiting."; exit 1; }
else
    echo "[INFO] Environment 'polishing_env' already exists."
fi
echo "[INFO] 'polishing_env' setup complete."
echo ""

# === 6. qc_env ===
echo "[6/9] Creating environment: qc_env"
if ! conda env list | grep -q "qc_env"; then
    conda create -y -n qc_env -c bioconda -c conda-forge \
        nanoplot fastqc qualimap bcftools samtools fastp || \
        { echo "[ERROR] Failed to create 'qc_env' Conda environment. Exiting."; exit 1; }
else
    echo "[INFO] Environment 'qc_env' already exists."
fi
echo "[INFO] 'qc_env' setup complete."
echo ""

# === 7. illuminareads_env ===
echo "[7/9] Creating environment: illuminareads_env"
if ! conda env list | grep -q "illuminareads_env"; then
    conda create -y -n illuminareads_env -c bioconda -c conda-forge \
        bwa pilon || \
        { echo "[ERROR] Failed to create 'illuminareads_env' Conda environment. Exiting."; exit 1; }
else
    echo "[INFO] Environment 'illuminareads_env' already exists."
fi
echo "[INFO] 'illuminareads_env' setup complete."
echo ""

# === 8. mummer_env ===
echo "[8/9] Creating environment: mummer_env"
if ! conda env list | grep -q "mummer_env"; then
    conda create -y -n mummer_env -c bioconda -c conda-forge \
        mummer bedtools blast || \
        { echo "[ERROR] Failed to create 'mummer_env' Conda environment. Exiting."; exit 1; }
else
    echo "[INFO] Environment 'mummer_env' already exists."
fi
echo "[INFO] 'mummer_env' setup complete."
echo ""

# === 9. fasta3_env ===
echo "[9/9] Creating environment: fasta3_env"
if ! conda env list | grep -q "fasta3_env"; then
    conda create -y -n fasta3_env -c conda-forge -c bioconda \
        pandas openpyxl biopython fasta3 || \
        { echo "[ERROR] Failed to create 'fasta3_env' Conda environment. Exiting."; exit 1; }
else
    echo "[INFO] Environment 'fasta3_env' already exists."
fi
echo "[INFO] 'fasta3_env' setup complete."
echo ""

echo "#####################################################"
echo "#           SETUP COMPLETE!                       #"
echo "#####################################################"
echo "All specified Conda environments have been created or verified."
echo "Manually installed tools (Dorado, Canu, Hifiasm, NextDenovo) should be in your PATH/PYTHONPATH."
echo "Please remember to restart your terminal or run 'source $SHELL_PROFILE' "
echo "to ensure all PATH and PYTHONPATH changes are active in new sessions."
echo ""
echo "Good luck with your work!"
