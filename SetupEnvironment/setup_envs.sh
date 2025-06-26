#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e
# Exit if any command in a pipeline fails.
set -o pipefail

echo "#################################################################"
echo "# Hybrid Comparative Genomic Pipeline Environment Setup Script #"
echo "################################################################"
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
# Note: The script will directly append to ~/.bashrc as per the user's latest input.
# This variable is still used for informational messages.
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
# Initialize Conda for the current shell session
eval "$(conda shell.bash hook)"
echo "Starting Conda environment setup..."

echo ""
echo "--- INSTALLING SOFTWARE AND CONDA ENVIRONMENTS ---"
echo "------------------------------------------------"

# === 1. preprocessing_env ===
echo "[1/9] Creating environment: preprocessing_env"
if ! conda env list | grep -q "preprocessing_env"; then
    conda create -y -n preprocessing_env -c bioconda -c conda-forge minimap2=2.28 || \
        { echo "[ERROR] Failed to create 'preprocessing_env' Conda environment. Exiting."; exit 1; }
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
# Prompt for Dorado installation directory
read -p "Path to the directory where Dorado should be installed (e.g., $HOME/opt): " DORADO_DIR
# Remove trailing slash if present
DORADO_DIR="${DORADO_DIR%/}"

DORADO_VERSION="0.9.6"
DORADO_TAR_GZ="dorado-$DORADO_VERSION-linux-x64.tar.gz"
DORADO_URL="https://cdn.oxfordnanoportal.com/software/analysis/$DORADO_TAR_GZ"
DORADO_BIN="$DORADO_DIR/dorado-$DORADO_VERSION-linux-x64/bin"
DORADO_DOWNLOAD_FILENAME="dorado.tar.gz"

if [[ -x "$DORADO_BIN/dorado" ]]; then
    echo "[INFO] Dorado is already installed in $DORADO_BIN, skipping download."
else
    echo "[INSTALLATION] Attempting to install Dorado v$DORADO_VERSION..."
    mkdir -p "$DORADO_DIR" || { echo "[ERROR] Failed to create directory $DORADO_DIR for Dorado. Exiting."; exit 1; }
    
    pushd "$DORADO_DIR" > /dev/null || { echo "[ERROR] Failed to enter directory $DORADO_DIR. Exiting."; exit 1; }
    
    echo "[INFO] Downloading Dorado from $DORADO_URL to $DORADO_DIR/$DORADO_DOWNLOAD_FILENAME..."
    if curl -fL "$DORADO_URL" -o "$DORADO_DOWNLOAD_FILENAME"; then
        if [[ -s "$DORADO_DOWNLOAD_FILENAME" ]]; then
            echo "[INFO] Extracting Dorado from $DORADO_DOWNLOAD_FILENAME..."
            if tar -xzf "$DORADO_DOWNLOAD_FILENAME"; then
                rm -f "$DORADO_DOWNLOAD_FILENAME" || echo "[WARNING] Failed to remove temporary Dorado tarball: $DORADO_DOWNLOAD_FILENAME"
                if [[ -x "$DORADO_BIN/dorado" ]]; then
                    echo "[INFO] Dorado v$DORADO_VERSION installed successfully."
                    echo "[INFO] Dorado version check: $("$DORADO_BIN/dorado" --version 2>&1)"
                else
                    echo "[ERROR] Dorado executable not found at $DORADO_BIN/dorado after extraction. Please check the extracted directory structure. Exiting.";
                    popd > /dev/null; exit 1
                fi
            else
                echo "[ERROR] Failed to extract Dorado from $DORADO_DOWNLOAD_FILENAME. Exiting.";
                rm -f "$DORADO_DOWNLOAD_FILENAME"
                popd > /dev/null; exit 1
            fi
        else
            echo "[ERROR] Downloaded Dorado file is empty or corrupted: $DORADO_DIR/$DORADO_DOWNLOAD_FILENAME. Exiting.";
            rm -f "$DORADO_DOWNLOAD_FILENAME"
            popd > /dev/null; exit 1
        fi
    else
        echo "[ERROR] Failed to download Dorado v$DORADO_VERSION from $DORADO_URL. Check network and URL. Exiting.";
        popd > /dev/null; exit 1
    fi
    popd > /dev/null || { echo "[ERROR] Failed to return from Dorado installation directory. Exiting."; exit 1; }
fi

# Add Dorado to PATH if not already added
if [[ ":$PATH:" != *":$DORADO_BIN:"* ]]; then
    echo "[INFO] Adding Dorado to PATH in $SHELL_PROFILE"
    {
        echo ""
        echo "# Added by Dorado installer"
        echo "export PATH=\"\$PATH:$DORADO_BIN\""
    } >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE" # Reload PATH for the current session
    echo "[INFO] PATH updated for the current session."
fi
# Ensure it's active for the rest of the script.
export PATH="$PATH:$DORADO_BIN"
echo ""

# === 2. kraken2_env ===
echo "[2/9] Creating environment: kraken2_env"
if ! conda env list | grep -q "kraken2_env"; then
    conda create -y -n kraken2_env -c bioconda -c conda-forge kraken2 seqtk || \
        { echo "[ERROR] Failed to create 'kraken2_env' Conda environment. Exiting."; exit 1; }
else
    echo "[INFO] Environment 'kraken2_env' already exists."
fi
echo "[INFO] Remember that Kraken2 requires a database. You will need to download it separately."
echo "Example: 'conda activate kraken2_env' then 'kraken2-build --download-library bacteria --db /path/to/your_kraken_db'"
echo "Or for a complete database: 'kraken2-build --download-taxonomy --db /path/to/your_kraken_db'"
echo "[INFO] 'kraken2_env' setup complete."
echo ""

# === 3. assembly_env ===
echo "[3/9] Creating environment: assembly_env"
if ! conda env list | grep -q "assembly_env"; then
    conda create -y -n assembly_env -c bioconda -c conda-forge \
        flye wtdbg paralleltask openjdk=8 || \
        { echo "[ERROR] Failed to create 'assembly_env' Conda environment. Exiting."; exit 1; }
else
    echo "[INFO] Environment 'assembly_env' already exists."
fi
echo "[INFO] 'assembly_env' setup complete."
echo ""

# === Canu install ===
echo "[INFO] Installing Canu v2.2..."
# Prompt for Canu installation directory
read -p "Path to the directory where Canu should be installed (e.g., $HOME/opt): " CANU_DIR
# Remove trailing slash if present
CANU_DIR="${CANU_DIR%/}"

CANU_VERSION="2.2"
CANU_ARCHIVE="canu-$CANU_VERSION.Linux-amd64.tar.xz"
CANU_URL="https://github.com/marbl/canu/releases/download/v$CANU_VERSION/$CANU_ARCHIVE"
CANU_BIN="$CANU_DIR/canu-$CANU_VERSION/bin"
CANU_DOWNLOAD_PATH="/tmp/$CANU_ARCHIVE" # Download to /tmp first

if [[ -x "$CANU_BIN/canu" ]]; then
    echo "[INFO] Canu is already installed in $CANU_BIN, skipping download."
else
    echo "[INSTALLATION] Attempting to install Canu v$CANU_VERSION..."
    mkdir -p "$CANU_DIR" || { echo "[ERROR] Failed to create directory $CANU_DIR for Canu. Exiting."; exit 1; }
    
    echo "[INFO] Downloading Canu from $CANU_URL to $CANU_DOWNLOAD_PATH..."
    if curl -fL "$CANU_URL" -o "$CANU_DOWNLOAD_PATH"; then
        if [[ -s "$CANU_DOWNLOAD_PATH" ]]; then
            echo "[INFO] Extracting Canu from $CANU_DOWNLOAD_PATH to $CANU_DIR..."
            if tar -xJf "$CANU_DOWNLOAD_PATH" -C "$CANU_DIR"; then
                rm -f "$CANU_DOWNLOAD_PATH" || echo "[WARNING] Failed to remove temporary Canu tarball: $CANU_DOWNLOAD_PATH"
                if [[ -x "$CANU_BIN/canu" ]]; then
                    echo "[INFO] Canu v$CANU_VERSION installed successfully."
                    echo "[INFO] Canu version check: $("$CANU_BIN/canu" --version 2>&1)"
                else
                    echo "[ERROR] Canu executable not found at $CANU_BIN/canu after extraction. Please check the extracted directory structure. Exiting.";
                    exit 1
                fi
            else
                echo "[ERROR] Failed to extract Canu from $CANU_DOWNLOAD_PATH. Exiting.";
                rm -f "$CANU_DOWNLOAD_PATH"
                exit 1
            fi
        else
            echo "[ERROR] Downloaded Canu file is empty or corrupted: $CANU_DOWNLOAD_PATH. Exiting.";
            rm -f "$CANU_DOWNLOAD_PATH"
            exit 1
        fi
    else
        echo "[ERROR] Failed to download Canu v$CANU_VERSION from $CANU_URL. Check network and URL. Exiting.";
        exit 1
    fi
fi

if [[ ":$PATH:" != *":$CANU_BIN:"* ]]; then
    echo "[INFO] Adding Canu to PATH in $SHELL_PROFILE"
    {
        echo ""
        echo "# Added by Canu installer"
        echo "export PATH=\"\$PATH:$CANU_BIN\""
    } >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE" # Reload PATH for the current session
    echo "[INFO] PATH updated for the current session."
fi
export PATH="$PATH:$CANU_BIN" # Ensure it's active for the rest of the script.
echo ""

# === Hifiasm install ===
echo "[INFO] Installing Hifiasm from source..."
# Prompt for Hifiasm installation directory
read -p "Path to the directory where Hifiasm should be installed (e.g., $HOME/opt): " HIFIASM_DIR
# Remove trailing slash if present
HIFIASM_DIR="${HIFIASM_DIR%/}"

HIFIASM_REPO="https://github.com/chhylp123/hifiasm.git"
HIFIASM_CLONE_DIR="$HIFIASM_DIR/hifiasm"
HIFIASM_BIN="$HIFIASM_CLONE_DIR/hifiasm"

if [[ -x "$HIFIASM_BIN" ]]; then
    echo "[INFO] Hifiasm is already installed in $HIFIASM_BIN, skipping build."
else
    echo "[INSTALLATION] Attempting to build Hifiasm from source..."
    mkdir -p "$HIFIASM_DIR" || { echo "[ERROR] Failed to create directory $HIFIASM_DIR for Hifiasm. Exiting."; exit 1; }
    
    if [[ -d "$HIFIASM_CLONE_DIR" ]]; then
        echo "[INFO] Hifiasm source directory already exists, pulling latest changes..."
        pushd "$HIFIASM_CLONE_DIR" > /dev/null || { echo "[ERROR] Failed to enter directory $HIFIASM_CLONE_DIR. Exiting."; exit 1; }
        git pull || { echo "[WARNING] Failed to pull latest Hifiasm changes. Proceeding with existing source."; }
        popd > /dev/null || { echo "[ERROR] Failed to return from Hifiasm source directory. Exiting."; exit 1; }
    else
        echo "[INFO] Cloning Hifiasm repository from $HIFIASM_REPO to $HIFIASM_CLONE_DIR..."
        git clone "$HIFIASM_REPO" "$HIFIASM_CLONE_DIR" || { echo "[ERROR] Failed to clone Hifiasm repository. Check git installation and network. Exiting."; exit 1; }
    fi

    pushd "$HIFIASM_CLONE_DIR" > /dev/null || { echo "[ERROR] Failed to enter Hifiasm source directory for compilation. Exiting."; exit 1; }
    echo "[INFO] Compiling Hifiasm. This might take a few minutes..."
    echo "[WARNING] Compiling Hifiasm requires 'make' and a C/C++ compiler (e.g., gcc)."
    echo "[WARNING] If you encounter errors, ensure 'build-essential' is installed on your system (e.g., sudo apt install build-essential)."
    make -j$(nproc) || { echo "[ERROR] Hifiasm compilation failed. Exiting."; popd > /dev/null; exit 1; }
    
    if [[ -x "$HIFIASM_BIN" ]]; then
        echo "[INFO] Hifiasm built successfully."
        echo "[INFO] Hifiasm version check: $("$HIFIASM_BIN" --version 2>&1)"
    else
        echo "[ERROR] Hifiasm executable not found at $HIFIASM_BIN after compilation. Exiting.";
        popd > /dev/null; exit 1
    fi
    popd > /dev/null || { echo "[ERROR] Failed to return from Hifiasm source directory. Exiting."; exit 1; }
fi

HIFIASM_PATH_DIR="$(dirname "$HIFIASM_BIN")"
if [[ ":$PATH:" != *":$HIFIASM_PATH_DIR:"* ]]; then
    echo "[INFO] Adding Hifiasm to PATH in $SHELL_PROFILE"
    {
        echo ""
        echo "# Added by Hifiasm installer"
        echo "export PATH=\"\$PATH:$HIFIASM_PATH_DIR\""
    } >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE" # Reload PATH for the current session
    echo "[INFO] PATH updated for the current session."
fi
export PATH="$PATH:$HIFIASM_PATH_DIR" # Ensure it's active for the rest of the script.
echo ""

# === 4. ndn_env (NextDenovo) ===
echo "[4/9] Creating environment: ndn_env (NextDenovo)"
# Updated to create environment with minimap2 and samtools directly
if ! conda env list | grep -q "ndn_env"; then
    echo "[INFO] Creating 'ndn_env' with NextDenovo dependencies (minimap2, samtools)..."
    conda create -y -n ndn_env python=3.10 -c bioconda -c conda-forge minimap2 samtools || \
        { echo "[ERROR] Failed to create 'ndn_env' Conda environment with dependencies. This might be due to conflicts with a global python pin (e.g., python=3.12). Exiting."; exit 1; }
else
    echo "[INFO] Environment 'ndn_env' already exists."
fi

# The following blocks for activating/installing dependencies are now redundant for minimap2, samtools
# as they are included in the conda create statement above.
# conda activate ndn_env || { echo "[ERROR] Failed to activate 'ndn_env'. Exiting."; exit 1; }
# echo "[INFO] Installing NextDenovo dependencies (minimap2, samtools) in 'ndn_env'..."
# conda install -y -c bioconda minimap2 samtools || \
#     { echo "[ERROR] Failed to install NextDenovo dependencies in 'ndn_env'. Exiting."; exit 1; }
# conda deactivate || { echo "[ERROR] Failed to deactivate 'ndn_env'. Exiting."; exit 1; }
echo "[INFO] NextDenovo dependencies setup handled during environment creation."

echo "[INFO] Cloning and setting up NextDenovo..."
# Prompt for NextDenovo installation directory
read -p "Path to the directory where NextDenovo should be installed (e.g., $HOME/opt): " NDN_DIR
# Remove trailing slash if present
NDN_DIR="${NDN_DIR%/}"

NEXTDENOVO_REPO="https://github.com/Nextomics/NextDenovo.git"
NEXTDENOVO_CLONE_DIR="$NDN_DIR/NextDenovo"
NDN_BIN="$NEXTDENOVO_CLONE_DIR/nextDenovo"

if [[ -f "$NDN_BIN" ]]; then # Check if the script file exists
    echo "[INFO] NextDenovo is already installed in $NDN_BIN, skipping clone."
else
    echo "[INSTALLATION] Attempting to clone NextDenovo..."
    mkdir -p "$NDN_DIR" || { echo "[ERROR] Failed to create directory $NDN_DIR for NextDenovo. Exiting."; exit 1; }
    
    if [[ -d "$NEXTDENOVO_CLONE_DIR" ]]; then
        echo "[INFO] NextDenovo source directory already exists, pulling latest changes..."
        pushd "$NEXTDENOVO_CLONE_DIR" > /dev/null || { echo "[ERROR] Failed to enter directory $NEXTDENOVO_CLONE_DIR. Exiting."; exit 1; }
        git pull || { echo "[WARNING] Failed to pull latest NextDenovo changes. Proceeding with existing source."; }
        popd > /dev/null || { echo "[ERROR] Failed to return from NextDenovo source directory. Exiting."; exit 1; }
    else
        echo "[INFO] Cloning NextDenovo repository from $NEXTDENOVO_REPO to $NEXTDENOVO_CLONE_DIR..."
        git clone "$NEXTDENOVO_REPO" "$NEXTDENOVO_CLONE_DIR" || { echo "[ERROR] Failed to clone NextDenovo repository. Check git installation and network. Exiting."; exit 1; }
    fi
    
    if [[ -f "$NDN_BIN" ]]; then
        echo "[INFO] NextDenovo cloned successfully."
    else
        echo "[ERROR] NextDenovo script not found at $NDN_BIN after cloning. Exiting.";
        exit 1
    fi
fi

NDN_LIB_PATH="$NEXTDENOVO_CLONE_DIR/lib"
if [[ ":$PYTHONPATH:" != *":$NDN_LIB_PATH:"* ]]; then
    echo "[INFO] Adding NextDenovo's paralleltask lib to PYTHONPATH in $SHELL_PROFILE"
    {
        echo ""
        echo "# Added by NextDenovo installer"
        echo "export PYTHONPATH=\"\$PYTHONPATH:$NDN_LIB_PATH\""
        echo "alias nextdenovo='python $NDN_BIN'"
    } >> "$SHELL_PROFILE"
    source "$SHELL_PROFILE" # Reload environment variables for the current session
    echo "[INFO] PYTHONPATH and 'nextdenovo' alias updated for the current session."
fi

# Apply changes to current session
export PYTHONPATH="$PYTHONPATH:$NDN_LIB_PATH"
# Note: Aliases cannot be directly exported for child processes.
# The `source` command above makes it available in the current shell.
# For new terminal sessions, the .bashrc change will apply.
alias nextdenovo="python $NDN_BIN"
echo "[INFO] 'ndn_env' and NextDenovo setup complete."
echo ""

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
        nanoplot fastqc qualimap bcftools samtools || \
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
        bwa pilon freebayes || \
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
