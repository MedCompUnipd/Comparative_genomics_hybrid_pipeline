#!/bin/bash
set -euo pipefail

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

# --- Load Conda Hook ---
echo "[INFO] Loading Conda hook..."
if command_exists conda; then
    eval "$(conda shell.bash hook)"
else
    echo "[ERROR] Conda not found. Please ensure Conda is installed and initialized in your shell profile."
    echo "Refer to the 'Requirements' section at the beginning of this script."
    exit 1
fi

echo ""
echo "--- INSTALLING SOFTWARE AND CONDA ENVIRONMENTS ---"
echo "------------------------------------------------"

# --- Environment: preprocessing_env (Script 1_dorado_bc_trim_int.sh, 4_hybrid.sh) ---
echo "[1/8] Creating/Activating environment: preprocessing_env"
if ! conda env list | grep -q "preprocessing_env"; then
    echo "[INFO] Creating 'preprocessing_env'..."
    conda create -y -n preprocessing_env -c bioconda -c conda-forge \
        python=3.9 ont-fast5-api=4.1.0 pod5=0.2.0 nanoq samtools minimap2 qualimap || { echo "[ERROR] Failed to create 'preprocessing_env' Conda environment. Exiting."; exit 1; }
else
    echo "[INFO] Environment 'preprocessing_env' already exists. Activating..."
fi
conda activate preprocessing_env || { echo "[ERROR] Failed to activate 'preprocessing_env' Conda environment. Exiting."; exit 1; }
echo "[INFO] Environment 'preprocessing_env' activated."

# Specific Dorado 0.9.6 installation
DORADO_VERSION="0.9.6"
DORADO_BIN_PATH=""

# Check if dorado is installed and is the correct version
if command_exists dorado; then
    CURRENT_DORADO_VERSION=$(dorado --version 2>&1 | grep -oP 'Version \K[0-9]+\.[0-9]+\.[0-9]+')
    if [[ "$CURRENT_DORADO_VERSION" == "$DORADO_VERSION" ]]; then
        echo "[INFO] Dorado v$DORADO_VERSION already installed."
        DORADO_BIN_PATH=$(dirname "$(command -v dorado)")
    else
        echo "[WARNING] Dorado is installed but not version $DORADO_VERSION (found $CURRENT_DORADO_VERSION)."
        read -p "Do you want to proceed with installing Dorado v$DORADO_VERSION? (y/n): " INSTALL_CHOICE
        if [[ "$INSTALL_CHOICE" =~ ^[Yy]$ ]]; then
            echo "[INSTALLATION] Installing Dorado v$DORADO_VERSION..."
            DORADO_INSTALL_DIR="$HOME/opt/dorado"
            mkdir -p "$DORADO_INSTALL_DIR" || { echo "[ERROR] Failed to create directory $DORADO_INSTALL_DIR for Dorado. Exiting."; exit 1; }
            if curl -L "https://cdn.oxfordnanopore.com/software/dorado/dorado-$DORADO_VERSION-linux-x64.tar.gz" | tar -xz -C "$DORADO_INSTALL_DIR"; then
                DORADO_EXTRACTED_DIR="$DORADO_INSTALL_DIR/dorado-$DORADO_VERSION-linux-x64"
                if [[ ! -d "$DORADO_EXTRACTED_DIR" ]]; then
                    echo "[ERROR] Dorado extracted directory not found at $DORADO_EXTRACTED_DIR. Tar extraction might have failed. Exiting."
                    exit 1
                fi
                DORADO_BIN_PATH="$DORADO_EXTRACTED_DIR/bin"
                echo "export PATH=\"$DORADO_BIN_PATH:\$PATH\"" >> "$SHELL_PROFILE"
                source "$SHELL_PROFILE" # Reload PATH for current session
                echo "[INFO] Dorado v$DORADO_VERSION installed in $DORADO_BIN_PATH and added to PATH."
            else
                echo "[ERROR] Failed to download or extract Dorado v$DORADO_VERSION. Check network and URL. Exiting."
                exit 1
            fi
        else
            echo "[INFO] Skipping Dorado v$DORADO_VERSION installation."
        fi
    fi
else
    echo "[INSTALLATION] Installing Dorado v$DORADO_VERSION..."
    DORADO_INSTALL_DIR="$HOME/opt/dorado"
    mkdir -p "$DORADO_INSTALL_DIR" || { echo "[ERROR] Failed to create directory $DORADO_INSTALL_DIR for Dorado. Exiting."; exit 1; }
    if curl -L "https://cdn.oxfordnanopore.com/software/dorado/dorado-$DORADO_VERSION-linux-x64.tar.gz" | tar -xz -C "$DORADO_INSTALL_DIR"; then
        DORADO_EXTRACTED_DIR="$DORADO_INSTALL_DIR/dorado-$DORADO_VERSION-linux-x64"
        if [[ ! -d "$DORADO_EXTRACTED_DIR" ]]; then
            echo "[ERROR] Dorado extracted directory not found at $DORADO_EXTRACTED_DIR. Tar extraction might have failed. Exiting."
            exit 1
        fi
        DORADO_BIN_PATH="$DORADO_EXTRACTED_DIR/bin"
        echo "export PATH=\"$DORADO_BIN_PATH:\$PATH\"" >> "$SHELL_PROFILE"
        source "$SHELL_PROFILE" # Reload PATH for current session
        echo "[INFO] Dorado v$DORADO_VERSION installed in $DORADO_BIN_PATH and added to PATH."
    else
        echo "[ERROR] Failed to download or extract Dorado v$DORADO_VERSION. Check network and URL. Exiting."
        exit 1
    fi
fi
conda deactivate || { echo "[ERROR] Failed to deactivate 'preprocessing_env'. Exiting."; exit 1; }
echo ""

# --- Environment: kraken_env (Script 2_tax_sel_qc.sh) ---
echo "[2/8] Creating/Activating environment: kraken_env"
if ! conda env list | grep -q "kraken_env"; then
    echo "[INFO] Creating 'kraken_env'..."
    conda create -y -n kraken_env -c bioconda -c conda-forge \
        kraken2 bracken seqtk || { echo "[ERROR] Failed to create 'kraken_env' Conda environment. Exiting."; exit 1; }
else
    echo "Environment 'kraken_env' already exists. Activating..."
fi
conda activate kraken_env || { echo "[ERROR] Failed to activate 'kraken_env' Conda environment. Exiting."; exit 1; }
echo "[INFO] Environment 'kraken_env' activated."
echo "[INFO] Remember that Kraken2 requires a database. You will need to download it separately."
echo "Example: kraken2-build --download-library bacteria --db /path/to/kraken_db"
echo "Or for a complete database: kraken2-build --download-taxonomy --db /path/to/kraken_db"
conda deactivate || { echo "[ERROR] Failed to deactivate 'kraken_env'. Exiting."; exit 1; }
echo ""

# --- Environment: qc_env (Script 2_tax_sel_qc.sh, 3_asm_pol.sh, 5_comparation.sh) ---
echo "[3/8] Creating/Activating environment: qc_env"
if ! conda env list | grep -q "qc_env"; then
    echo "[INFO] Creating 'qc_env'..."
    conda create -y -n qc_env -c bioconda -c conda-forge \
        nanoplot=1.40.0 quast bcftools htslib samtools mummer || { echo "[ERROR] Failed to create 'qc_env' Conda environment. Exiting."; exit 1; }
else
    echo "Environment 'qc_env' already exists. Activating..."
fi
conda activate qc_env || { echo "[ERROR] Failed to activate 'qc_env' Conda environment. Exiting."; exit 1; }
echo "[INFO] Environment 'qc_env' activated."
conda deactivate || { echo "[ERROR] Failed to deactivate 'qc_env'. Exiting."; exit 1; }
echo ""

# --- Environment: assembly_env (Script 3_asm_pol.sh) ---
echo "[4/8] Creating/Activating environment: assembly_env"
if ! conda env list | grep -q "assembly_env"; then
    echo "[INFO] Creating 'assembly_env'..."
    conda create -y -n assembly_env -c bioconda -c conda-forge \
        flye wtdbg openjdk=8 paralleltask || { echo "[ERROR] Failed to create 'assembly_env' Conda environment. Exiting."; exit 1; }
else
    echo "Environment 'assembly_env' already exists. Activating..."
fi
conda activate assembly_env || { echo "[ERROR] Failed to activate 'assembly_env' Conda environment. Exiting."; exit 1; }
echo "[INFO] Environment 'assembly_env' activated."

# Specific Canu 2.2 installation
CANU_VERSION="2.2"
CANU_BIN_PATH=""

if command_exists canu; then
    CURRENT_CANU_VERSION=$(canu --version 2>&1 | grep -oP 'v\K[0-9]+\.[0-9]+')
    if [[ "$CURRENT_CANU_VERSION" == "$CANU_VERSION" ]]; then
        echo "[INFO] Canu v$CANU_VERSION already installed."
        CANU_BIN_PATH=$(dirname "$(command -v canu)")
    else
        echo "[WARNING] Canu is installed but not version $CANU_VERSION (found $CURRENT_CANU_VERSION)."
        read -p "Do you want to proceed with installing Canu v$CANU_VERSION? (y/n): " INSTALL_CHOICE
        if [[ "$INSTALL_CHOICE" =~ ^[Yy]$ ]]; then
            echo "[INSTALLATION] Downloading and installing Canu v$CANU_VERSION..."
            CANU_INSTALL_DIR="$HOME/opt/canu"
            mkdir -p "$CANU_INSTALL_DIR" || { echo "[ERROR] Failed to create directory $CANU_INSTALL_DIR for Canu. Exiting."; exit 1; }
            if curl -L "https://github.com/marbl/canu/releases/download/v$CANU_VERSION/canu-$CANU_VERSION.tar.gz" | tar -xz -C "$CANU_INSTALL_DIR"; then
                CANU_EXTRACTED_DIR="$CANU_INSTALL_DIR/canu-$CANU_VERSION"
                if [[ ! -d "$CANU_EXTRACTED_DIR" ]]; then
                    echo "[ERROR] Canu extracted directory not found at $CANU_EXTRACTED_DIR. Tar extraction might have failed. Exiting."
                    exit 1
                fi
                CANU_BIN_PATH="$CANU_EXTRACTED_DIR/bin"
                echo "export PATH=\"$CANU_BIN_PATH:\$PATH\"" >> "$SHELL_PROFILE"
                source "$SHELL_PROFILE" # Reload PATH for current session
                echo "[INFO] Canu v$CANU_VERSION installed in $CANU_BIN_PATH and added to PATH."
            else
                echo "[ERROR] Failed to download or extract Canu v$CANU_VERSION. Check network and URL. Exiting."
                exit 1
            fi
        else
            echo "[INFO] Skipping Canu v$CANU_VERSION installation."
        fi
    fi
else
    echo "[INSTALLATION] Downloading and installing Canu v$CANU_VERSION..."
    CANU_INSTALL_DIR="$HOME/opt/canu"
    mkdir -p "$CANU_INSTALL_DIR" || { echo "[ERROR] Failed to create directory $CANU_INSTALL_DIR for Canu. Exiting."; exit 1; }
    if curl -L "https://github.com/marbl/canu/releases/download/v$CANU_VERSION/canu-$CANU_VERSION.tar.gz" | tar -xz -C "$CANU_INSTALL_DIR"; then
        CANU_EXTRACTED_DIR="$CANU_INSTALL_DIR/canu-$CANU_VERSION"
        if [[ ! -d "$CANU_EXTRACTED_DIR" ]]; then
            echo "[ERROR] Canu extracted directory not found at $CANU_EXTRACTED_DIR. Tar extraction might have failed. Exiting."
            exit 1
        fi
        CANU_BIN_PATH="$CANU_EXTRACTED_DIR/bin"
        echo "export PATH=\"$CANU_BIN_PATH:\$PATH\"" >> "$SHELL_PROFILE"
        source "$SHELL_PROFILE" # Reload PATH for current session
        echo "[INFO] Canu v$CANU_VERSION installed in $CANU_BIN_PATH and added to PATH."
    else
        echo "[ERROR] Failed to download or extract Canu v$CANU_VERSION. Check network and URL. Exiting."
        exit 1
    fi
fi
conda deactivate || { echo "[ERROR] Failed to deactivate 'assembly_env'. Exiting."; exit 1; }
echo ""

# --- Environment: polishing_env (Script 3_asm_pol.sh, 4_hybrid.sh) ---
echo "[5/8] Creating/Activating environment: polishing_env"
if ! conda env list | grep -q "polishing_env"; then
    echo "[INFO] Creating 'polishing_env'..."
    conda create -y -n polishing_env -c bioconda -c conda-forge \
        racon medaka bwa pilon || { echo "[ERROR] Failed to create 'polishing_env' Conda environment. Exiting."; exit 1; }
else
    echo "Environment 'polishing_env' already exists. Activating..."
fi
conda activate polishing_env || { echo "[ERROR] Failed to activate 'polishing_env' Conda environment. Exiting."; exit 1; }
echo "[INFO] Environment 'polishing_env' activated."
conda deactivate || { echo "[ERROR] Failed to deactivate 'polishing_env'. Exiting."; exit 1; }
echo ""

# --- Installing Hifiasm (Manual Download) ---
echo "[INSTALLATION] Hifiasm"
HIFIASM_VERSION="0.19.8" # Verify the latest stable version
HIFIASM_INSTALL_DIR="$HOME/opt/hifiasm"
HIFIASM_BINARY_CONTENTS_DIR="$HIFIASM_INSTALL_DIR/hifiasm-$HIFIASM_VERSION-Linux" # The directory containing the hifiasm executable after extraction

if ! command_exists hifiasm || [[ "$(hifiasm --version 2>&1 | grep -oP 'hifiasm v\K[0-9]+\.[0-9]+\.[0-9]+')" != "$HIFIASM_VERSION" ]]; then
    echo "[INSTALLATION] Hifiasm v$HIFIASM_VERSION needs to be downloaded and installed."
    echo "[INFO] Downloading and extracting Hifiasm to $HIFIASM_INSTALL_DIR..."
    mkdir -p "$HIFIASM_INSTALL_DIR" || { echo "[ERROR] Failed to create directory $HIFIASM_INSTALL_DIR for Hifiasm. Exiting."; exit 1; }
    if curl -L "https://github.com/chhylp123/hifiasm/releases/download/$HIFIASM_VERSION/hifiasm-$HIFIASM_VERSION-Linux.tar.gz" | tar -xz -C "$HIFIASM_INSTALL_DIR"; then
        echo "[INFO] Hifiasm downloaded and extracted successfully."
        if [[ ! -d "$HIFIASM_BINARY_CONTENTS_DIR" ]]; then
            echo "[ERROR] Extracted Hifiasm directory not found at $HIFIASM_BINARY_CONTENTS_DIR. Tar extraction might have failed. Exiting."
            exit 1
        fi
        echo "export PATH=\"$HIFIASM_BINARY_CONTENTS_DIR:\$PATH\"" >> "$SHELL_PROFILE"
        source "$SHELL_PROFILE" # Reload PATH for current session
        echo "[INFO] Hifiasm v$HIFIASM_VERSION installed and added to PATH. Please restart your terminal or 'source $SHELL_PROFILE' for changes to take full effect."
    else
        echo "[ERROR] Failed to download or extract Hifiasm v$HIFIASM_VERSION. Check network and URL. Exiting."
        exit 1
    fi
else
    echo "[INFO] Hifiasm v$HIFIASM_VERSION appears to be already installed and is the correct version."
fi
echo ""

# --- Installing NextDenovo (Manual Download) ---
echo "[INSTALLATION] NextDenovo"
NEXTDENOVO_VERSION="2.5.0" # Verify the latest stable version
NEXTDENOVO_INSTALL_DIR="$HOME/opt/nextdenovo"
NEXTDENOVO_EXTRACTED_ROOT_DIR="$NEXTDENOVO_INSTALL_DIR/NextDenovo" # The root directory after unzipping
NEXTDENOVO_BIN_PATH="$NEXTDENOVO_EXTRACTED_ROOT_DIR/bin" # The directory containing NextDenovo executable

if ! command_exists NextDenovo; then
    echo "[INSTALLATION] NextDenovo v$NEXTDENOVO_VERSION needs to be downloaded and installed."
    echo "[INFO] Downloading and extracting NextDenovo to $NEXTDENOVO_INSTALL_DIR..."
    mkdir -p "$NEXTDENOVO_INSTALL_DIR" || { echo "[ERROR] Failed to create directory $NEXTDENOVO_INSTALL_DIR for NextDenovo. Exiting."; exit 1; }
    # Download to /tmp first, then unzip
    TEMP_NEXTDENOVO_ZIP="/tmp/NextDenovo-$NEXTDENOVO_VERSION.zip"
    if curl -L "https://github.com/NextDenovo/NextDenovo/releases/download/v$NEXTDENOVO_VERSION/NextDenovo.zip" -o "$TEMP_NEXTDENOVO_ZIP"; then
        if unzip -q "$TEMP_NEXTDENOVO_ZIP" -d "$NEXTDENOVO_INSTALL_DIR"; then # -q for quiet unzip
            echo "[INFO] NextDenovo downloaded and extracted successfully."
            if [[ ! -d "$NEXTDENOVO_EXTRACTED_ROOT_DIR" ]]; then
                echo "[ERROR] Extracted NextDenovo root directory not found. Expected at $NEXTDENOVO_EXTRACTED_ROOT_DIR. Unzip might have failed. Exiting."
                exit 1
            fi
            echo "export PATH=\"$NEXTDENOVO_BIN_PATH:\$PATH\"" >> "$SHELL_PROFILE"
            source "$SHELL_PROFILE" # Reload PATH for current session
            echo "[INFO] NextDenovo v$NEXTDENOVO_VERSION installed and added to PATH. Please restart your terminal or 'source $SHELL_PROFILE' for changes to take full effect."
        else
            echo "[ERROR] Failed to unzip NextDenovo v$NEXTDENOVO_VERSION. Check zip file integrity. Exiting."
            exit 1
        fi
        rm "$TEMP_NEXTDENOVO_ZIP" # Clean up downloaded zip
    else
        echo "[ERROR] Failed to download NextDenovo v$NEXTDENOVO_VERSION. Check network and URL. Exiting."
        exit 1
    fi
else
    echo "[INFO] NextDenovo appears to be already installed."
fi
echo ""

# --- Environment: ndn_env (Script 3_asm_pol.sh) - Likely redundant if NextDenovo is standalone ---
echo "[6/8] Creating/Activating environment: ndn_env (for script compatibility, if NextDenovo has no Conda dependencies)"
if ! conda env list | grep -q "ndn_env"; then
    echo "[INFO] Creating 'ndn_env'..."
    conda create -y -n ndn_env python=3.9 || { echo "[ERROR] Failed to create 'ndn_env' Conda environment. Exiting."; exit 1; }
else
    echo "Environment 'ndn_env' already exists. Activating..."
fi
conda activate ndn_env || { echo "[ERROR] Failed to activate 'ndn_env' Conda environment. Exiting."; exit 1; }
echo "[INFO] Environment 'ndn_env' activated."
conda deactivate || { echo "[ERROR] Failed to deactivate 'ndn_env'. Exiting."; exit 1; }
echo ""

# --- Environment: canu_env (Script 3_asm_pol.sh) ---
echo "[7/8] Creating/Activating environment: canu_env (for script compatibility, if Canu has no Conda dependencies)"
if ! conda env list | grep -q "canu_env"; then
    echo "[INFO] Creating 'canu_env'..."
    conda create -y -n canu_env python=3.9 || { echo "[ERROR] Failed to create 'canu_env' Conda environment. Exiting."; exit 1; }
else
    echo "Environment 'canu_env' already exists. Activating..."
fi
conda activate canu_env || { echo "[ERROR] Failed to activate 'canu_env' Conda environment. Exiting."; exit 1; }
echo "[INFO] Environment 'canu_env' activated."
conda deactivate || { echo "[ERROR] Failed to deactivate 'canu_env'. Exiting."; exit 1; }
echo ""

# --- Environment: fasta3_env (Script 10_extract_align.sh, 11_nuctoaa.sh) ---
echo "[8/8] Creating/Activating environment: fasta3_env"
if ! conda env list | grep -q "fasta3_env"; then
    echo "[INFO] Creating 'fasta3_env'..."
    conda create -y -n fasta3_env python=3.9 || { echo "[ERROR] Failed to create 'fasta3_env' Conda environment. Exiting."; exit 1; }
else
    echo "Environment 'fasta3_env' already exists. Activating..."
fi
conda activate fasta3_env || { echo "[ERROR] Failed to activate 'fasta3_env' Conda environment. Exiting."; exit 1; }
echo "[INFO] Environment 'fasta3_env' activated."
conda deactivate || { echo "[ERROR] Failed to deactivate 'fasta3_env'. Exiting."; exit 1; }
echo ""


# --- Installing additional Python packages (for 7_rbh_lists.py, 8_extraction_fastagenes.py, 9_masking.py) ---
echo ""
echo "--- Installing Additional Python Packages ---"
echo "Installing pandas and Biopython into the 'base' Conda environment (or an environment of your choice)."
echo "These are required for Python scripts 7_rbh_lists.py, 8_extraction_fastagenes.py, and 9_masking.py."

read -p "Do you want to install pandas and Biopython into the 'base' Conda environment? (y/n): " INSTALL_PYTHON_LIBS
if [[ "$INSTALL_PYTHON_LIBS" =~ ^[Yy]$ ]]; then
    conda activate base || { echo "[ERROR] Failed to activate 'base' environment for Python packages. Exiting."; exit 1; }
    pip install pandas biopython || { echo "[ERROR] Failed to install pandas or biopython. Check network and pip. Exiting."; exit 1; }
    conda deactivate || { echo "[ERROR] Failed to deactivate 'base' environment after Python package installation. Exiting."; exit 1; }
    echo "[INFO] pandas and biopython installed in the 'base' environment."
else
    echo "[INFO] Skipping installation of pandas and biopython. Ensure they are available in the environment where you run your Python scripts."
fi
echo ""

# --- Installing FASTA package (for glsearch36 in 6_rbh.sh and 10_extract_align.sh, and fasty36 in 11_nuctoaa.sh) ---
echo "[INSTALLATION] FASTA package (for glsearch36 and fasty36)"
if ! command_exists glsearch36 || ! command_exists fasty36; then
    read -p "Do you want to install the FASTA package (which includes glsearch36 and fasty36)? (y/n): " INSTALL_FASTA_PKG
    if [[ "$INSTALL_FASTA_PKG" =~ ^[Yy]$ ]]; then
        echo "[WARNING] Compiling FASTA package requires 'make' and a C/C++ compiler (e.g., gcc)."
        echo "[WARNING] If you encounter errors, ensure 'build-essential' is installed on your system (e.g., sudo apt install build-essential)."
        echo "[INFO] Downloading FASTA package (v36.3.8g, recommended version)..."
        FASTA_VERSION="36.3.8g"
        FASTA_DIR="$HOME/opt/fasta_package"
        mkdir -p "$FASTA_DIR" || { echo "[ERROR] Failed to create directory $FASTA_DIR for FASTA package. Exiting."; exit 1; }

        # Use curl to download and pipe to tar for extraction
        if curl -L "https://fasta.bioch.virginia.edu/fasta_tf/fasta-$FASTA_VERSION.tar.gz" | tar -xz -C "$FASTA_DIR"; then
            echo "[INFO] Compiling FASTA package. This might take a few minutes..."
            FASTA_SRC_DIR="$FASTA_DIR/fasta-$FASTA_VERSION/src"
            if [[ ! -d "$FASTA_SRC_DIR" ]]; then
                echo "[ERROR] FASTA package source directory not found at $FASTA_SRC_DIR. Tar extraction might have failed. Exiting."
                exit 1
            fi
            pushd "$FASTA_SRC_DIR" > /dev/null || { echo "[ERROR] Failed to enter FASTA src directory for compilation. Exiting."; exit 1; }
            make -f ../make/Makefile.linux64 || { echo "[ERROR] FASTA compilation failed. Please ensure 'make' and a C/C++ compiler (like 'gcc') are installed. Try installing 'build-essential' on Debian/Ubuntu systems. Exiting."; popd > /dev/null; exit 1; }
            popd > /dev/null || { echo "[ERROR] Failed to return from FASTA src directory. Exiting."; exit 1; }

            FASTA_BIN_PATH="$FASTA_DIR/fasta-$FASTA_VERSION/bin"
            echo "export PATH=\"$FASTA_BIN_PATH:\$PATH\"" >> "$SHELL_PROFILE"
            source "$SHELL_PROFILE" # Reload PATH for current session
            echo "[INFO] glsearch36 and fasty36 binaries installed and added to PATH."
        else
            echo "[ERROR] Failed to download or extract FASTA package v$FASTA_VERSION. Check network and URL. Exiting."
            exit 1
        fi
    else
        echo "[WARNING] glsearch36 and fasty36 will not be installed automatically. You will need to install them manually to run 6_rbh.sh, 10_extract_align.sh and 11_nuctoaa.sh."
    fi
else
    echo "[INFO] glsearch36 and fasty36 appear to be already installed."
fi
echo ""

echo "#####################################################"
echo "#           SETUP COMPLETE!                         #"
echo "#####################################################"
echo "For manual installations (Hifiasm, NextDenovo, and FASTA package if not automatically compiled),"
echo "ensure that their binaries are in your PATH. You might need to restart your terminal or"
echo "run 'source $SHELL_PROFILE' to ensure all PATH changes are active."
echo ""
echo "Good luck with your work!"
