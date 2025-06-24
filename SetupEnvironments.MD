# Installation Setup Script

This script will help you set up your working environment by installing all necessary software for your bioinformatics pipeline. It uses conda for environment and package management where possible and provides instructions for manual installations for tools that require them.

Requirements:

Miniconda3 or Anaconda3 must already be installed and initialized. If not, please follow these instructions before running this script:
Download Miniconda: wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
Run the installer: bash Miniconda3-latest-Linux-x86_64.sh (accept terms, choose path, and answer yes to initialization).
Restart your terminal or run source ~/.bashrc (or source ~/.zshrc if you use zsh).
Basic build tools (make, gcc/g++) should be installed on your system for compiling tools like FASTA package, if you choose to install it. If not, you might need to install them manually (e.g., sudo apt install build-essential on Ubuntu/Debian).
This script assumes a Debian/Ubuntu-based Linux distribution. Adjust PATH modifications for other shells (e.g., ~/.zshrc for zsh) if needed.
# How to run the script:

Save the following content to a file named setup_envs.sh.
Make the script executable: chmod +x setup_envs.sh
Execute the script: ./setup_envs.sh
