#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO

# ==== Conda Activation (optional, handled by calling shell script) ====
# Python scripts typically rely on the calling shell to activate the Conda environment.
# Direct activation via subprocess.run is generally discouraged within Python for general use,
# as it can lead to environment conflicts.
# The original script had `activate_conda_env`, but it's commented out here
# because the calling shell script `10_extract_align.sh` already activates `fasta3_env`.

print("--- FASTA Genes Extraction Script ---")
print("-----------------------------------")

# ==== User Input ====
fasta_path = input("Enter path to the de novo genome FASTA file (e.g., your_contig.fasta): ").strip()
tsv_path = input("Enter path to the TSV file (e.g., RBH_all.tsv from script 7): ").strip()
output_path = input("Enter desired output FASTA file path for extracted genes: ").strip()

# ==== Input Validation ====
if not os.path.exists(fasta_path):
    sys.exit(f"[ERROR] FASTA file not found: {fasta_path}")
if not os.path.exists(tsv_path):
    sys.exit(f"[ERROR] TSV file not found: {tsv_path}")

# Ensure output directory exists
output_dir = os.path.dirname(output_path)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    print(f"[INFO] Created output directory: {output_dir}")

# ==== Load Contig Sequence ====
print(f"[INFO] Loading contig sequence from: {fasta_path}")
contig_seq = ""
try:
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        sys.exit("[ERROR] No sequences found in the FASTA file.")
    if len(records) > 1:
        print("[WARNING] Multiple sequences found in FASTA file. Using the first sequence as the contig.")
    contig_seq = str(records[0].seq)
except Exception as e:
    sys.exit(f"[ERROR] An error occurred while loading the FASTA file: {e}")

seq_len = len(contig_seq)
print(f"[INFO] Loaded contig sequence of length {seq_len}.")

# ==== Extract Genes based on TSV ====
print(f"[INFO] Extracting gene regions based on TSV file: {tsv_path}")
import csv # Import csv here, as it's used later

written_count = 0
try:
    with open(tsv_path, newline="") as tsv_file, open(output_path, "w") as out_fasta:
        reader = csv.DictReader(tsv_file, delimiter="\t")
        for row in reader:
            try:
                gene_id = row["Gene"]
                # Convert 1-based coordinates from TSV to 0-based for Python slicing
                start_1_based = int(row["Contig_start"])
                end_1_based = int(row["Contig_end"])

                # Ensure start is always less than or equal to end for slicing
                start_0_based = min(start_1_based, end_1_based) - 1
                end_0_based = max(start_1_based, end_1_based)

            except (KeyError, ValueError) as e:
                print(f"[WARNING] Skipping invalid row due to parsing error: {row}. Error: {e}", file=sys.stderr)
                continue

            # Validate coordinates against sequence length
            if start_0_based < 0 or end_0_based > seq_len or start_0_based >= end_0_based:
                print(f"[WARNING] Coordinates out of bounds or invalid for gene {gene_id}: {start_1_based}-{end_1_based} (sequence length: {seq_len}). Skipping.", file=sys.stderr)
                continue

            extracted_seq = contig_seq[start_0_based:end_0_based]

            # Write to output FASTA
            out_fasta.write(f">{gene_id}_ContigRegion_{start_1_based}-{end_1_based}\n")
            # Wrap sequence if desired (e.g., 60 characters per line)
            for i in range(0, len(extracted_seq), 60):
                out_fasta.write(extracted_seq[i:i+60] + "\n")
            written_count += 1

except Exception as e:
    sys.exit(f"[ERROR] An error occurred during gene extraction: {e}")

print(f"[INFO] Successfully extracted {written_count} gene regions to: {output_path}")
print("\n-----------------------------------")
print("FASTA gene extraction completed.")
