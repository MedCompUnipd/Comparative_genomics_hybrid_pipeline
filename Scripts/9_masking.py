#!/usr/bin/env python3

import os
import sys
import csv
from Bio import SeqIO

print("\n--- FASTA Genes Extraction Script (RBH TSV) ---")
print("------------------------------------------------")

# ==== Input ====
fasta_path = input("Enter path to the de novo genome FASTA file (e.g., your_contig.fasta): ").strip()
tsv_path = input("Enter path to the TSV file (e.g., RBH_all.tsv from script 7): ").strip()
output_path = input("Enter desired output FASTA file path for extracted genes: ").strip()

# ==== Validate ====
if not os.path.exists(fasta_path):
    sys.exit(f"[ERROR] FASTA file not found: {fasta_path}")
if not os.path.exists(tsv_path):
    sys.exit(f"[ERROR] TSV file not found: {tsv_path}")

# Ensure output directory exists
output_dir = os.path.dirname(output_path)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)
    print(f"[INFO] Created output directory: {output_dir}")

# ==== Load contig sequence ====
print(f"[INFO] Loading contig sequence from: {fasta_path}")
try:
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        sys.exit("[ERROR] No sequences found in the FASTA file.")
    if len(records) > 1:
        print("[WARNING] Multiple sequences found. Using the first.")
    contig_seq = str(records[0].seq)
except Exception as e:
    sys.exit(f"[ERROR] Failed to load FASTA: {e}")

seq_len = len(contig_seq)
print(f"[INFO] Contig sequence loaded, length = {seq_len} bp.")

# ==== Extract sequences ====
written_count = 0
try:
    with open(tsv_path, newline="") as tsv_file, open(output_path, "w") as out_fasta:
        reader = csv.DictReader(tsv_file, delimiter="\t")
        for row in reader:
            try:
                gene_id = row["Gene"]
                start = int(row["Contig_Start"])
                end = int(row["Contig_End"])

                # Convert to 0-based indexing for slicing
                start_0 = min(start, end) - 1
                end_0 = max(start, end)

                if start_0 < 0 or end_0 > seq_len:
                    print(f"[WARNING] Skipping {gene_id}: out-of-bounds ({start}-{end})", file=sys.stderr)
                    continue

                seq = contig_seq[start_0:end_0]
                header = f">{gene_id}_ContigRegion_{start}-{end}"
                out_fasta.write(header + "\n")

                # Wrap sequence 60 characters per line
                for i in range(0, len(seq), 60):
                    out_fasta.write(seq[i:i+60] + "\n")

                written_count += 1

            except (KeyError, ValueError) as e:
                print(f"[WARNING] Skipping row: {row}. Error: {e}", file=sys.stderr)
                continue

except Exception as e:
    sys.exit(f"[ERROR] Failed during extraction: {e}")

print(f"[INFO] Extracted {written_count} gene regions.")
print(f"[INFO] Output written to: {output_path}")
print("------------------------------------------------")
