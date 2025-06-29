#!/usr/bin/env python3

import csv
from decimal import Decimal, InvalidOperation
import os
import sys

print("--- Genome Masking Script ---")
print("-----------------------------")

# ==== User Input ====
print("Please provide the required paths:")
tsv_file = input("Enter path to TSV file (e.g., RBH_all.tsv or RBH_100.tsv): ").strip()
fasta_file = input("Enter path to the input FASTA file (the contig to be masked): ").strip()
output_dir = input("Enter path to output directory for masked FASTA files: ").strip()

# ==== Input Validation ====
if not os.path.exists(tsv_file):
    sys.exit(f"[ERROR] TSV file not found: {tsv_file}")
if not os.path.exists(fasta_file):
    sys.exit(f"[ERROR] FASTA file not found: {fasta_file}")

# Extract prefix from FASTA file name (without extension)
prefix = os.path.basename(fasta_file)
prefix = os.path.splitext(prefix)[0]

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)
print(f"[INFO] Output directory: {output_dir}")

# ==== Helper Functions ====
def parse_tsv(tsv_path, min_id=Decimal("0.0"), max_id=Decimal("100.0")):
    regions = []
    try:
        with open(tsv_path, newline='') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            for row in reader:
                try:
                    identity = Decimal(row['Identity'].strip())
                    start_1b = int(row['Contig_Start'])
                    end_1b = int(row['Contig_End'])

                    start = min(start_1b, end_1b) - 1
                    end = max(start_1b, end_1b)

                    if min_id <= identity <= max_id:
                        regions.append((start, end))
                except (InvalidOperation, ValueError, KeyError) as e:
                    print(f"[WARNING] Skipping row due to parsing error: {row}. Error: {e}", file=sys.stderr)
    except Exception as e:
        sys.exit(f"[ERROR] An error occurred while parsing the TSV file: {e}")
    return regions


def load_fasta(fasta_path):
    header = ""
    sequence_lines = []
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header:
                        print("[WARNING] Multiple sequences found in FASTA. Only the first will be processed.", file=sys.stderr)
                        break
                    header = line[1:].strip()
                else:
                    sequence_lines.append(line)
        if not header or not sequence_lines:
            sys.exit(f"[ERROR] No valid FASTA sequence found in {fasta_path}")
    except Exception as e:
        sys.exit(f"[ERROR] An error occurred while loading the FASTA file: {e}")
    return header, ''.join(sequence_lines)


def mask_sequence(sequence, regions):
    masked = list(sequence)
    for start, end in regions:
        start = max(0, start)
        end = min(len(sequence), end)
        if start < end:
            masked[start:end] = ['N'] * (end - start)
    return "".join(masked)


def write_fasta(header, sequence, output_path, wrap_len=60):
    try:
        with open(output_path, 'w') as f:
            f.write(f">{header}\n")
            for i in range(0, len(sequence), wrap_len):
                f.write(sequence[i:i+wrap_len] + '\n')
    except Exception as e:
        print(f"[ERROR] Failed to write FASTA file {output_path}: {e}", file=sys.stderr)

# ==== MAIN EXECUTION ====
print(f"[INFO] Loading FASTA file: {fasta_file}")
header, sequence = load_fasta(fasta_file)
print(f"[INFO] FASTA sequence loaded. Length: {len(sequence)}")

# Mask 100% identity
print("[INFO] Masking regions with 100% identity...")
regions_100 = parse_tsv(tsv_file, min_id=Decimal("100.0"), max_id=Decimal("100.0"))
masked_100 = mask_sequence(sequence, regions_100)
output_100 = os.path.join(output_dir, f"{prefix}_masked_100.fasta")
write_fasta(header, masked_100, output_100)
print(f"[INFO] Saved: {output_100}")

# Mask 95â€“99.9% identity
print("[INFO] Masking regions with 95.0% to 99.9% identity...")
regions_95_999 = parse_tsv(tsv_file, min_id=Decimal("95.0"), max_id=Decimal("99.9"))
masked_95_999 = mask_sequence(sequence, regions_95_999)
output_95_999 = os.path.join(output_dir, f"{prefix}_masked_95_99.9.fasta")
write_fasta(header, masked_95_999, output_95_999)
print(f"[INFO] Saved: {output_95_999}")

# Mask ALL aligned regions
print("[INFO] Masking ALL aligned regions (regardless of identity)...")
regions_all = parse_tsv(tsv_file, min_id=Decimal("0.0"), max_id=Decimal("100.0"))
masked_all = mask_sequence(sequence, regions_all)
output_all = os.path.join(output_dir, f"{prefix}_masked_all.fasta")
write_fasta(header, masked_all, output_all)
print(f"[INFO] Saved: {output_all}")

print("\n-----------------------------")
print("Genome masking completed.")
