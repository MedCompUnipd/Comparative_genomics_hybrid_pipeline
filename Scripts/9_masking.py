#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from decimal import Decimal, InvalidOperation
import os
import sys
import io # Import io module for StringIO

print("--- Genome Masking Script ---")
print("-----------------------------")

# ==== User Input ====
print("Please provide the required paths:")
tsv_file = input("Enter path to TSV alignment mummer file (e.g., ref_denovo.tsv): ").strip()
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
            lines = tsvfile.readlines()
            
            # Check if there are enough lines in the file
            # Now expecting at least 4 lines (3 to skip + 1 header line)
            if len(lines) < 4: 
                sys.exit(f"[ERROR] TSV file '{tsv_path}' does not contain enough lines (expected at least 4 including header).")

            # The header is on the 4th line (index 3)
            header_line = lines[3].strip()
            
            # --- NEW: Strip brackets from header_line ---
            header_line = header_line.replace('[', '').replace(']', '')
            # --- END NEW ---
            
            # The data starts from the 5th line (index 4) onwards
            data_lines = [line.strip() for line in lines[4:] if line.strip()] # Filter out empty lines

            if not header_line:
                sys.exit(f"[ERROR] The header line (4th line) in TSV file '{tsv_path}' is empty.")
            
            # Combine header and data lines into a single string to be read by DictReader
            full_tsv_content = header_line + '\n' + '\n'.join(data_lines)
            
            # Use io.StringIO to treat the string as a file-like object
            tsv_reader_input = io.StringIO(full_tsv_content)
            
            # Utilizziamo le intestazioni del tuo file TSV
            # Le colonne di interesse sono '% IDY', 'S1' e 'E1'
            reader = csv.DictReader(tsv_reader_input, delimiter='\t') 
            
            # Controllo per assicurarsi che le colonne necessarie siano presenti
            required_columns = ['% IDY', 'S1', 'E1']
            if not all(col in reader.fieldnames for col in required_columns):
                sys.exit(f"[ERROR] TSV file is missing one or more required columns: {', '.join(required_columns)}. Found: {', '.join(reader.fieldnames)}")

            for row in reader:
                try:
                    # Modificato per usare i nomi delle colonne dal tuo TSV
                    identity = Decimal(row['% IDY'].strip())
                    start_1b = int(row['S1'])
                    end_1b = int(row['E1'])

                    # Mantiene la logica per gestire allineamenti in entrambe le direzioni
                    # Converti a 0-based start (Python slicing � esclusivo sull'end)
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


# Mask ALL aligned regions
print("[INFO] Masking ALL aligned regions (regardless of identity)...")
regions_all = parse_tsv(tsv_file, min_id=Decimal("0.0"), max_id=Decimal("100.0"))
masked_all = mask_sequence(sequence, regions_all)
output_all = os.path.join(output_dir, f"{prefix}_masked_all.fasta")
write_fasta(header, masked_all, output_all)
print(f"[INFO] Saved: {output_all}")

print("\n-----------------------------")
print("Genome masking completed.")
