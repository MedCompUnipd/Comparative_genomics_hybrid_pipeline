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
    """
    Parses a TSV file to extract sequence regions based on identity thresholds.

    Args:
        tsv_path (str): Path to the TSV file.
        min_id (Decimal): Minimum identity threshold (inclusive).
        max_id (Decimal): Maximum identity threshold (inclusive).

    Returns:
        list: A list of tuples, where each tuple contains (start, end) of a region (0-based).
    """
    regions = []
    try:
        with open(tsv_path, newline='') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            for row in reader:
                try:
                    identity = Decimal(row['Identity'].strip())
                    # Convert 1-based coordinates from TSV to 0-based for Python slicing
                    start_1_based = int(row['Contig_start'])
                    end_1_based = int(row['Contig_end'])

                    # Ensure start is always less than or equal to end for slicing
                    start_0_based = min(start_1_based, end_1_based) - 1
                    end_0_based = max(start_1_based, end_1_based)

                    if min_id <= identity <= max_id:
                        regions.append((start_0_based, end_0_based))
                except (InvalidOperation, ValueError, KeyError) as e:
                    print(f"[WARNING] Skipping row due to parsing error: {row}. Error: {e}", file=sys.stderr)
                    continue
    except FileNotFoundError:
        sys.exit(f"[ERROR] TSV file not found: {tsv_path}")
    except Exception as e:
        sys.exit(f"[ERROR] An error occurred while parsing the TSV file: {e}")
    return regions


def load_fasta(fasta_path):
    """
    Loads a single FASTA sequence from a file.

    Args:
        fasta_path (str): Path to the FASTA file.

    Returns:
        tuple: (header, sequence_string)
    """
    header = ""
    sequence_lines = []
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header: # Already read a sequence, assuming single sequence file
                        print("[WARNING] Multiple sequences found in FASTA. Only the first will be processed.", file=sys.stderr)
                        break
                    header = line[1:].strip() # Remove '>' and leading/trailing whitespace
                else:
                    sequence_lines.append(line)
        if not header or not sequence_lines:
            sys.exit(f"[ERROR] No valid FASTA sequence found in {fasta_path}")
    except FileNotFoundError:
        sys.exit(f"[ERROR] FASTA file not found: {fasta_path}")
    except Exception as e:
        sys.exit(f"[ERROR] An error occurred while loading the FASTA file: {e}")
    return header, ''.join(sequence_lines)


def mask_sequence(sequence, regions):
    """
    Masks regions of a sequence with 'N's.

    Args:
        sequence (str): The original sequence string.
        regions (list): A list of (start, end) tuples (0-based) to mask.

    Returns:
        str: The masked sequence string.
    """
    masked_list = list(sequence)
    for start, end in regions:
        # Ensure coordinates are within bounds
        start = max(0, start)
        end = min(len(sequence), end)
        if start < end: # Only mask if region is valid
            masked_list[start:end] = ['N'] * (end - start)
    return "".join(masked_list)


def write_fasta(header, sequence, output_path, wrap_len=60):
    """
    Writes a FASTA sequence to a file.

    Args:
        header (str): The FASTA header (without '>').
        sequence (str): The sequence string.
        output_path (str): Path to the output FASTA file.
        wrap_len (int): Length at which to wrap sequence lines.
    """
    try:
        with open(output_path, 'w') as out_f:
            out_f.write(f">{header}\n")
            for i in range(0, len(sequence), wrap_len):
                out_f.write(sequence[i:i+wrap_len] + "\n")
    except Exception as e:
        print(f"[ERROR] Failed to write FASTA file {output_path}: {e}", file=sys.stderr)


# --- Main execution ---
print(f"[INFO] Loading FASTA file: {fasta_file}")
header, sequence = load_fasta(fasta_file)
print(f"[INFO] FASTA sequence loaded. Length: {len(sequence)}")

# Process and mask for 100% identity regions
print("[INFO] Processing regions with 100% identity...")
regions_100 = parse_tsv(tsv_file, min_id=Decimal("100.0"), max_id=Decimal("100.0"))
masked_100 = mask_sequence(sequence, regions_100)
output_path_100 = os.path.join(output_dir, f"{prefix}_masked_100.fasta")
write_fasta(header, masked_100, output_path_100)
print(f"[INFO] Masked 100% identity sequences saved to: {output_path_100}")

# Process and mask for 95.0% to 99.9% identity regions
print("[INFO] Processing regions with 95.0% to 99.9% identity...")
regions_95_999 = parse_tsv(tsv_file, min_id=Decimal("95.0"), max_id=Decimal("99.9")) # Use 99.9 for strict upper bound
masked_95_999 = mask_sequence(sequence, regions_95_999)
output_path_95_999 = os.path.join(output_dir, f"{prefix}_masked_95_99.9.fasta") # Filename reflects change
write_fasta(header, masked_95_999, output_path_95_999)
print(f"[INFO] Masked 95.0% to 99.9% identity sequences saved to: {output_path_95_999}")

print("\n-----------------------------")
print("Genome masking completed.")
