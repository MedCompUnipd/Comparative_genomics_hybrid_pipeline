#!/usr/bin/env python3

import os
import re
import pandas as pd
import sys # Import sys for exit

# ==== Conda Activation (optional, usually handled by calling shell script) ====
# This part is generally handled by the calling shell script using `conda activate`.
# If you run this script directly, ensure your environment is active.
# No direct conda activation commands inside the Python script for flexibility.

# ==== Input paths from user ====
print("--- Reciprocal Best Hit (RBH) List Generation Script ---")
print("-------------------------------------------------------")

contig_vs_genes_dir = input("Enter path to 'contig_vs_genes' directory: ").strip()
genes_vs_contig_dir = input("Enter path to 'genes_vs_contig' directory: ").strip()
output_dir = input("Enter path to output directory for RBH results: ").strip()

# Validate input directories
if not os.path.isdir(contig_vs_genes_dir):
    print(f"[ERROR] 'contig_vs_genes' directory not found: {contig_vs_genes_dir}", file=sys.stderr)
    sys.exit(1)
if not os.path.isdir(genes_vs_contig_dir):
    print(f"[ERROR] 'genes_vs_contig' directory not found: {genes_vs_contig_dir}", file=sys.stderr)
    sys.exit(1)

os.makedirs(output_dir, exist_ok=True)
print(f"[INFO] Output directory: {output_dir}")

# ==== Function to parse glsearch output ====
def parse_glsearch(filepath):
    """
    Parses a glsearch output file to extract hit name, identity, and coordinates.
    Assumes -m 0 output format.
    """
    identity_pattern = re.compile(r'([\d\.]+)% identity')
    coords_pattern = re.compile(r'overlap \((\d+)-(\d+):(\d+)-(\d+)\)')
    hit_name_pattern = re.compile(r'^>>(\S+)') # Captures the first non-whitespace word after ">>"

    identity = None
    query_coords = None # Coordinates for the query sequence
    hit_coords = None   # Coordinates for the hit sequence
    hit_name = None

    try:
        with open(filepath, 'r', encoding="utf-8", errors="ignore") as f:
            for line in f:
                if hit_name is None and line.startswith(">>"):
                    hit_match = hit_name_pattern.match(line)
                    if hit_match:
                        hit_name = hit_match.group(1)
                
                if identity is None:
                    id_match = identity_pattern.search(line)
                    if id_match:
                        identity = float(id_match.group(1))
                
                if coords_pattern.search(line): # There can be multiple overlap lines, take the first one found
                    coord_match = coords_pattern.search(line)
                    if coord_match:
                        query_start = int(coord_match.group(1))
                        query_end = int(coord_match.group(2))
                        hit_start = int(coord_match.group(3))
                        hit_end = int(coord_match.group(4))
                        query_coords = f"{query_start}-{query_end}"
                        hit_coords = f"{hit_start}-{hit_end}"
                        # Once we found the coordinates, we can stop searching for them
                        # This assumes the first overlap line is the primary one
                        if identity is not None and hit_name is not None:
                            break # Found all required info
    except Exception as e:
        print(f"[ERROR] Failed to parse {filepath}: {e}", file=sys.stderr)
        return None, None, None, None, None # Return None for all values on error

    return hit_name, identity, query_coords, hit_coords

# ==== Process glsearch output files ====
print("[INFO] Processing glsearch output files...")
contig_vs_genes_results = {} # key: contig, value: {gene: (identity, contig_coords, gene_coords)}
genes_vs_contig_results = {} # key: gene, value: {contig: (identity, gene_coords, contig_coords)}

# Process contig_vs_genes directory
for filename in os.listdir(contig_vs_genes_dir):
    if filename.endswith(".glsearch.out"):
        filepath = os.path.join(contig_vs_genes_dir, filename)
        # In contig_vs_genes, query is the contig, hit is the gene
        # The filename usually indicates the gene name (e.g., GENE_A.glsearch.out)
        # The hit_name in glsearch output will be the contig name (usually 'denovo_contig' from 6_rbh.sh)
        
        gene_name_from_file = filename.replace(".glsearch.out", "")
        # The parse_glsearch function returns (hit_name, identity, query_coords, hit_coords)
        # Here, query is contig, hit is gene. But the glsearch output (m 0) has the query header as the first line with '>>'
        # The parsing function gets hit_name from the '>>' line.
        # Let's adjust parsing function's return to be more generic, and then map them here.
        # For glsearch36 -m 0, the first '>>' is the hit. The query is implicit from input file.
        # Re-evaluating parse_glsearch: it parses the *hit* name and its coords.
        # For `glsearch36 query_file hit_file`, hit_name is from hit_file.
        # So in contig_vs_genes: query=contig, hit=gene.
        # And in genes_vs_contig: query=gene, hit=contig.

        # Let's refine parse_glsearch to be more explicit about query vs hit
        # The -m 0 format has query and hit coordinates
        
        # New approach: Parse the actual query and hit names from the file.
        # For -m 0 format:
        # >><hit_name>
        # ...
        # optimal alignment of: <query_name> to <hit_name>
        # ...
        # query: <query_start>-<query_end>  hit: <hit_start>-<hit_end>

        # Redefine parse_glsearch for better clarity
        def parse_glsearch_v2(filepath):
            query_name = None
            hit_name = None
            identity = None
            query_coords = None
            hit_coords = None

            query_pattern = re.compile(r'optimal alignment of:\s*(\S+)\s+to')
            hit_pattern = re.compile(r'^>>(\S+)')
            identity_pattern = re.compile(r'([\d\.]+)% identity')
            coords_line_pattern = re.compile(r'query:\s*(\d+)-(\d+)\s+hit:\s*(\d+)-(\d+)')

            try:
                with open(filepath, 'r', encoding="utf-8", errors="ignore") as f:
                    content = f.read() # Read entire content for multiple matches

                # Extract hit name from '>>' line (first one)
                hit_match = hit_pattern.search(content)
                if hit_match:
                    hit_name = hit_match.group(1)
                
                # Extract query name from "optimal alignment of" line
                query_match = query_pattern.search(content)
                if query_match:
                    query_name = query_match.group(1)

                # Extract identity
                id_match = identity_pattern.search(content)
                if id_match:
                    identity = float(id_match.group(1))

                # Extract coordinates
                coords_match = coords_line_pattern.search(content)
                if coords_match:
                    query_coords = f"{coords_match.group(1)}-{coords_match.group(2)}"
                    hit_coords = f"{coords_match.group(3)}-{coords_match.group(4)}"
            except Exception as e:
                print(f"[ERROR] Failed to parse {filepath}: {e}", file=sys.stderr)
                return None, None, None, None, None

            return query_name, hit_name, identity, query_coords, hit_coords


# Process contig_vs_genes directory (query is contig, hit is gene)
for filename in os.listdir(contig_vs_genes_dir):
    if filename.endswith(".glsearch.out"):
        filepath = os.path.join(contig_vs_genes_dir, filename)
        query_name, hit_name, identity, query_coords, hit_coords = parse_glsearch_v2(filepath)
        if all([query_name, hit_name, identity is not None, query_coords, hit_coords]):
            # The query is the contig, the hit is the gene.
            # In contig_vs_genes_results, key by contig name, then gene name.
            # E.g., contig_vs_genes_results['denovo_contig']['geneA'] = (identity, contig_coords, gene_coords)
            if query_name not in contig_vs_genes_results:
                contig_vs_genes_results[query_name] = {}
            contig_vs_genes_results[query_name][hit_name] = (identity, query_coords, hit_coords)

# Process genes_vs_contig directory (query is gene, hit is contig)
for filename in os.listdir(genes_vs_contig_dir):
    if filename.endswith(".glsearch.out"):
        filepath = os.path.join(genes_vs_contig_dir, filename)
        query_name, hit_name, identity, query_coords, hit_coords = parse_glsearch_v2(filepath)
        if all([query_name, hit_name, identity is not None, query_coords, hit_coords]):
            # The query is the gene, the hit is the contig.
            # In genes_vs_contig_results, key by gene name, then contig name.
            # E.g., genes_vs_contig_results['geneA']['denovo_contig'] = (identity, gene_coords, contig_coords)
            if query_name not in genes_vs_contig_results:
                genes_vs_contig_results[query_name] = {}
            genes_vs_contig_results[query_name][hit_name] = (identity, query_coords, hit_coords)

# ==== Identify Reciprocal Best Hits (RBH) ====
print("[INFO] Identifying Reciprocal Best Hits (RBH)...")
rbh_records = []

# Find best hit for each gene in contig_vs_genes_results (Contig as query -> Gene as hit)
# This means: for a given contig, what gene does it align best to?
contig_best_hits = {} # {contig_name: (best_gene_name, max_identity)}
for contig, hits_to_genes in contig_vs_genes_results.items():
    if not hits_to_genes: continue
    best_gene = None
    max_identity = -1
    for gene, (identity, _, _) in hits_to_genes.items():
        if identity > max_identity:
            max_identity = identity
            best_gene = gene
    if best_gene:
        contig_best_hits[contig] = (best_gene, max_identity)

# Find best hit for each gene in genes_vs_contig_results (Gene as query -> Contig as hit)
# This means: for a given gene, what contig does it align best to?
gene_best_hits = {} # {gene_name: (best_contig_name, max_identity)}
for gene, hits_to_contigs in genes_vs_contig_results.items():
    if not hits_to_contigs: continue
    best_contig = None
    max_identity = -1
    for contig, (identity, _, _) in hits_to_contigs.items():
        if identity > max_identity:
            max_identity = identity
            best_contig = contig
    if best_contig:
        gene_best_hits[gene] = (best_contig, max_identity)

# Compare best hits to find RBHs
for contig, (best_hit_gene, contig_to_gene_identity) in contig_best_hits.items():
    if best_hit_gene in gene_best_hits:
        best_hit_contig_for_gene, gene_to_contig_identity = gene_best_hits[best_hit_gene]
        if contig == best_hit_contig_for_gene:
            # It's an RBH!
            # Get coordinates and identity for the RBH
            
            # Identity and coordinates from Contig -> Gene alignment
            identity_c_g, contig_coords_c_g, gene_coords_c_g = contig_vs_genes_results[contig][best_hit_gene]

            # Identity and coordinates from Gene -> Contig alignment
            identity_g_c, gene_coords_g_c, contig_coords_g_c = genes_vs_contig_results[best_hit_gene][contig]
            
            # For RBH, identities should be very similar. Take one, e.g., from contig->gene
            avg_identity = round((identity_c_g + identity_g_c) / 2, 2)
            
            rbh_records.append({
                "Contig": contig,
                "Gene": best_hit_gene,
                "Identity": avg_identity,
                "Contig_start": contig_coords_c_g.split('-')[0],
                "Contig_end": contig_coords_c_g.split('-')[1],
                "Gene_start": gene_coords_c_g.split('-')[0],
                "Gene_end": gene_coords_c_g.split('-')[1]
            })
            print(f"[RBH FOUND] Contig: {contig}, Gene: {best_hit_gene}, Identity: {avg_identity}%")


# Sort RBH records by identity (descending) then by gene name
rbh_records.sort(key=lambda x: (x["Identity"], x["Gene"]), reverse=True)


# ==== Create DataFrame and filter ===
df = pd.DataFrame(rbh_records)

if df.empty:
    print("[WARNING] No Reciprocal Best Hits found. Output files will be empty or not created.")
    sys.exit(0) # Exit gracefully if no RBHs

df_100 = df[df["Identity"] == 100.0]
df_95_99 = df[(df["Identity"] >= 95.0) & (df["Identity"] < 99.99)] # Corrected to < 99.99 for strict range
df_below_95 = df[df["Identity"] < 95.0]

# ==== Write Excel file ====
excel_path = os.path.join(output_dir, "RBH_results.xlsx")
try:
    with pd.ExcelWriter(excel_path) as writer:
        df.to_excel(writer, sheet_name="All_RBH", index=False)
        df_100.to_excel(writer, sheet_name="Identity_100", index=False)
        df_95_99.to_excel(writer, sheet_name="Identity_95_99", index=False)
        df_below_95.to_excel(writer, sheet_name="Identity_below_95", index=False)
    print(f"[INFO] Excel output written to: {excel_path}")
except Exception as e:
    print(f"[ERROR] Failed to write Excel file: {e}", file=sys.stderr)

# ==== Write TSV files ====
try:
    df.to_csv(os.path.join(output_dir, "RBH_all.tsv"), sep="\t", index=False)
    df_100.to_csv(os.path.join(output_dir, "RBH_100.tsv"), sep="\t", index=False)
    df_95_99.to_csv(os.path.join(output_dir, "RBH_95_99.tsv"), sep="\t", index=False)
    df_below_95.to_csv(os.path.join(output_dir, "RBH_below_95.tsv"), sep="\t", index=False)
    print(f"[INFO] TSV outputs written to: {output_dir}")
except Exception as e:
    print(f"[ERROR] Failed to write TSV files: {e}", file=sys.stderr)

print("\n-------------------------------------------------------")
print("RBH list generation completed.")
