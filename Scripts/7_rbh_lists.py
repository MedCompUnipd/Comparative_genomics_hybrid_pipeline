# -*- coding: utf-8 -*-
import os
import re
import pandas as pd

# ==== Input paths from user ====
contig_vs_genes_dir = input("Enter path to 'contig_vs_genes' directory: ").strip()
genes_vs_contig_dir = input("Enter path to 'genes_vs_contig' directory: ").strip()
output_dir = input("Enter path to output directory: ").strip()
os.makedirs(output_dir, exist_ok=True)

# ==== Function to parse glsearch output ====
def parse_glsearch(filepath):
    identity_pattern = re.compile(r'([\d\.]+)% identity')
    overlap_pattern = re.compile(r'overlap \((\d+)-(\d+):(\d+)-(\d+)\)')
    hit_name_pattern = re.compile(r'^>>(\S+)')

    identity = None
    query_coords = None
    subject_coords = None
    hit_name = None

    try:
        with open(filepath, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if hit_name is None and line.startswith(">>"):
                    hit_match = hit_name_pattern.match(line)
                    if hit_match:
                        hit_name = hit_match.group(1)

                if identity is None:
                    id_match = identity_pattern.search(line)
                    if id_match:
                        identity = float(id_match.group(1))

                if query_coords is None and subject_coords is None:
                    overlap_match = overlap_pattern.search(line)
                    if overlap_match:
                        q_start = int(overlap_match.group(1))
                        q_end = int(overlap_match.group(2))
                        query_coords = (min(q_start, q_end), max(q_start, q_end))

                        s_start = int(overlap_match.group(3))
                        s_end = int(overlap_match.group(4))
                        subject_coords = (min(s_start, s_end), max(s_start, s_end))

                if identity is not None and query_coords is not None and subject_coords is not None and hit_name is not None:
                    break
    except Exception as e:
        print(f"Error reading {filepath}: {e}")

    if identity is None or query_coords is None or subject_coords is None or hit_name is None:
        print(f"Warning: Missing data in {filepath} (identity={identity}, query_coords={query_coords}, subject_coords={subject_coords}, hit={hit_name})")

    return identity, query_coords, subject_coords, hit_name

# ==== Helper functions ====
def overlaps(coords1, coords2, tol=5):
    return not (coords1[1] + tol < coords2[0] or coords2[1] + tol < coords1[0])

def is_rbh(c_identity, g_identity, contig_actual_coords, gene_actual_coords, identity_tol=1.0, coord_tol=5):
    if None in (c_identity, g_identity, contig_actual_coords, gene_actual_coords):
        return False
    identity_close = abs(c_identity - g_identity) <= identity_tol
    coords_overlap = overlaps(contig_actual_coords, gene_actual_coords, tol=coord_tol)
    return identity_close and coords_overlap

# ==== Process files ====
rbh_records = []
total_files = 0
skipped_files = 0

for filename in os.listdir(contig_vs_genes_dir):
    if not filename.endswith(".glsearch.out"):
        continue

    total_files += 1
    contig_file = os.path.join(contig_vs_genes_dir, filename)
    gene_file = os.path.join(genes_vs_contig_dir, filename)

    if not os.path.exists(gene_file):
        print(f"Warning: Matching gene file not found for {filename}, skipping.")
        skipped_files += 1
        continue

    basename = filename.replace(".glsearch.out", "")
    identity_contig_query, contig_q_coords, gene_s_coords_from_contig_file, _ = parse_glsearch(contig_file)
    identity_gene_query, gene_q_coords, contig_s_coords_from_gene_file, _ = parse_glsearch(gene_file)

    if is_rbh(identity_contig_query, identity_gene_query, contig_q_coords, gene_q_coords):
        # Compare coordinates between CvsG and GvsC: must be identical
        if (contig_q_coords == gene_q_coords) and (gene_s_coords_from_contig_file == contig_s_coords_from_gene_file):
            gene_start = contig_q_coords[0]
            gene_end = contig_q_coords[1]
            contig_start = gene_s_coords_from_contig_file[0]
            contig_end = gene_s_coords_from_contig_file[1]

            gene_len = abs(gene_end - gene_start) + 1
            contig_len = abs(contig_end - contig_start) + 1
            query_coverage = (contig_len / gene_len) * 100 if gene_len > 0 else 0.0

            avg_identity = (identity_contig_query + identity_gene_query) / 2

            rbh_records.append({
                "Gene": basename,
                "Identity": avg_identity,
                "Gene_Start": gene_start,
                "Gene_End": gene_end,
                "Contig_Start": contig_start,
                "Contig_End": contig_end,
                "Query_Coverage": round(query_coverage, 2)
            })
        else:
            print(f"Skipped {basename}: coordinates differ between CvsG and GvsC")
    else:
        print(f"Skipped non-RBH: {basename} "
              f"(Identities: {identity_contig_query}, {identity_gene_query}; "
              f"Contig Query Coords: {contig_q_coords}, Gene Query Coords: {gene_q_coords})")

# ==== Create and filter DataFrame ====
df = pd.DataFrame(rbh_records)
desired_order = ["Gene", "Identity", "Gene_Start", "Gene_End", "Contig_Start", "Contig_End", "Query_Coverage"]
df = df[desired_order]

df_100 = df[df["Identity"] >= 99.99]
df_95_99 = df[(df["Identity"] >= 95.0) & (df["Identity"] < 99.99)]
df_below_95 = df[df["Identity"] < 95.0]

# ==== Write Excel ====
excel_path = os.path.join(output_dir, "RBH_results.xlsx")
with pd.ExcelWriter(excel_path) as writer:
    df.to_excel(writer, sheet_name="All_RBH", index=False)
    df_100.to_excel(writer, sheet_name="Identity_100", index=False)
    df_95_99.to_excel(writer, sheet_name="Identity_95_99", index=False)
    df_below_95.to_excel(writer, sheet_name="Identity_below_95", index=False)

# ==== Write TSV ====
df.to_csv(os.path.join(output_dir, "RBH_all.tsv"), sep="\t", index=False)
df_100.to_csv(os.path.join(output_dir, "RBH_100.tsv"), sep="\t", index=False)
df_95_99.to_csv(os.path.join(output_dir, "RBH_95_99.tsv"), sep="\t", index=False)
df_below_95.to_csv(os.path.join(output_dir, "RBH_below_95.tsv"), sep="\t", index=False)

# ==== Summary ====
print(f"Processed {total_files} files.")
print(f"Skipped {skipped_files} files due to missing pairs.")
print(f"Total RBHs found: {len(df)}")
print(f" - 100% identity (≥99.99): {len(df_100)}")
print(f" - 95-99% identity: {len(df_95_99)}")
print(f" - Below 95% identity: {len(df_below_95)}")
print(f"Excel file saved to: {excel_path}")
print("TSV files saved in directory:", output_dir)
