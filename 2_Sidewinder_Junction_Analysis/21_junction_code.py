import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

"""
# This script analyzes junction connections from sequencing reads of DNA
# assembly products. It:
#   1) Builds all possible fragment-junctions (+/- x bp from ends),
#   2) BLASTs junctions vs junctions to pick a bitscore threshold,
#   3) BLASTs junctions vs reads,
#   4) Generates a junction heatmap and exports misassemblies to FASTA.
"""

# -------------------------
# 1) INPUTS / PARAMETERS
# -------------------------

# TODO: Change these paths to your inputs
fastq_reads = "../data/mScar.fastq"
csv_fragments = "fragments.csv"

# Junction window size (per end). Junction length = 2 * x
x = 25

# Outputs
junction_path_name = "junctions.fasta"
heatmap_name = "mScar_heatmap.png"
excel_name = "mScar_heatmap.xlsx"

# -------------------------
# 2) LOAD FRAGMENTS
# -------------------------

df = pd.read_csv(csv_fragments, header=None)
fragments = [str(s).upper() for s in df[0].tolist()]
num_fragments = len(fragments)

# -------------------------
# 3) GENERATE JUNCTIONS
# -------------------------

def generate_x_lists(fragments, x=25):
    """
    Generate three lists from the given n number of fragments (f1..fn):
      1) left_ends: last x bases of each fragment.
      2) right_ends: first x bases of each fragment.
      3) all_junctions: for every i,j, left_ends[i] + right_ends[j].
    """
    for i, seq in enumerate(fragments, start=1):
        if len(seq) < x:
            raise ValueError(f"Fragment f{i} is shorter than x={x} bases.")

    left_ends = [seq[-x:] for seq in fragments]
    right_ends = [seq[:x] for seq in fragments]

    all_junctions = []
    for i in range(len(fragments)):
        for j in range(len(fragments)):
            all_junctions.append(left_ends[i] + right_ends[j])

    return left_ends, right_ends, all_junctions


left_ends, right_ends, all_junctions = generate_x_lists(fragments, x)

print(f"Left ends (last {x} bases of f1..{len(fragments)}):")
for i, seq in enumerate(left_ends, 1):
    print(f"F{i}_last{x} = {seq}")

print(f"Right ends (first {x} bases of f1..{len(fragments)}):")
for i, seq in enumerate(right_ends, 1):
    print(f"F{i}_first{x} = {seq}")

print(f"\nAll possible junctions ({len(fragments) * len(fragments)} total):")
print(all_junctions[:6], "...")

def write_junctions_to_fasta(all_junctions, outfile="junctions.fasta"):
    """
    Write junctions to FASTA with headers like >junction_1-1 ... >junction_n-n.
    """
    n = int(np.sqrt(len(all_junctions)))
    with open(outfile, "w") as f:
        idx = 0
        for i in range(n):
            for j in range(n):
                f.write(f">junction_{i+1}-{j+1}\n{all_junctions[idx]}\n")
                idx += 1

write_junctions_to_fasta(all_junctions, junction_path_name)

# -------------------------
# 4) BLAST HELPERS
# -------------------------

def run_cmd(cmd_list, cwd=None):
    subprocess.run(cmd_list, check=True, cwd=cwd)

def fastq_to_fasta(fastq_file):
    """Convert FASTQ -> FASTA (simple 4-line FASTQ blocks)."""
    base_name = os.path.splitext(fastq_file)[0]
    fasta_file = base_name + ".fasta"
    with open(fastq_file, 'r') as fq, open(fasta_file, 'w') as fa:
        while True:
            header = fq.readline().strip()
            if not header:
                break
            sequence = fq.readline().strip()
            fq.readline()  # '+'
            fq.readline()  # quality
            fa.write(">" + header[1:].split()[0] + "\n")
            fa.write(sequence + "\n")
    print(f"Converted {fastq_file} -> {fasta_file}")
    return fasta_file

def parse_blast_outfmt6(filename):
    """Parse BLAST tabular (outfmt 6) file into a pandas DataFrame."""
    colnames = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    return pd.read_csv(filename, sep="\t", names=colnames, comment='#', engine='python')

# -------------------------
# 5) JUNCTION vs JUNCTION
# -------------------------

# Self-BLAST to estimate bitscore threshold
run_cmd([
    "blastn",
    "-task", "blastn-short",
    "-query", junction_path_name,
    "-subject", junction_path_name,
    "-out", "junction_against_junctions.csv",
    "-outfmt", "6",
    "-word_size", "7",
    "-reward", "1",
    "-penalty", "-3",
    "-gapopen", "5",
    "-gapextend", "2",
    "-dust", "no",
    "-evalue", "1e-3",
    "-max_target_seqs", "40000000"
])

junction_against_junction = parse_blast_outfmt6('junction_against_junctions.csv')
junction_against_junction['pairs'] = (
    junction_against_junction['qseqid'] + '+' + junction_against_junction['sseqid']
)
junction_against_junction_filtered = (
    junction_against_junction.loc[junction_against_junction.groupby('pairs')['bitscore'].idxmax()]
    .copy()
)
junction_against_junction_filtered[['i', 'j']] = (
    junction_against_junction_filtered['qseqid']
    .str.extract(r'junction_(\d+)-(\d+)').astype(int)
)

plt.figure(figsize=(8, 5))
plt.hist(junction_against_junction_filtered['bitscore'],
         bins=int(max(junction_against_junction_filtered.bitscore) / 2))
plt.title('Junction vs Junction Bitscore Distribution')
plt.xlabel('Bitscore')
plt.ylabel('Count')
plt.tight_layout()
#plt.show()

#Select an appropriate cutoff based on junction vs. junction blast result

# -------------------------
# 6) JUNCTIONS vs READS
# -------------------------


fasta_reads = fastq_to_fasta(fastq_reads)

# BLAST junctions vs reads (use -subject on the FASTA; no DB needed)
run_cmd([
    "blastn",
    "-task", "blastn-short",
    "-query", junction_path_name,
    "-subject", fasta_reads,
    "-out", "junction_against_reads.csv",
    "-outfmt", "6",
    "-word_size", "7",
    "-reward", "1",
    "-penalty", "-3",
    "-gapopen", "5",
    "-gapextend", "2",
    "-dust", "no",
    "-evalue", "1e-3",
    "-max_target_seqs", "40000000"
])

junction_against_reads = parse_blast_outfmt6("junction_against_reads.csv")
junction_against_reads['pairs'] = (
    junction_against_reads['qseqid'] + '+' + junction_against_reads['sseqid']
)
junction_against_reads_filtered = (
    junction_against_reads.loc[junction_against_reads.groupby('pairs')['bitscore'].idxmax()]
    .copy()
)
junction_against_reads_filtered[['i', 'j']] = (
    junction_against_reads_filtered['qseqid']
    .str.extract(r'junction_(\d+)-(\d+)').astype(int)
)

plt.figure(figsize=(10, 6))
plt.hist(junction_against_reads_filtered['bitscore'], bins=50, alpha=0.7,
         edgecolor='black')
plt.title('Junction vs Reads Bitscore Distribution')
plt.xlabel('Bitscore')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)
plt.tight_layout()
#plt.show()

# -------------------------
# 7) HEATMAP
# -------------------------

def generate_heatmap(heatmap_name, num_fragments, df_hits,
                     bitscore_cutoff=52, length_cutoff=38,
                     trace=None, excel_name=None):
    """Build a junction heatmap and export counts to Excel."""
    df = df_hits.copy()
    if trace:
        df = df[df['sseqid'] == trace]
    df['length'] = df['length'].astype(int)
    df['bitscore'] = df['bitscore'].astype(float)
    df = df[(df['bitscore'] > bitscore_cutoff) & (df['length'] > length_cutoff)]

    filtered_counts = df['qseqid'].value_counts()
    junction_counts = {}
    for junction, count in filtered_counts.items():
        i, j = map(int, junction.replace('junction_', '').split('-'))
        junction_counts[(i, j)] = int(count)

    heatmap_matrix = np.zeros((num_fragments, num_fragments), dtype=int)
    for (i, j), count in junction_counts.items():
        if 1 <= i <= num_fragments and 1 <= j <= num_fragments:
            heatmap_matrix[i - 1, j - 1] = count

    plt.figure(figsize=(5, 5) if trace else (20, 20))
    sns.heatmap(
        heatmap_matrix, annot=True, fmt="g", cmap="viridis",
        xticklabels=range(1, num_fragments + 1),
        yticklabels=range(1, num_fragments + 1),
        annot_kws={"size": 25 if not trace else 12}
    )
    plt.title(trace if trace else "Junction Heatmap")
    plt.xlabel("Fragment j")
    plt.ylabel("Fragment i")
    plt.savefig(heatmap_name, dpi=300)

    # Excel export
    if excel_name is None:
        excel_name = heatmap_name[:-4] + ".xlsx" if heatmap_name.lower().endswith(".png") else f"{heatmap_name}.xlsx"
    labels = list(range(1, num_fragments + 1))
    matrix_df = pd.DataFrame(heatmap_matrix, index=labels, columns=labels)
    matrix_df.index.name = "Fragment i"
    matrix_df.columns.name = "Fragment j"

    long_df = (
        matrix_df.reset_index()
                 .melt(id_vars="Fragment i", var_name="Fragment j", value_name="count")
                 .sort_values(["Fragment i", "Fragment j"], ignore_index=True)
    )

    with pd.ExcelWriter(excel_name, engine="xlsxwriter") as writer:
        matrix_df.to_excel(writer, sheet_name="matrix")
        long_df.to_excel(writer, sheet_name="long_format", index=False)
        workbook = writer.book
        int_fmt = workbook.add_format({"num_format": "0"})
        ws = writer.sheets["matrix"]
        _, ncols = matrix_df.shape
        ws.set_column(1, ncols, 10, int_fmt)

    return heatmap_matrix

bitscore_cutoff = 52
length_cutoff = 38

heatmap = generate_heatmap(heatmap_name, num_fragments, junction_against_reads_filtered,
                     bitscore_cutoff=bitscore_cutoff, length_cutoff=length_cutoff,
                     trace=None, excel_name=excel_name)

# -------------------------
# 8) EXPORT MISASSEMBLIES
# -------------------------

# Misassemblies where j != i + 1 (not adjacent fragments in expected order)
misassembly_filtered = junction_against_reads_filtered[
    junction_against_reads_filtered['j'] != (junction_against_reads_filtered['i'] + 1)
]
blast_misassemblies = list(
    misassembly_filtered[(misassembly_filtered['bitscore'] > bitscore_cutoff) & (misassembly_filtered['length'] > length_cutoff)].sseqid
)
output_file = "misassemblies.fasta"
found_sequences = {}
count = 0

# Directly parse the reads FASTA we created earlier
for record in SeqIO.parse(fasta_reads, "fasta"):
    if record.id in blast_misassemblies:
        count += 1
        print(f"Matched sequence #{count}: {record.id}")
        found_sequences[record.id] = str(record.seq)

with open(output_file, 'w') as f:
    for header, sequence in found_sequences.items():
        f.write(f">{header}\n{sequence}\n")

print(f"Extracted {len(found_sequences)} sequences. Output written to {output_file}.")
