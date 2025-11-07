import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

"""
This script conducts fragment analysis of Sidewinder assembly reads .
It:
1) Runs BLAST alignment of every fragment to every read.
2) Parses BLAST hits and classifies reads as: confirmed_full_correct, confirmed_partial_correct, to_be_checked, or unusable,
3) Saves classifications and summary CSVs,
4) Exports reads needing review to FASTA.
"""


# ----------------------------
# Inputs / Parameters
# ----------------------------
csv_fragments = "fragments.csv"
df = pd.read_csv(csv_fragments, header=None)
fragments = [str(s).upper() for s in df[0].tolist()]


reads_fastq = "../data/ApoE_GelExtracted.fastq"
fragments_fasta = "fragments.fasta"
blast_out = "fragments_against_reads.csv"

# Output FASTA files (by class)
out_dir = "./"
os.makedirs(out_dir, exist_ok=True)
out_remaining = os.path.join(out_dir, "Sequences_to_be_checked.fasta")

# ----------------------------
# Helpers
# ----------------------------

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

reads_fasta = fastq_to_fasta(reads_fastq)

def write_fragments_to_fasta(fragments, outfile="fragments.fasta"):
    """Write fragments to FASTA with headers >fragment_1..N."""
    with open(outfile, "w") as f:
        for i, seq in enumerate(fragments, start=1):
            f.write(f">fragment_{i}\n{seq}\n")

def run_cmd(cmd_list):
    """Run a shell command, raising an error if it fails."""
    subprocess.run(cmd_list, check=True)

def parse_blast_outfmt6(filename):
    """Parse BLAST outfmt 6 to a pandas DataFrame."""
    colnames = [
        "qseqid", "sseqid", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore"
    ]
    return pd.read_csv(filename, sep="\t", names=colnames, comment="#", engine="python")

def parse_blast_results(df):
    """Normalize BLAST hits and apply minimal length filter."""
    df = df.copy()
    df["fragment_num"] = df["qseqid"].str.extract(r"fragment_(\d+)").astype(int)
    df["smin"] = df[["sstart", "send"]].min(axis=1)
    df["smax"] = df[["sstart", "send"]].max(axis=1)
    df = df[df["length"] > 40]  # keep reasonably long matches
    return df

def classify_read_hits(group):
    """Classify a read based on the ordered list of fragment hits along the read."""
    ordered_hits = group.sort_values("smin")["fragment_num"].tolist()
    if not ordered_hits:
        return "unusable"

    def has_full_subsequence(lst, target):
        for i in range(len(lst) - len(target) + 1):
            if lst[i:i + len(target)] == target:
                return True
        return False

    # exact or embedded exact sequence 1..12 or 12..1
    ascending = list(range(1, 13))
    descending = ascending[::-1]
    if ordered_hits == ascending or ordered_hits == descending:
        return "confirmed_full_correct"
    if has_full_subsequence(ordered_hits, ascending) or has_full_subsequence(ordered_hits, descending):
        return "confirmed_full_correct"

    def check_for_jump(lst):
        return any(1 < abs(b - a) for a, b in zip(lst, lst[1:]))

    def has_direction_flip(lst):
        if len(lst) < 2:
            return False
        direction = None  # 'up' or 'down'
        for a, b in zip(lst, lst[1:]):
            if b == a:
                return True  # duplicate fragment number indicates an issue
            if b > a:
                if direction == "down":
                    return True
                direction = "up"
            else:  # b < a
                if direction == "up":
                    return True
                direction = "down"
        return False

    if check_for_jump(ordered_hits):
        return "to_be_checked"
    if has_direction_flip(ordered_hits):
        return "to_be_checked"

    return "confirmed_partial_correct"

def analyze_blast(blast_df, reads_fasta):
    """Return per-read classification and a summary table."""
    df = parse_blast_results(blast_df)
    grouped = df.groupby("sseqid")
    classification = grouped.apply(classify_read_hits)
    classification.name = "assembly_class"

    # any reads with no hits → unusable
    all_read_ids = {rec.id for rec in SeqIO.parse(reads_fasta, "fasta")}
    classified_ids = set(classification.index)
    for read_id in (all_read_ids - classified_ids):
        classification.loc[read_id] = "unusable"

    summary = classification.value_counts().rename_axis("assembly_class").reset_index(name="count")
    return classification, summary

def write_subset_fasta(src_fasta, dest_fasta, wanted_ids):
    """Write subset of records (by id) to a FASTA."""
    wanted = set(wanted_ids)
    with open(dest_fasta, "w") as out_f:
        for record in SeqIO.parse(src_fasta, "fasta"):
            if record.id in wanted:
                SeqIO.write(record, out_f, "fasta")

def sort_fasta_by_header(input_file, output_file):
    """Sort a FASTA file numerically by header id (assumes ids are integers)."""
    records = list(SeqIO.parse(input_file, "fasta"))
    try:
        records.sort(key=lambda r: int(r.id))
    except ValueError:
        # Fallback to lexicographic sort if IDs are not integers
        records.sort(key=lambda r: r.id)
    SeqIO.write(records, output_file, "fasta")

# ----------------------------
# Main flow
# ----------------------------
if __name__ == "__main__":
    # 1) Write fragments to FASTA
    write_fragments_to_fasta(fragments, fragments_fasta)

    # 2) Run BLAST: fragments vs reads
    run_cmd([
        "blastn",
        "-query", fragments_fasta,
        "-subject", reads_fasta,
        "-word_size", "7",
        "-out", blast_out,
        "-outfmt", "6",
        "-max_target_seqs", "4000000",
        "-dust", "no",
    ])

    # 3) Parse BLAST and classify reads
    blast_df = parse_blast_outfmt6(blast_out)
    classification, summary = analyze_blast(blast_df, reads_fasta)

    print("Summary of read classifications:")
    print(summary.to_string(index=False))

    # Save classifications
    classification.to_csv("read_classification.csv")
    summary.to_csv("read_classification_summary.csv", index=False)

    # 4) Export remaining FASTA subsets to be categorized manually afterwards
    headers_remaining= classification[classification.eq("to_be_checked")].index.tolist()

    write_subset_fasta(reads_fasta, out_remaining, headers_remaining)

    # 6) Quick bar plot from computed summary
    plt.figure(figsize=(8, 6))
    plt.bar(summary["assembly_class"], summary["count"])
    for i, v in enumerate(summary["count"]):
        plt.text(i, v, str(v), ha="center", va="bottom")
    plt.ylabel("Read Count")
    plt.title("Read Classification Summary")
    plt.xticks(rotation=15)
    plt.tight_layout()
    #plt.show()
