import os
import pickle
from Bio.Seq import Seq

"""
This script automates preprocessing of PacBio sequencing data in FASTQ format.

Workflow:
1) Discover input files by scanning a directory for names beginning with 'Lib_' and ending with '.fastq'.
2) For each file, read sequences and associated quality strings.
3) Filter reads to retain only those whose lengths fall within ±n bases of a given reference sequence length.
4) Store the filtered reads and their corresponding quality strings as pickled lists for downstream analysis.
5) Generate a bash script for submitting a HPC task for alignment analysis.

"""

sample_overhead_path = "./../data/"
sample_names = [
    file for file in os.listdir(sample_overhead_path)
    if file.endswith('.fastq') and file.startswith('Lib_')
]

reference_seq = "GCACTGAAGGTCCTCAATCGCACTGGAAACATCAAGGTCGNNNNNNNNNNNTTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGCaaagaggagaaaggatctatgGTCagtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattctctgtcagtggagagggtgaaggtgatgcaacatacggaaaacttacccttaaatttatttgcactactggaaagctacctgttccatggccaacacttgtcactactttgRStYRKggtVWtcaatgctttKcaagatacccagatcatatgaaacagcatgactttttcaagagtgccatgcccgaaggttatgtacaggaaagaactatattttacaaagatgacgggaactacaaatcacgtgctgaagtcaagtttgaaggtgataccctcgttaatagaRttgagttaaaaggtattgattttaaagaagatggaaacattcttggacacaaaatggaatacaacYataWctcacRcaatgtatacatcaYggcagacaaacaaaagaatggaatcMaagYtaacttcaaaRttagacacaacMttgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttccaaagatcccaaSgaaaagagagatcacatgatccttcttgagtttgtaacagctgctgggattacacatggcatggatgaactatacRaataaGTTTCCGTCTACGAACTCC"

n = 10

updated_file_paths = [f"{sample_overhead_path}{file_name}" for file_name in sample_names]

bash_template = """#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --ntasks=1 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8  
#SBATCH --mem=8G         
#SBATCH -J "{sample_name}"                                                                

module load emboss/6.6.0-gcc-11.3.1-ahzdg56

reference_sequence="reference.fa"
read_sequence="{sample_name}"
matrix="./CUSTOM"

water -asequence $reference_sequence -bsequence $read_sequence -datafile $matrix -gapopen 10 -gapextend 0.5 -outfile ${{read_sequence}}.txt
"""


def readFastq(filename):
    """Reads FASTQ file and removes special characters."""
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip()  # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def write_fasta(sequences, output_file):
    """Writes a FASTA file with reads and their reverse complements."""
    with open(output_file, 'w') as fasta_file:
        for i, seq in enumerate(sequences):
            header = f">Read_{i+1}"
            header_rc = f">Read_{i+1}_rc"
            seq_rc = str(Seq(seq).reverse_complement())
            fasta_file.write(f"{header}\n{seq}\n{header_rc}\n{seq_rc}\n")


for i in range(len(sample_names)):
    name = sample_names[i]
    path_name = updated_file_paths[i]
    sequences, qualities = readFastq(path_name)
    filtered_reads, seq_qualities = zip(*[
        (read, qual)
        for read, qual in zip(sequences, qualities)
        if len(reference_seq) - n <= len(read) <= len(reference_seq) + n
    ])
    filtered_reads = list(filtered_reads)
    seq_qualities = list(seq_qualities)

    seq_list_name = f"{name}_filtered_reads_list"
    score_list_name = f"{name}_quality_score_list"

    with open(seq_list_name, "wb") as fp:
        pickle.dump(filtered_reads, fp)

    with open(score_list_name, "wb") as fp:
        pickle.dump(seq_qualities, fp)

    output_file = f"{name}_filter_reads.fasta"
    write_fasta(filtered_reads, output_file)

    script_content = bash_template.format(sample_name=output_file)
    script_name = f"{name}_runEMBOSS.sh"
    with open(script_name, 'w') as script_file:
        script_file.write(script_content)
    # Make the script executable
    os.chmod(script_name, 0o755)
