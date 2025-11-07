"""
Batch post-processing of EMBOSS 'water' alignments to build nucleotide and insertion
profile matrices from high-quality reads.

Workflow:
1) Discover EMBOSS output files: scans '../' for files named 'Lib_*.txt' that are
   the paired forward/reverse alignments produced by 'water'.
2) Derive companion paths for each file's original FASTQ-derived objects:
   - '<...>filtered_reads_list' (list of read strings)
   - '<...>quality_score_list' (list of ASCII quality strings)
3) Parse EMBOSS files:
   - Read per-alignment scores from '# Score:' lines.
   - Iterate alignments in pairs (forward, reverse) to select the better orientation
     with a custom threshold score > 3000 and both ref/read gaps < 20.
4) Build two position-wise matrices:
   - profile_matrix: counts of A/T/C/G and gaps ('-') observed at each reference position
     (only when read base has a corresponding 'I' quality char at that position,
      indicating a high-confidence base in this pipeline’s convention).
   - insertion_matrix: counts of A/T/C/G and gaps for insertions relative to the reference.
   Reads are kept only if their mean Phred score ≥ 39.5.
5) Save the matrices to CSV.

Note: 
start_gap = 9: Read positions before this are trimmed off when building matrices. This is a dataset-specific offset that can be changed.
"""

import os
import re
import pickle
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq

# Discover EMBOSS output files and derive companion object paths
sample_overhead_path = "./"
txt_files = [
    file for file in os.listdir(sample_overhead_path)
    if file.endswith('.txt') and file.startswith('Lib_')
]
updated_txt_paths = [f"{sample_overhead_path}{file_name}" for file_name in txt_files]
updated_fastq_read_paths = [
    file.replace('filter_reads.fasta.txt', 'filtered_reads_list')
    for file in updated_txt_paths
]
updated_score_paths = [
    file.replace('filter_reads.fasta.txt', 'quality_score_list')
    for file in updated_txt_paths
]

# Reference used to position profile matrices
reference_sequence = "GCACTGAAGGTCCTCAATCGCACTGGAAACATCAAGGTCGNNNNNNNNNNNTTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGCaaagaggagaaaggatctatgGTCagtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattctctgtcagtggagagggtgaaggtgatgcaacatacggaaaacttacccttaaatttatttgcactactggaaagctacctgttccatggccaacacttgtcactactttgRStYRKggtVWtcaatgctttKcaagatacccagatcatatgaaacagcatgactttttcaagagtgccatgcccgaaggttatgtacaggaaagaactatattttacaaagatgacgggaactacaaatcacgtgctgaagtcaagtttgaaggtgataccctcgttaatagaRttgagttaaaaggtattgattttaaagaagatggaaacattcttggacacaaaatggaatacaacYataWctcacRcaatgtatacatcaYggcagacaaacaaaagaatggaatcMaagYtaacttcaaaRttagacacaacMttgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttccaaagatcccaaSgaaaagagagatcacatgatccttcttgagtttgtaacagctgctgggattacacatggcatggatgaactatacRaataaGTTTCCGTCTACGAACTCC".upper()


def circular_index(string, sub):
    """Find the index of 'sub' in 'string' treating 'string' as circular; return -1 if absent."""
    index = string.find(sub)
    if index == -1:
        return (string + string).find(sub)
    return index


def extract_scores_from_water(file_path):
    """Extract float alignment scores from lines starting with '# Score:' in an EMBOSS file."""
    scores = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("# Score:"):
                score = float(line.split(":")[1].strip())
                scores.append(score)
    return scores


def parse_emboss(output_file):
    """
    Parse an EMBOSS 'water' output file in 'emboss' format.

    Returns:
        best_alignments_ref (list[str]): chosen ref alignments (uppercased) per read pair
        best_alignments_read (list[str]): chosen read alignments (uppercased) per read pair
        forward_indicator (list[bool]): True if forward chosen, False if reverse or 'wrong'
    """
    best_alignments_ref, best_alignments_read, forward_indicator = [], [], []
    current_pair = []
    emboss_scores = extract_scores_from_water(output_file)
    emboss_score_counter = 0

    print("Parsing", output_file)

    # Iterate in pairs: (forward, reverse)
    for alignment in AlignIO.parse(output_file, "emboss"):
        current_pair.append(alignment)
        if len(current_pair) == 2:
            forward_alignment, reverse_alignment = current_pair
            forward_score = emboss_scores[emboss_score_counter]
            reverse_score = emboss_scores[emboss_score_counter + 1]
            for_ref = str(forward_alignment[0].seq)
            for_read = str(forward_alignment[1].seq)
            rev_ref = str(reverse_alignment[0].seq)
            rev_read = str(reverse_alignment[1].seq)

            if (forward_score > reverse_score and forward_score > 3000 and
                    for_ref.count('-') < 20 and for_read.count('-') < 20):
                # forward read selected
                best_alignments_ref.append(for_ref.upper())
                best_alignments_read.append(for_read.upper())
                forward_indicator.append(True)

            elif (reverse_score > forward_score and reverse_score > 3000 and
                  rev_ref.count('-') < 20 and rev_read.count('-') < 20):
                # reverse read selected
                best_alignments_ref.append(rev_ref.upper())
                best_alignments_read.append(rev_read.upper())
                forward_indicator.append(False)

            else:
                # low-quality pair
                best_alignments_ref.append("wrong")
                best_alignments_read.append("wrong")
                forward_indicator.append(False)

            current_pair = []
            emboss_score_counter += 2

    return [best_alignments_ref, best_alignments_read, forward_indicator]


def ascii_to_phred(quality_string):
    """Convert an ASCII-encoded Phred quality string to a list of integer scores."""
    return [ord(char) - 33 for char in quality_string]


def make_profile_matrix(ref_aligned, read_aligned, reference_sequence, start_gap,
                        forward_indicator, score_list, fastq_read_list):
    """
    Build per-position nucleotide and insertion count matrices from chosen alignments.

    Args:
        ref_aligned (list[str]): aligned reference strings (with '-' gaps) per read.
        read_aligned (list[str]): aligned read strings (with '-' gaps) per read.
        reference_sequence (str): ungapped reference sequence.
        start_gap (int): number of leading reference positions to skip (trim) in output.
        forward_indicator (list[bool]): True if read maps forward, False if reverse.
        score_list (list[str]): ASCII quality strings (same order as reads).
        fastq_read_list (list[str]): raw read strings (same order).

    Returns:
        (pd.DataFrame, pd.DataFrame): (profile_matrix, insertion_matrix),
        both with columns ["A", "T", "C", "G", "-"] and rows across reference positions.
    """
    aligned_length = len(reference_sequence) - start_gap
    base_idx = {"A": 0, "T": 1, "C": 2, "G": 3, "-": 4}
    profile_matrix = np.zeros((aligned_length, 5), dtype=int)
    insertion_matrix = np.zeros((aligned_length, 5), dtype=int)

    for i in range(len(ref_aligned)):
        if (i % 5000) == 0:
            print("Processing Read #", i)

        ref = ref_aligned[i]
        read = read_aligned[i]
        isforward = forward_indicator[i]
        score = score_list[i]
        fastq_read = fastq_read_list[i]

        # Keep only high-quality reads (mean Phred >= 39.5)
        if np.mean(ascii_to_phred(score)) < 39.5:
            continue

        # Orient to forward if needed
        if isforward is False:
            fastq_read = str(Seq(fastq_read).reverse_complement())
            score = score[::-1]

        aligned_read_no_gap = read.replace("-", "").upper()
        aligned_ref_no_gap = ref.replace("-", "").upper()

        # Locate alignment segment in circular read and reference
        index_startOf_read_alignment = circular_index(fastq_read, aligned_read_no_gap)
        if index_startOf_read_alignment == -1 :
            continue

        fastq_read_substring = fastq_read[int(index_startOf_read_alignment):]
        fastq_score_substring = score[int(index_startOf_read_alignment):]

        ref_arg = circular_index(reference_sequence, aligned_ref_no_gap[:10])
        if index_startOf_read_alignment != -1:
            real_position = ref_arg - start_gap
            count_pointer = 0

            # Walk the pairwise alignment, updating counts
            for ref_char, read_char in zip(ref, read):
                if ref_char != '-':
                    # Base aligned to reference position
                    if read_char == "-":
                        profile_matrix[real_position, base_idx.get(read_char, 4)] += 1
                    else:
                        # Only count when ASCII char indicates 'I' (pipeline-specific gate)
                        if fastq_score_substring[count_pointer] == "I":
                            profile_matrix[real_position, base_idx.get(read_char, 4)] += 1
                            count_pointer += 1
                        else:
                            count_pointer += 1
                    real_position += 1
                else:
                    # Insertion relative to reference
                    if fastq_score_substring[count_pointer] == "I":
                        insertion_matrix[real_position, base_idx.get(read_char, 4)] += 1
                    count_pointer += 1

    df = pd.DataFrame(profile_matrix, columns=["A", "T", "C", "G", "-"])
    df_insert = pd.DataFrame(insertion_matrix, columns=["A", "T", "C", "G", "-"])
    return df, df_insert


def process_file(i):
    """
    Process a single EMBOSS file by index:
    - Load companion score and read lists
    - Parse alignments and construct profile matrices
    - Save matrices to CSV under 'scores/'
    """
    txt_file = txt_files[i]
    txt_file_path = updated_txt_paths[i]
    score_file = updated_score_paths[i]
    fastq_read_file = updated_fastq_read_paths[i]

    # Ensure required companion objects exist
    if (not os.path.exists(txt_file_path) or
            not os.path.exists(score_file) or
            not os.path.exists(fastq_read_file)):
        print(f"One of the files for {txt_file} is missing, skipping...")
        return

    with open(score_file, "rb") as fp:
        score_list = pickle.load(fp)
    with open(fastq_read_file, "rb") as fpp:
        fastq_read_list = pickle.load(fpp)

    # Extract range key like '123_to_456' from filename
    match = re.search(r'Lib_(\d+_to_\d+)\.fastq_filter_reads\.fasta\.txt', txt_file)
    if match:
        key = match.group(1)
        best_alignments_ref, best_alignments_read, forward_indicator = parse_emboss(txt_file_path)

        # Build matrices (start_gap = 9 as in original)
        df, df_insert = make_profile_matrix(
            best_alignments_ref,
            best_alignments_read,
            reference_sequence,
            9,
            forward_indicator,
            score_list,
            fastq_read_list
        )

        # Persist outputs
        os.makedirs("scores", exist_ok=True)
        df.to_csv(f"scores/{key}_Q40.csv", index=False)
        df_insert.to_csv(f"scores/{key}_insertion.csv", index=False)

    return f"Finished processing {txt_file}"


if __name__ == "__main__":
    # Parallelize over files; adjust to your HPC resources
    num_processes = 16
    with Pool(processes=num_processes) as pool:
        results = pool.map(process_file, range(len(updated_txt_paths)))

    for result in results:
        print(result)
