## 1. PacBio Read Alignment and Post-Processing Pipeline 


This folder contains a 3-step Python workflow for automated preprocessing, alignment, and post-analysis of PacBio long-read sequencing data, as well as a sample dataset named `Lib_0_to_1000.fastq` located in the `data` folder


## Overview

| Step | Script | Purpose | Output |
|-|-|-|-|
| **1** | `11_generate_hpc_commands.py` | Pre-filters raw FASTQ files, generates `.fasta` files and companion pickle objects, and creates batch scripts for HPC alignment jobs. | Filtered FASTA and pickle files, plus corresponding HPC job scripts. |
| **2** | `12_Lib_0_to_1000_Seq.fastq_runEMBOSS_local.sh` | Executes EMBOSS `water` pairwise alignments for each filtered library. | A `.txt` file in EMBOSS format containing both forward and reverse alignments | 
| **3** | `13_parse_EMBOSS.py` | Parses alignment output files, determines best alignments, and constructs per-position nucleotide count matrices. | CSV files under `/scores/` | 

---

## Example Run

```bash
# Step 1: Preprocess and generate jobs
python 11_generate_hpc_commands.py

# Step 2: Submit generated alignment jobs
sbatch Lib_0_to_1000.fastq_runEMBOSS.sh

## OR

bash 12_Lib_0_to_1000_Seq.fastq_runEMBOSS_local.sh

# Step 3: Parse alignments and build matrices
python 13_parse_EMBOSS.py
```

## 2. Sidewinder Junction Analysis and Misassembly Detection

This folder contains a Python workflow for **automated analysis of Sidewinder junctions** in assembly product sequencing reads, with an example by the ** mScarlet ** sample mentioned in the paper. 
## 1. Junction Analysis and Misassembly Detection Pipeline

This folder contains a Python workflow for **automated analysis of fragment junctions** in DNA assembly products. 

## Overview

| Step | Script | Purpose | Output |
|-|-|-|-|
| **1** | `21_junction_code.py` | Generates all possible Sidewinder junctions, performs BLAST alignment against sequencing reads, visualizes with heatmap, outputs misassembled reads. | `mScar_heatmap.png`, `mScar_heatmap.xlsx`, and `misassemblies.fasta` |

## Example Run

```
python 21_junction_analysis.py
```

## 3. Fragment-to-Read Assembly Validation Pipeline

This folder contains a Python workflow for conducting fragment-level analysis by aligning expected fragments against sequencing reads (e.g., from PacBio or Nanopore) and evaluating assembly accuracy.

---

## Overview

| Step | Script | Purpose | Output |
|-|-|-|-|
| **1** | `31_Fragment_Analysis.py` | Converts reads from FASTQ → FASTA, writes fragment FASTA, runs BLAST, parses hits, and classifies reads as full, partial, to-be-checked, or unusable. | `read_classification.csv`, `read_classification_summary.csv`, and `Sequences_to_be_checked.fasta` |

---

## Example Run

```
python 31_Fragment_Analysis.py
```

---

## Dependencies

- Python ≥ 3.8  
- Biopython ≥ 1.80
- NumPy ≥ 1.20  
- Pandas ≥ 1.3  
- EMBOSS ≥ 6.6.0
- Matplotlib ≥ 3.5
- Seaborn ≥ 0.11
- BLAST+ ≥ 2.12.0


