## 1. PacBio Read Alignment and Post-Processing Pipeline 


This folder contains a 3-step Python workflow for automated preprocessing, alignment, and post-analysis of PacBio long-read sequencing data, as well as a sample dataset named `Lib_0_to_1000.fastq` located in the `data` folder


## Overview

| Step | Script | Purpose | Output |
|-|-|-|-|
| **1** | `11_generate_hpc_commands.py` | Pre-filters raw FASTQ files, generates `.fasta` files and companion pickle objects, and creates batch scripts for HPC alignment jobs. | Filtered FASTA and pickle files, plus corresponding HPC job scripts. |
| **2** | `12_Lib_0_to_1000_Seq.fastq_runEMBOSS_local.sh` | Executes EMBOSS `water` pairwise alignments for each filtered library. | A `.txt` file in EMBOSS format containing both forward and reverse alignments | 
| **3** | `13_parse_EMBOSS.py` | Parses alignment output files, determines best alignments, and constructs per-position nucleotide count matrices. | CSV files under `/scores/` | 

---

## Dependencies

- Python ≥ 3.8  
- [Biopython](https://biopython.org/)  
- NumPy ≥ 1.20  
- Pandas ≥ 1.3  
- EMBOSS ≥ 6.6.0  

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
python3 13_parse_EMBOSS.py
```


