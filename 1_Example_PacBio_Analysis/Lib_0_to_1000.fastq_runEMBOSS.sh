#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --ntasks=1 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8  
#SBATCH --mem=8G         
#SBATCH -J "Lib_0_to_1000.fastq_filter_reads.fasta"                                                                

module load emboss/6.6.0-gcc-11.3.1-ahzdg56

reference_sequence="reference.fa"
read_sequence="Lib_0_to_1000.fastq_filter_reads.fasta"
matrix="./CUSTOM"

water -asequence $reference_sequence -bsequence $read_sequence -datafile $matrix -gapopen 10 -gapextend 0.5 -outfile ${read_sequence}.txt
