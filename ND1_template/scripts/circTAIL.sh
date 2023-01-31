#!/bin/bash

#SBATCH --time=00:90:00       # The duration in HH:MM:SS format
#SBATCH --cpus-per-task=1     # The number of cores
#SBATCH --mem=6000M          # Total memory for this task


module load StdEnv/2020 gcc/9.3.0
module load python/3.8.10
module load scipy-stack
module load blast+/2.12.0

set -e
touch StartingAnalysis

python s0_*
python s1_*
python s2_*
python s3_*

python s4_Pruning.py ../intermediates/L_alignment.txt
python s4_Pruning.py ../intermediates/R_alignment_reversed.txt

python s5_*
python s6_*
python s7_*
python s8_*


makeblastdb -in ../input/MaskedTranscriptome.fan -dbtype nucl -out ../intermediates/blastdb/Genome
blastn -task blastn -word_size 7 -query ../intermediates/BetUTRs.fa -evalue 0.001 -db ../intermediates/blastdb/Genome -outfmt 6 -out ../intermediates/BlastMatches.txt



python s9_*

python s10_*
python s11_*
python s12_*
python s13_*
python s14_*
python s15_*
python s16_*
python s17_*
python s18_*
python s19_*
python s20_*
touch EndingAnalysis
