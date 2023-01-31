#!/bin/bash

#SBATCH --time=02:00:00       # The duration in HH:MM:SS format
#SBATCH --cpus-per-task=1     # The number of cores
#SBATCH --mem=6000M          # Total memory for this task


module load StdEnv/2020 gcc/9.3.0
module load python/3.8.10
module load scipy-stack
module load blast+/2.12.0

set -e
touch StartingAnalysis

dos2unix ../input/*
python s0_* ../input/*.fasta
./circTailAlignerGapLess.o ../input/*.fasta ../input/*_L.fa  0 10 20 > ../intermediates/L_alignment_BeforeSelection_BeforeTgapCorrection_BeforeOverLapCorrection.txt

python s2_*
./circTailAlignerZeroTGap.o ../intermediates/*_reversed.fasta ../intermediates/*_R_reversed.fa  0 10 20 > ../intermediates/R_alignment_reversed_BeforeSelection_BeforeTgapCorrection_BeforeOverLapCorrection.txt


python TgapCorrector.py ../intermediates/L_alignment_BeforeSelection_BeforeTgapCorrection_BeforeOverLapCorrection.txt ../intermediates/L_alignment_BeforeSelection_BeforeOverLapCorrection.txt 
python TgapCorrector.py ../intermediates/R_alignment_reversed_BeforeSelection_BeforeTgapCorrection_BeforeOverLapCorrection.txt ../intermediates/R_alignment_reversed_BeforeSelection_BeforeOverLapCorrection.txt


python s4_Pruning.py ../intermediates/L_alignment_BeforeSelection_BeforeOverLapCorrection.txt \
../intermediates/L_alignment_BeforeSelection_BeforeOverLapCorrection_pruned.txt

python s4_Pruning.py ../intermediates/R_alignment_reversed_BeforeSelection_BeforeOverLapCorrection.txt \
../intermediates/R_alignment_reversed_BeforeSelection_BeforeOverLapCorrection_pruned.txt

python OverlappingReadsFinder.py ../intermediates/L_alignment_BeforeSelection_BeforeOverLapCorrection_pruned.txt \
../intermediates/R_alignment_reversed_BeforeSelection_BeforeOverLapCorrection_pruned.txt \
../intermediates/ReadsWithOverlaps.txt ../intermediates/ReadsWithoutOverlapIDs.txt

python ConflictResolutionByAliScore.py ../intermediates/ReadsWithOverlaps.txt \
../intermediates/Overlapping_L_alignment_Overlap_Resolved.txt \
../intermediates/Overlapping_R_alignment_reversed_Overlap_Resolved.txt

#### Merging Overlap Resolved cases with the rest ####
python mergingOverLapResolvedWithTheRest.py ../intermediates/Overlapping_L_alignment_Overlap_Resolved.txt \
../intermediates/L_alignment_BeforeSelection_BeforeOverLapCorrection_pruned.txt \
../intermediates/L_alignment_BeforeSelection_pruned.txt

python mergingOverLapResolvedWithTheRest.py ../intermediates/Overlapping_R_alignment_reversed_Overlap_Resolved.txt \
../intermediates/R_alignment_reversed_BeforeSelection_BeforeOverLapCorrection_pruned.txt \
../intermediates/R_alignment_reversed_BeforeSelection_pruned.txt


python s5_*
python s6_*
python s7_*
python s8_*
python s9_*
python s10_*

makeblastdb -in ../input/MaskedTranscriptome.fan -dbtype nucl -out ../intermediates/blastdb/Genome
blastn -task blastn -strand plus -max_hsps 1 -max_target_seqs 1 -word_size 7 -query ../intermediates/BetUTRsNonEmptyMoreThan2CG.fa -evalue 0.001 -db ../intermediates/blastdb/Genome -outfmt 6 -out ../intermediates/BlastMatches.txt

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
python s21_*
touch EndingAnalysis
