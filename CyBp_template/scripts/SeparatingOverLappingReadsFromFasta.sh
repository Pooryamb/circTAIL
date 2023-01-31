#module load seqtk

seqtk subseq ../input/*_renamed.fasta ../intermediates/IDsOfOverlappingReads.txt > ../intermediates/ReadWithOverLap.fasta
