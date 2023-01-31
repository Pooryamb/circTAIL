#circTAIL analyzer

Prerequisites:

•	Python

•	gcc

•	Numpy

•	Pandas

•	Blast

For never-edited libraries and pre-edited libraries, the aligner in circTailAlignerGapLess 
and circTailAlignerZeroTGap must be compiled, respectively. The compiled aligners should be 
copied to the “scripts” directory of each library. 

The fasta file containing the merged reads must be placed in the input directory with the 
".fasta" extension. The masked transcriptome sequence should be inside the input directory 
named as "MaskedTranscriptome.fan".

The circTAIL analysis will start by executing the circTAIL.sh script.
