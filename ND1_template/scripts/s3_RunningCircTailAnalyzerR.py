import os
import glob
PosFile = "../input/CodingLenInTemplate.txt"
FileLines = open(PosFile).read().strip().split("\n")

LocDict= {}
for line in FileLines:
    LocDict[line.split()[0]] = line.split()[1]

R_stop_codon = LocDict["R"]

ReadFile = glob.glob("../intermediates/*_reversed.fasta")
R_temp = glob.glob("../intermediates/*_R_reversed.fa")

os.system("./circTailAlignerGapLess.o " + ReadFile[0] + " " + R_temp[0] +  " -1000 1000 1000 " +" > ../intermediates/R_alignment_reversed.txt")
