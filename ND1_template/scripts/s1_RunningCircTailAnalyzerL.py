import os
import glob

PosFile = "../input/CodingLenInTemplate.txt"
FileLines = open(PosFile).read().strip().split("\n")

LocDict= {}
for line in FileLines:
    LocDict[line.split()[0]] = line.split()[1]

L_stop_codon = LocDict["L"]


ReadFile = glob.glob("../input/*.fasta")
L_temp = glob.glob("../input/*_L.fa")

os.system("./circTailAlignerGapLess.o " + ReadFile[0] + " " + L_temp[0] +  " 0 10 20 " +" > ../intermediates/L_alignment.txt")
