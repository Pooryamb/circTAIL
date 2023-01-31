import numpy as np
import pandas as pd
import glob
np.random.seed(1000)
SampleSize = 100
file = open("../output/TailFile.fa")
Lines = file.readlines()
NumOfReads = len(Lines)//2

def FileContet2Dict(content):
    SeqDict = {}
    Segments = content.strip().split(">")[1:]
    for seg in Segments:
        SeqDict[seg.split("\n")[0]] = seg.split("\n")[1]
    return SeqDict

SelectedReads = []
for i in range(NumOfReads):
    if len(Lines[2*i + 1].strip()) >0:
        SelectedReads.append(int(Lines[2*i].strip().strip(">")))

if len(SelectedReads) <SampleSize:
    RandomReadIDs = SelectedReads
else:
    RandomReadIDs = np.random.choice(SelectedReads, size = SampleSize, replace=False)


df_Aln = pd.read_csv("../intermediates/R_alignment_reversed_pruned.txt", sep="\t", header=None)

ReadsFile = glob.glob("../intermediates/*_renamed_reversed.fasta")[0]
SeqDict = FileContet2Dict(open(ReadsFile).read())

temp = open(glob.glob("../intermediates/*_R_reversed.fa")[0]).readlines()[1].strip()

df_Aln_sel = df_Aln[df_Aln[0].isin(RandomReadIDs)]
sorted_IDs = df_Aln_sel[0].to_list()
#print(sorted_IDs)

LINES = []

for ID in sorted_IDs:
    R_sde = df_Aln_sel[df_Aln_sel[0] == ID]  
    LINES.append("The alignment for read "+ str(ID) +" starts\n\n")
    qstart = R_sde[1].to_list()[0]
    qend = R_sde[2].to_list()[0]
    sstart = R_sde[3].to_list()[0]
    send = R_sde[4].to_list()[0]
    qseq_aln = R_sde[8].to_list()[0]
    sseq_aln = R_sde[9].to_list()[0]
    qseq = SeqDict[str(ID)]
    LBA = max(qstart, sstart)  #This is for calculating length before alignment
    
    L0 = (LBA - qstart)*"-" +  qseq[:qstart-1] + qseq_aln + qseq[qend:]
    L1 = (LBA - 1)*" " + qseq_aln
    L2 = (LBA - 1)*" " + sseq_aln
    L3 = (LBA - sstart)*"-" +  temp[:sstart-1] + sseq_aln + temp[send:]
    L_max = max(len(L0), len(L1), len(L2), len(L3))
    
    L0 = L0 + (L_max - len(L0))* " "
    L1 = L1 + (L_max - len(L1))* " "
    L2 = L2 + (L_max - len(L2))* " "
    L3 = L3 + (L_max - len(L3))* " "
    
    LINES = LINES + [L0[::-1], L1[::-1], L2[::-1], L3[::-1],"\n"]


#R_sde = df_Aln_R_sel[df_Aln_R_sel[0] == ID]
VisualizedLeftSide = open("../output/VisualizedR_side_NonZero_LenReads.txt",'w')
VisualizedLeftSide.write("\n".join(LINES))
VisualizedLeftSide.close()


