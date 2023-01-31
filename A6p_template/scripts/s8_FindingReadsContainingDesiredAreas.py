import pandas as pd
import os

dfL = pd.read_csv("../intermediates/L_alignment_pruned.txt", sep="\t", header=None)
dfR = pd.read_csv("../intermediates/R_alignment_pruned.txt", sep="\t", header=None)
dfL_Sel = dfL[[0,1,2,3,4]]
dfR_Sel = dfR[[0,1,2,3,4]]
df_BothSides = dfL_Sel.merge(dfR_Sel, on=[0])
df_BothSides.to_csv("../intermediates/Align_R_L_summary.txt",sep="\t", header=None, index=None)


CovFile=open("../intermediates/Coverage.txt")
CovLines = CovFile.readlines()
CovFile.close()

L_Min = int(CovLines[-1].strip().split("\t")[1])
L_Max = int(CovLines[-1].strip().split("\t")[2])
R_Min = int(CovLines[-2].strip().split("\t")[1])
R_Max = int(CovLines[-2].strip().split("\t")[2])


ReadsList = open("../intermediates/Align_R_L_summary.txt")
WellConsReads = open("../intermediates/WellConservedReadsAndUTRLocs.txt",'w')

for line in ReadsList:
    parts = line.strip().split("\t")
    ReadName = parts[0]
    L_Q_start,L_Q_end,L_T_start,L_T_end,R_Q_start,R_Q_end,R_T_start,R_T_end = [int(x) for x in parts[1:]]
    if ((L_T_start <= L_Min) and ( L_T_end >= L_Max)) and ((R_T_start <= R_Min) and ( R_T_end >= R_Max)):
        WellConsReads.write("\t".join([ReadName, str(L_Q_end), str(R_Q_start), str(L_T_end),str(R_T_start)]) + "\n")

WellConsReads.close()
ReadsList.close()
