import os
import pandas as pd
import glob

Start_Stop_CodonPosFile="../input/CodingLenInTemplate.txt"
TailPosFile="../intermediates/WellConservedReadsAndUTRLocs.txt"


GenomeMappedReads = open("../intermediates/GenomeMapped.txt").read().strip().split()

lines = open(Start_Stop_CodonPosFile).readlines()
LSide=int(lines[0].strip().split()[1])
RSide=int(lines[1].strip().split()[1])

LENGTH = {}
side = "R"
#add input to the next line
Ref = glob.glob("../input/*_" + side + ".fa")[0]

#add input to the next line
file = open(Ref)
content = file.readlines()
file.close()
LENGTH[side] = len("".join(content[1:]).replace("\n","").strip())

LocsOnAllReads= pd.read_csv(TailPosFile,sep="\t", header=None)
LocsOnAllReads = LocsOnAllReads[~LocsOnAllReads[0].isin(GenomeMappedReads)]

TailFull = LocsOnAllReads[LocsOnAllReads[2]-LocsOnAllReads[1] > 1]
TailLess = LocsOnAllReads[LocsOnAllReads[2]-LocsOnAllReads[1] <= 1]
DictOfPos = {"TailFull": TailFull, "TailLess": TailLess}

for key, TailPos in DictOfPos.items():
    TailPos[5] = TailPos[3]- LSide
    L_len = TailPos[5]
    L_Len = pd.DataFrame({0:L_len})
    Res_L = L_Len.groupby(0)[0].count()
    Res_L.to_csv("../output/L_len" + key +".csv",sep="\t",header=None)

    TailPos[6] = LENGTH["R"] - TailPos[4] - RSide +1
    R_len = TailPos[6]
    R_Len = pd.DataFrame({0:R_len})
    Res_R = R_Len.groupby(0)[0].count()
    Res_R.to_csv("../output/R_len" + key +".csv",sep="\t",header=None)



