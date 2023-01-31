#This part is for reading the coding len
import glob
LenDict = {}
file = open("../input/CodingLenInTemplate.txt")
for line in file:
    LenDict[line.strip().split()[0]] = int(line.strip().split()[1])
file.close()

sides = ["L","R"]

for side in sides:
    FileName2Inspect = glob.glob("../intermediates/" + side + "_alignment_*BeforeSelection_pruned.txt")[0]
    LINES = []
    CodingLen = LenDict[side]
    file = open(FileName2Inspect)
    OutLines = []
    for line in file:
        ReadID, qstart, qend, sstart,send, score, qlen, slen, qseq , sseq = line.strip("\n").split("\t")
        #MaxInd = max(0, CodingLen - int(sstart) + 1)
        if (int(qstart) >= int(qend)):
            continue
        #Aligner is penalized 2 units if it starts aligns one seqeunce with a gap, however it is not penalized if 
        # it aligns a T with a gap in between.
        MinInd = max(10 - int(sstart),0)
        if (qseq.lstrip("-")[MinInd:].count("-") == 0) and (sseq.lstrip("-")[MinInd:].count("-") == 0):
            LINES.append(line)
        else:
            OutLines.append(line)
    newFileNameForSelected = FileName2Inspect.replace("_BeforeSelection","")
    newFileNameForUnSelected = FileName2Inspect.replace("_BeforeSelection","_UnSelected")

    outFile1 = open(newFileNameForSelected,'w')
    outFile1.write("".join(LINES))
    outFile1.close()

    outFile2 = open(newFileNameForUnSelected,'w')
    outFile2.write("".join(OutLines))
    outFile2.close()
    file.close()
