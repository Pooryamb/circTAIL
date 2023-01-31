#qseq = 'TGTTTATTTTGTTTTGCATTTGTTGTTGGTATGTTTTTAGTGTATGCGTTTATTTTTG---TATA-'
#sseq = 'TGTTTATTTTGTTTTGCATTTGTTGTTGGTATGTTTTTAGTGTATGCGTTTT-TTTTGTTTTATAT'
#qseq = "TTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGGTTTTAA-A--------A------A-AA-" 
#sseq = "TTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTTAATATTTTTTTTATTTTTTATAAT"
#TTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTTATA-A--------A------A-AA-        
#TTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTTA-ATATTTTTTTTATTTTTTATAAT

import sys
import numpy as np
#file_name = sys.argv[1]
#file = open(sys.argv[1])


GapScore = 1
MismatchScore = 1
CGpenalty = -3

PruningMatrix = {'A':{'A':0, "C":1 , "G":1 , "T":1}, 
                 'C':{'A':1, "C":-3, "G":1 , "T":1},
                 'G':{'A':1, "C":1 , "G":-3, "T":1},
                 'T':{'A':1, "C":1 , "G":1 , "T":0}}

def Pruner(record):
    read_ID, qstart, qend, sstart, send, score, qlen, slen, qseq, sseq = record.strip("\n").split("\t")

    m = len(qseq)
    scores = m*[0] 
    CumScore = 0
    MaxScoreByFar = 0
    MaxInd = 10000
    Mismatches = m * [-1] 
    Matches = m * [-1]
    MismatchesByFar = 0
    MatchesByFar = 0
    
    for i in range(1,m+1):
        ind = -1*i
        if (qseq[ind] == "-") or (sseq[ind] == "-"):
            CumScore += GapScore
        else:
            CumScore += PruningMatrix[qseq[ind]][sseq[ind]]
            if (qseq[ind] == sseq[ind]):
                MatchesByFar += 1
            else:
                MismatchesByFar += 1
    
        scores[ind] = CumScore
        Mismatches[ind] = MismatchesByFar
        Matches[ind] = MatchesByFar
        
        if scores[ind] > MaxScoreByFar:
            MaxScoreByFar = scores[ind]
            MaxInd = ind
#    maxArg = len(scores) - np.argmax(scores[::-1]) -1
    if MaxInd == 10000:
        return record.strip("\n")
    Tnum = 0
    RemQ = qseq[MaxInd:].replace("-",'')
    RemS = sseq[MaxInd:].replace("-",'')
    for i in range(min(len(RemQ), len(RemS))):
        if (RemQ[i]==RemS[i]) and (RemS[i]=="T"):
            Tnum +=1
        else:
            break
    newQend = int(qend) -(-1*MaxInd - qseq[MaxInd:].count("-")) + Tnum
    newSend = int(send) -(-1*MaxInd - sseq[MaxInd:].count("-")) + Tnum
    newScore = int(score) - Matches[MaxInd] + 3 * Mismatches[MaxInd] + Tnum
    return "\t".join([read_ID, qstart, str(newQend), sstart, str(newSend), str(newScore), qlen, slen, qseq[:MaxInd] + Tnum*"T", sseq[:MaxInd]+ Tnum*"T"])
    
    
line = "280922  1       71      2       88      69      177     164     TTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTTATA-A--------A------A-AA-        TTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTTA-ATATTTTTTTTATTTTTTATAAT"


inputFile = sys.argv[1]
if len(inputFile.split(".")) > 1:
    outPutFile = ".".join(inputFile.split(".")[:-1]) +"_pruned." + inputFile.split(".")[-1]
else:
    outPutFile = inputFile + "_pruned"

file = open(inputFile)
LINES = []
for line in file:
    LINES.append(Pruner(line))
outFile = open(outPutFile,'w')
outFile.write("\n".join(LINES))
outFile.close()
