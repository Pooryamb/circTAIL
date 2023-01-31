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


GapScore = 2



PruningMatrix = {'A':{'A':0, "C":3 , "G":3 , "T":3}, 
                 'C':{'A':3, "C":-3, "G":3 , "T":3},
                 'G':{'A':3, "C":3 , "G":-3, "T":3},
                 'T':{'A':3, "C":3 , "G":3 , "T":0}}

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
    
    
#line = "280922  1       71      2       88      69      177     164     TTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTTATA-A--------A------A-AA-        TTTACCCGGTGTAAAGTATTATACACGTATTGTAAGTTAGATTTAGATATAAGATATGTTTTTA-ATATTTTTTTTATTTTTTATAAT"

if __name__=="__main__":
    inputFile = sys.argv[1]
    
    
    file = open(inputFile)
    LINES = []
    for line in file:
        LINES.append(Pruner(line))
    
    outputFile = sys.argv[2]
    outFile = open(sys.argv[2],'w')
    outFile.write("\n".join(LINES))
    outFile.close()
