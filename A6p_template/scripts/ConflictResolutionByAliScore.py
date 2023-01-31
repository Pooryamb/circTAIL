from s4_Pruning import Pruner
import numpy as np
import sys

CodonFileLines = open("../input/CodingLenInTemplate.txt").readlines()
L_Len_cod = int(CodonFileLines[0].strip().split()[1])
R_Len_cod = int(CodonFileLines[1].strip().split()[1])

def ConvertRecord2Nums(record):
    Read_id, Q_start, Q_end, T_start, T_end, score, Qlen, Tlen, Qseq, Tseq = record.strip("\n").split("\t")
    return [Read_id, int(Q_start), int(Q_end), int(T_start), int(T_end), int(score), int(Qlen), int(Tlen), Qseq, Tseq]
    
def NoOfLets2ReadFromEnd(read, overlap):
    NoOfNonGapNucs = 0
    NumOfNucs2ReadFromEnd = 0
    while NoOfNonGapNucs != overlap:
        NumOfNucs2ReadFromEnd += 1
        if read[-NumOfNucs2ReadFromEnd]!= "-":
            NoOfNonGapNucs += 1
    return NumOfNucs2ReadFromEnd

def RemovingNucsFromRecord(record, NumOfNucs2Remove):
    if NumOfNucs2Remove==0:
        return record
    Read_id, Q_start, Q_end, T_start, T_end, score, Qlen, Tlen, Qseq, Tseq = record.strip("\n").split("\t")
    GapsInQ = Qseq[-NumOfNucs2Remove:].count("-")
    GapsInT = Tseq[-NumOfNucs2Remove:].count("-")
    Q_end = str(int(Q_end) - NumOfNucs2Remove + GapsInQ)
    T_end = str(int(T_end) - NumOfNucs2Remove + GapsInT)    
    ScoreReduction = 0
    for i in range(1,NumOfNucs2Remove +1):
        if (Qseq[-i] == Tseq[-i]):
            ScoreReduction += 1
        elif (Qseq[-i] == "-" or Tseq[-i]=="-"):
            continue
        else:
            ScoreReduction -= 3
    score = str(int(score) - ScoreReduction)
    Qseq = Qseq[:-NumOfNucs2Remove]
    Tseq = Tseq[:-NumOfNucs2Remove]
    return "\t".join([Read_id, Q_start, Q_end, T_start, T_end, score, Qlen, Tlen, Qseq, Tseq])
    
def ScoreVecCalculator(record, PrunedRecord, overlap):
    RecData = ConvertRecord2Nums(record)
    VecLen = ConvertRecord2Nums(record)[2] - ConvertRecord2Nums(PrunedRecord)[2]
    NoOfLetter2Check = NoOfLets2ReadFromEnd(ConvertRecord2Nums(record)[8], VecLen)
    qseq_sel = RecData[8][-NoOfLetter2Check:]
    tseq_sel = RecData[9][-NoOfLetter2Check:]
    ScoreVec = [0] * VecLen
    CumScore = 0
    ScoreVecCurser = 0
    NoOfNonGaps = 0
    for i in range(NoOfLetter2Check):
        if qseq_sel[i] == tseq_sel[i]:
            CumScore += 1
            ScoreVec[ScoreVecCurser] = CumScore
            ScoreVecCurser += 1
        elif (qseq_sel[i] =="-"):
            CumScore += 0  
        elif (tseq_sel[i] =="-"):
            CumScore += 0
            ScoreVec[ScoreVecCurser] = CumScore
            ScoreVecCurser += 1
        else: #qseq_sel[i] != tseq_sel[i] (and none of them is gap)
            CumScore -= 3
            ScoreVec[ScoreVecCurser] = CumScore
            ScoreVecCurser += 1            
    return ScoreVec[-overlap:]


def BorderFinder(VecL, VecR):
    VecL = [0] + VecL
    VecR = [0] + VecR
    VecR = VecR[::-1]
    TotalScoreVec = np.array(VecL) + np.array(VecR)
    L_side = TotalScoreVec.argmax()
    R_side = len(VecR) -1 - L_side #I added one unit to it in the beginning, so I remove it.
    #first number shows the length that must be removed from R and the second one the length to be removed from L
    return L_side, R_side
    


L_aln_lines = []
R_aln_lines = []
file = open(sys.argv[1])
for line in file:    
    Read_id, LQ_start, LQ_end, LT_start, LT_end, Lscore, QlenL, TlenL, LQseq, LTseq, RQ_start, RQ_end, RT_start, RT_end, Rscore, QlenR, TlenR, RQseq, RTseq = line.strip("\n").split("\t")
    L_record = "\t".join([Read_id, LQ_start, LQ_end, LT_start, LT_end, Lscore, QlenL, TlenL, LQseq, LTseq])
    R_record = "\t".join([Read_id, RQ_start, RQ_end, RT_start, RT_end, Rscore, QlenR, TlenR, RQseq, RTseq])        
    overlap = int(LQ_end) + int(RQ_end)- int(QlenL)

    if (len(RQseq) < R_Len_cod or len(LQseq) < L_Len_cod):
        continue    
    NumOfNucs2ReadFromEndR = NoOfLets2ReadFromEnd(RQseq,overlap)
    NumOfNucs2ReadFromEndL = NoOfLets2ReadFromEnd(LQseq,overlap)

    L_record_pruned = Pruner(RemovingNucsFromRecord(L_record, NumOfNucs2ReadFromEndL))
    R_record_pruned = Pruner(RemovingNucsFromRecord(R_record, NumOfNucs2ReadFromEndR))
    #I have to prune until the start of the alignment. And then find the number of corrections needed t

    Rvec = ScoreVecCalculator(R_record, R_record_pruned, overlap)
    Lvec = ScoreVecCalculator(L_record, L_record_pruned, overlap)

    Len2RemoveR, Len2RemoveL = BorderFinder(Lvec, Rvec)

    L_rec_final = Pruner(RemovingNucsFromRecord(L_record, Len2RemoveL))
    R_rec_final = Pruner(RemovingNucsFromRecord(R_record, Len2RemoveR))
    L_aln_lines.append(L_rec_final)
    R_aln_lines.append(R_rec_final)
    
L_aln_corrected = open(sys.argv[2],'w')
L_aln_corrected.write("\n".join(L_aln_lines))
L_aln_corrected.close()

R_aln_corrected = open(sys.argv[3], 'w')
R_aln_corrected.write("\n".join(R_aln_lines))
R_aln_corrected.close()
