import re
import sys

def Corrector(seq1, seq2):
    hits = list(re.finditer(r"-+T+", seq1))
    TempSeq1 = seq1
    for hit in hits:
        SpanOfHit = hit.span()
        NumberOfGaps = hit.group().count("-")
        TotalNumberOfTs = hit.group().count("T")
        StartOnSeq2 = NumberOfGaps + SpanOfHit[0]
        InspectedAreaOnSeq2 = seq2[StartOnSeq2:SpanOfHit[1]]
        NoOfTsMatching2Ts = len(InspectedAreaOnSeq2) - len(InspectedAreaOnSeq2.lstrip("T"))
        TempSeq1 = TempSeq1[:SpanOfHit[0]] + NoOfTsMatching2Ts * "T" + NumberOfGaps*"-" + (TotalNumberOfTs - NoOfTsMatching2Ts)*"T" + TempSeq1[SpanOfHit[1]:]
    return TempSeq1



#Files = ["../intermediates/L_alignment_BeforeSelection_BeforeTgapCorrection.txt", "../intermediates/R_alignment_reversed_BeforeSelection_BeforeTgapCorrection.txt"] #, "../intermediates/R_aln_ov_rev.txt"]

if __name__ == "__main__":
    fileName = sys.argv[1]
    file = open(fileName)
    OutLines = []
    for line in file:
        items = line.strip("\n").split("\t")
        if len(items)==0:
            continue
        if len(items) < 10:
            print(line)
        CorrectedReadSeq = Corrector(items[8], items[9])
        CorrectedTempSeq = Corrector(items[9], items[8])
        #if items[8] != CorrectedReadSeq or items[9]!= CorrectedTempSeq:
         #   print(line)
        newLine = "\t".join(items[:8] + [Corrector(items[8], items[9]), Corrector(items[9], items[8])])
        OutLines.append(newLine)
    outFile = open(sys.argv[2], 'w')
    outFile.write("\n".join(OutLines))
    outFile.close()

