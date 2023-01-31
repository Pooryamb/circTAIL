import os
import glob


ReadFileName = glob.glob("../input/*.fasta")[0]

ReadFile = open(ReadFileName)
Sequences = ReadFile.read().strip().split(">")[1:]
ReadFile.close()

SeqDict= {}
for Seq in Sequences:
    ReadName = Seq.strip().split("\n")[0]
    ReadSeq = "".join(Seq.strip().split("\n")[1:])
    SeqDict[ReadName] = ReadSeq
del Sequences

TailFile = open("../intermediates/BetUTRs.fa",'w')

LocationFile = open("../intermediates/WellConservedReadsAndUTRLocs.txt")

for line in LocationFile:
    ReadName,start,end,T_start,T_end = line.strip().split("\t")
    start = int(start)
    end = int(end)
    TailFile.write(">" + ReadName + "\n")
    TailFile.write(SeqDict[ReadName][start:end-1] + "\n")

TailFile.close()
LocationFile.close()