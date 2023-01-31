TailFile = "../output/TailFile.fa"

file = open(TailFile)
SeqDict={}
for line in file:
    if line[0]==">":
        continue
    SeqDict[line.strip()] = SeqDict.get(line.strip(),0) + 1

OutPut = "../output/TailCounts.txt"

outF = open(OutPut,'w')

for item in SeqDict.keys():
    outF.write(item + '\t' + str(SeqDict[item]) + '\n')

outF.close()

