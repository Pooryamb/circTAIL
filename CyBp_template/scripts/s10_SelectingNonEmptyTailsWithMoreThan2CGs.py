from Fasta2Dict import Fasta2Dict

FileContent = open("../intermediates/BetUTRs.fa").read()


FAdict = Fasta2Dict(FileContent)
MoreThan2CGs = [x for x in FAdict.keys() if ((FAdict[x].count("C") + FAdict[x].count("G")) >2 and len(FAdict[x]) > 1)]
outFile = open("../intermediates/BetUTRsNonEmptyMoreThan2CG.fa",'w')
LINES = []
for item in MoreThan2CGs:
    LINES.append(">" + item)
    LINES.append(FAdict[item])
outFile.write("\n".join(LINES))
outFile.close()