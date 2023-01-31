from Fasta2Dict import Fasta2Dict

FileContent = open("../intermediates/BetUTRs.fa").read()


FAdict = Fasta2Dict(FileContent)
MoreThan2CGs = [x for x in FAdict.keys() if (FAdict[x].count("C") + FAdict[x].count("G")) >2]
outFile = open("../intermediates/BetUTRsMoreThan2CG.txt",'w')

outFile.write("\n".join(MoreThan2CGs))
outFile.close()