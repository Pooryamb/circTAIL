import glob
import os

InFile = glob.glob("../input/*.fasta")[0]
if len(glob.glob("../input/*_renamed.fasta")) != 0 :
    quit()

OutFile = InFile.replace(".fasta", "_renamed.fasta")
Lines = []
IDorSeq = 0 
ID_num = 0
FileLines = open(InFile).readlines()
for line in FileLines:
    if (IDorSeq%2 == 0):
        Lines.append(">" + str(ID_num))
        ID_num += 1
    else:
        Lines.append(line.strip())
    IDorSeq += 1


os.system("mv " + InFile + " " + InFile.replace(".fasta", ".fasta_main"))
outF = open(OutFile, 'w')

outF.write("\n".join(Lines))

outF.close()
