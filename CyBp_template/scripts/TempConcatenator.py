import glob
import os
def TempSeqExtractor(fileAdd):
    return "".join(open(fileAdd).read().strip().split("\n")[1:])

L_temp_add = glob.glob("../input/*_L.fa")[0]
R_temp_add = glob.glob("../input/*_R.fa")[0]

L_Cont = TempSeqExtractor(L_temp_add)
R_Cont = TempSeqExtractor(R_temp_add)

ConCat = L_Cont + R_Cont

outFile = open("../intermediates/ConCatenatedTemps.fa",'w')
outFile.write(">Temp\n" + ConCat + "\n")
outFile.close()

os.system("needleall -auto -stdout -asequence ../intermediates/ReadWithOverLap.fasta -bsequence ../intermediates/ConCatenatedTemps.fa -datafile EDNAFULL -gapopen 10.0 -gapextend 0.5 -aformat3 fasta -snucleotide1 -snucleotide2 -outfile ../intermediates/ConcatenatedAlignment.fan")