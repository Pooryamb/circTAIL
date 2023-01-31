from Fasta2Dict import Fasta2Dict


TailDict = Fasta2Dict(open("../output/TailFile.fa").read().strip())

GC = 0
GC1 = 0
GC2 = 0

for tail in TailDict.values():
    GCnum = tail.count("C") + tail.count("G")
    if GCnum >0:
        if GCnum==1:
            GC1 +=1
        elif GCnum==2:
            GC2 +=1
        GC +=1
file = open("../output/NumOfTailsWithGC.txt",'w')
file.write("total\t"+ str(len(TailDict)) + "\n"+ "GC\t" + str(GC) + "\n" + "GC1\t" + str(GC1) + "\n" + "GC2\t" + str(GC2) )