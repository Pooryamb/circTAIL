import pandas as pd
from Fasta2Dict import Fasta2Dict

GenomeMappeds = open("../intermediates/GenomeMapped.txt").read().strip().split()


BetMappedRegions = "../intermediates/BetUTRs.fa"

SeqDict = Fasta2Dict(open(BetMappedRegions).read().strip())

CandidateIDs = set(SeqDict.keys()) - set(GenomeMappeds)

GenomeFreeTails = open("../output/TailFile.fa", 'w')

for ID in CandidateIDs:
    GenomeFreeTails.write(">" + ID + "\n" +SeqDict[ID] + "\n" )

GenomeFreeTails.close()




