import pandas as pd


GenomeMappeds = open("../intermediates/GenomeMapped.txt").read().strip().split()


BetMappedRegions = "../intermediates/BetUTRs.fa"
SeqSegments = open(BetMappedRegions).read().strip().split(">")[1:]
DictOfSeq = {}

GenomeFreeTails = open("../output/TailFile.fa", 'w')

for Seg in SeqSegments:
    Read_ID = Seg.strip().split("\n")[0]
    if Read_ID not in GenomeMappeds:
        GenomeFreeTails.write(">" + Seg)

GenomeFreeTails.close()




