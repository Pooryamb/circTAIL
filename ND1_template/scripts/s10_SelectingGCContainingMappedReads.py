import pandas as pd

blastdb = pd.read_csv("../intermediates/BlastMatches.txt", sep="\t", header=None)
HighGC = pd.read_csv("../intermediates/BetUTRsMoreThan2CG.txt", header=None, sep="\t")

GenomeMatches = blastdb[blastdb[0].isin(HighGC[0])][[0]].drop_duplicates()

GenomeMatches.to_csv("../intermediates/GenomeMapped.txt", sep="\t", header=None, index=None)
