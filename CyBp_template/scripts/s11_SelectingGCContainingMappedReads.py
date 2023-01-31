import pandas as pd

blastdb = pd.read_csv("../intermediates/BlastMatches.txt", sep="\t", header=None)

GenomeMatches = blastdb[[0]].drop_duplicates()

GenomeMatches.to_csv("../intermediates/GenomeMapped.txt", sep="\t", header=None, index=None)
