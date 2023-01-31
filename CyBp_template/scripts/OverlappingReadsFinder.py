import numpy as np
import pandas as pd
import sys



#dfL = pd.read_csv("../intermediates/L_alignment_BeforeSelection_pruned.txt", sep="\t", header=None)
dfL = pd.read_csv(sys.argv[1], sep="\t", header=None)
#dfR = pd.read_csv("../intermediates/R_alignment_reversed_BeforeSelection_pruned.txt", sep="\t", header=None)
dfR = pd.read_csv(sys.argv[2], sep="\t", header=None)

dfAll = dfL.merge(dfR, on=0)


OverlappingDF    = dfAll[dfAll['2_x'] + dfAll['2_y'] >  dfAll['6_x']]
OverlappingDF.to_csv("../intermediates/ReadsWithOverlaps.txt", sep="\t", header=None, index =None)
OverlappingDF.to_csv(sys.argv[3], sep="\t", header=None, index =None)

NonOverlappingDF_IDs = dfAll[dfAll['2_x'] + dfAll['2_y'] <= dfAll['6_x']][[0]]
NonOverlappingDF_IDs.to_csv("../intermediates/ReadsWithoutOverlapIDs.txt", sep="\t", header=None, index =None)
NonOverlappingDF_IDs.to_csv(sys.argv[4], sep="\t", header=None, index =None)
