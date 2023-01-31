import pandas as pd
import sys

df1 = pd.read_csv(sys.argv[1], sep="\t", header=None)
df2 = pd.read_csv(sys.argv[2], sep="\t", header=None)
df2 = df2[~df2[0].isin(df1[0])]

dfAll = pd.concat([df1,df2])
dfAll.sort_values(by=0, inplace=True)

dfAll.to_csv(sys.argv[3], sep="\t", header=None, index=None)
