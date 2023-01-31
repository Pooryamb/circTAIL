import glob
L_tempName = glob.glob("../input/*_L.fa")[0]
L_len = len("".join(open(L_tempName).read().strip().split("\n")[1:]))


AlnRawFile = "../intermediates/ConcatenatedAlignment.fan"

AlnData = open(AlnRawFile).read().strip()

Parts = AlnData.split(">")[1:]
IDs = [x.split("\n")[0].strip() for x in Parts][0::2]
ReadSeq = ["".join(x.split("\n")[1:]) for x in Parts][0::2]
TempSeq = ["".join(x.split("\n")[1:]) for x in Parts][1::2]

DictOfAln = {a:(b,c) for (a,b,c) in zip(IDs, ReadSeq, TempSeq)}

def LocationOfNthNuc(read, loc):
    NoOfNonGapNucs = 0
    NumOfLets2Read = 0
    while NoOfNonGapNucs != loc:        
        if read[NumOfLets2Read]!= "-":
            NoOfNonGapNucs += 1
        NumOfLets2Read += 1
    return NumOfLets2Read



print(DictOfAln["18"][1] + "\n\n")
print(LocationOfNthNuc(DictOfAln["18"][1], 12))
