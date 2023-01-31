def Fasta2Dict(text):
    SeqDict = {}
    parts = text.strip().split(">")[1:]
    for part in parts:
        ID = part.split("\n")[0]
        Seq = "".join(part.strip().split("\n")[1:])
        SeqDict[ID] = Seq
    return SeqDict
