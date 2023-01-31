import glob

TempFile = glob.glob("../input/*_R.fa")[0]
ReadFile = glob.glob("../input/*.fasta")[0]

DictOfFiles = {"TempFile":TempFile, "ReadFile":ReadFile}

def ReversingFasta(LinesOfFile):
    Reversed = []
    N = len(LinesOfFile)
    for i in range(N):
        if i%2==0:
            Reversed.append(LinesOfFile[i].strip() + "\n")
        else:
            Reversed.append(LinesOfFile[i].strip()[::-1].strip()+ "\n")
    return Reversed

for fileName, fileAdd in DictOfFiles.items():
    OnlyName = fileAdd.replace("../input/","")
    fileLines = open(fileAdd).readlines()
    outFile = open("../intermediates/" + ".".join(OnlyName.split(".")[:-1]) + "_reversed" + "." + OnlyName.split(".")[-1],'w')
    outFile.write("".join(ReversingFasta(fileLines)))
    outFile.close()

