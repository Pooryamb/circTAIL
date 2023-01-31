FoldAdd = "../intermediates/"

FileName = "R_alignment_reversed_pruned.txt"
outFile = "R_alignment_pruned.txt"

file = open(FoldAdd + FileName)
outLines = []
for line in file:
    parts = line.strip().split("\t")
    start_on_seq1 = int(parts[6]) - int(parts[2]) + 1
    start_on_seq2 = int(parts[7]) - int(parts[4]) + 1
    end_on_seq1 = int(parts[6]) - int(parts[1]) + 1
    end_on_seq2 = int(parts[7]) - int(parts[3]) + 1
    qseq = parts[-2][::-1]
    sseq = parts[-1][::-1]
    outLines.append("\t".join([parts[0], str(start_on_seq1), str(end_on_seq1), 
    str(start_on_seq2), str(end_on_seq2), parts[5], parts[6], qseq, sseq]))

file.close()

outFile = open(FoldAdd +outFile,'w')
outFile.write("\n".join(outLines))
outFile.close()
