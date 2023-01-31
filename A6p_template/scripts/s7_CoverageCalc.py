import os
import numpy as np

LENGTH = {}

sides = ["R","L"]
for side in sides:
    files = os.listdir("../input/")
    file = [x for x in files if x.endswith("_" + side + ".fa")]
    Ref = file[0]
    file = open("../input/" + Ref)
    content = file.readlines()
    file.close()
    LENGTH[side] = len("".join(content[1:]).replace("\n","").strip())


#for each side, the following code will calculate the number of times each nucleotide
#has been appeared in the reads, in other words, it calculates the coverage for every
#position of the reference sequence. The last element of the Counter list is dedicated
#to counting the total number of reads. There will be one empty element between the position
#counters and the total counter

PositionCounter = {}
for side in sides:
    Counter = np.zeros(LENGTH[side] + 2, dtype=int)
    AlignOut = "../intermediates/" + side + "_alignment_pruned.txt" 
    file = open(AlignOut)
#The header of Blat result is 5 lines But sam2psl does not produce header
#    for i in range(5):
#        dumb = file.readline()
    for line in file:
        parts = line.split("\t")
        sstart= int(parts[3]) -1
        send = int(parts[4])
        Counter[sstart:send] +=1
        Counter[-1] += 1
    file.close()
    PositionCounter[side] = Counter

#In the following script, we find the positions that had been observed in at least
#90 percent of the reads
RangesOf90PercentMapped = []

positionsOfStartAndEnd = {}
for side in sides:
    Counts = PositionCounter[side]
    NumOfMappedReads = Counts[-1]
    CountsAbove90Per = (Counts[:-1] > 0.9*NumOfMappedReads) * np.array(list(range(1,len(Counts))))
    Nonzeros = [x for x in CountsAbove90Per if x!=0]
    IndMin=min(Nonzeros)
    IndMax = max(Nonzeros)
    positionsOfStartAndEnd[side] = (IndMin, IndMax)
#os.makedirs("PositionCount",'w')
file=open("../intermediates/Coverage.txt",'w')
file.write("L: \n" + str(PositionCounter["L"])+ "\n\n")
file.write("R: \n" + str(PositionCounter["R"]) + "\n\n")


# As I want to make sure the coding part is always in the UTRs, I make the following adjustments
CodonFileLines = open("../input/CodingLenInTemplate.txt").readlines()
L_limit = int(CodonFileLines[0].strip().split()[1])
R_limit = LENGTH["R"] - int(CodonFileLines[1].strip().split()[1]) + 1

#The written positions start from 1 (not 0)
side = "R"
file.write(side + "\t" + str(min(positionsOfStartAndEnd[side][0], R_limit)) + "\t" + 
    str(positionsOfStartAndEnd[side][1])+"\n")

side = "L"
file.write(side + "\t" + str(positionsOfStartAndEnd[side][0]) + "\t" + 
    str(max(positionsOfStartAndEnd[side][1], L_limit))+"\n")
file.close()