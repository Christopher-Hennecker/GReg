from get_GQs import get_GQs
from check_seq import check_seq
from statistics import median

seq_path = input("what is the path for the input file?")
file = open(seq_path, "r")
with open(seq_path, "r") as infh:
    seq = infh.read().strip().replace("\n", "")
max_loop = int(input("What is the maximum number of loop length allowed?"))
max_bulge = int(input("What is the maximum number of bulge allowed?"))
min_temp = int(input("What is the minimum temperature allowed?"))

GQs, scores = get_GQs(max_loop, max_bulge, min_temp)
multi, GRegs = check_seq(seq, GQs, scores, max_loop)

for info in GRegs:
    sequence = seq[info[1]:info[2]]
    each_multi = str(info[0])
    start_index = str(info[1]+1)
    end_index = str(info[2])
    length = str(info[2]-info[1])
    G_number = 0
    for nucleotide in sequence:
        if nucleotide == "G":
            G_number = G_number + 1
    G_content = str(round((G_number / int(length) * 100)))
    Ntot = str(int(sum(info[0])/12))
    Ntand = str(info[4])
    Tm_min = str(round(min(info[3]),1))
    Tm_median = str(round(median(info[3]), 1))
    Tm_max = str(round(max(info[3]), 1))

    position = start_index + "-" + end_index

    output_line = "> position: " + position + " | number of nucleotides: " + length + " | percentage G content: " + G_content + " | total number of quadruplex(es) could form " + Ntot + " | number of tandem repeat(s) " + Ntand + " | maximum estimated melting temperature: " + Tm_max + " | median estimated melting temperature: " + Tm_median + " | minimum estimated melting temperature:" + Tm_min
    
    print(sequence)
    print(each_multi)
    print(output_line)
    print()
