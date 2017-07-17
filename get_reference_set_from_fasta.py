
import sys
#add option for notation: ensemble (1,2,3) o ucsc (chr1, chr2,..)
from utils import find_repeats

chr = sys.argv[1]

fi = "chrs_fa/chr"+str(chr)+".fa"

#load fasta file:
# we skip first line
with open(fi, "rb") as myfile:
    sequence = ''.join([line.strip("\n").upper() for line in myfile.readlines()[1:]])


o = find_repeats(sequence, 0, chr, min_score=5)

## XX add to REFINE THE SET
## add to find location

f = open("reference_set_"+str(chr)+"_up_to_5rep_in_flank.txt",'w')

for ms in o:
    repeats_flanking_left = find_repeats(
        sequence[ms[1]-10:ms[1]],
        0,
        chr,
        min_score=4
    )
    repeats_flanking_right = find_repeats(
        sequence[ms[2]:ms[2]+10],
        0,
        chr,
        min_score=4
    )
    if len(repeats_flanking_left) == 0 and len(repeats_flanking_right) == 0:
        cmd = "{}\t{}\t{}\t{}\t{}\t{}".format(
            ms[0], ms[1], ms[2], ms[3], ms[4], ms[5]
        )
        f.write(cmd + "\n")

f.close()
