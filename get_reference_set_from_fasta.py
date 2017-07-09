
import sys
from sputnik_reference import find_repeats
#add option for notation: ensemble (1,2,3) o ucsc (chr1, chr2,..)


#chr=int(sys.argv[1]) ##XX
chr=sys.argv[1]
#notation=sys.argv[2])
#if notation not in ["Ensemble","UCSC"]:
#    print "The second argument, correspoding to the chr notation, needs to be Ensemble or UCSC"
#    raise

fi = "chrs_fa/chr"+str(chr)+".fa"

#load fasta file:
# we skip first line
with open(fi, "rb") as myfile:
    sequence = ''.join([line.strip("\n").upper() for line in myfile.readlines()[1:]])


o = find_repeats(sequence,0,chr,min_score=5)
#del sequence

## XX add to REFINE THE SET
## add to find location

f = open("reference_set_"+str(chr)+"_up_to_5rep_in_flank.txt",'w') 

for ms in o:
    repeats_flanking_left = find_repeats(sequence[ms[1]-10:ms[1]],0,chr,min_score=4)
    repeats_flanking_right = find_repeats(sequence[ms[2]:ms[2]+10],0,chr,min_score=4)
    if len(repeats_flanking_left) == 0 and len(repeats_flanking_right) == 0:
        cmd = str(ms[0]) + "\t" +str(ms[1]) + "\t" +str(ms[2])+ "\t"+str(ms[3]) + "\t"+str(ms[4])+ "\t" +str(ms[5]) #+ "\t"+sequence[ms[1]-8:ms[1]] + "\t" + sequence[ms[2]:ms[2]+8]
        f.write(cmd+"\n")

f.close()
