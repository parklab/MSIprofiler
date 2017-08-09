import argparse
from reps_from_reference  import find_repeats

parser = argparse.ArgumentParser(description='usage: get_reference_set_from_fasta.py --chr 1 --notation Ensemble --min_score 5')
parser.add_argument('--chr', help='Chromosome',required=False,default=22)
# Optional arguments
parser.add_argument('--notation', help='Chromosome notation. The two options are: Ensemble (e.g. 1,2,3,..) or UCSC (e.g. chr1, chr2, chr3,..)',required=False, default="Ensemble",choices=["Ensemble","UCSC"])
parser.add_argument('--min_score', help='Minimum score to detect microsatellites in the flanking regions. Default is 4.',required=False,default=5,type=int)
parser.add_argument('--max_length', help='Maximum microsatellite length. Default is 60.',required=False,default=60,type=int)
args = parser.parse_args()

fi = "chrs_fa/chr"+str(args.chr)+".fa"

#load fasta file:
# we skip first line, as it contains the chr name
with open(fi, "rb") as myfile:
    sequence = ''.join([line.strip("\n").upper() for line in myfile.readlines()[1:]])


o = find_repeats(sequence,0,args.chr,min_score=args.min_score)
#del sequence

f = open("reference_set_"+str(args.chr)+".txt",'w') 

for ms in o:
    repeats_flanking_left = find_repeats(sequence[ms[1]-10:ms[1]],0,args.chr,min_score=args.min_score)
    repeats_flanking_right = find_repeats(sequence[ms[2]:ms[2]+10],0,args.chr,min_score=args.min_score)
    if len(repeats_flanking_left) == 0 and len(repeats_flanking_right) == 0:
        # get the GC contents
        region = sequence[ms[1]-1000:ms[2]+1000] 
        GC = region.count("G") + region.count("C")
        GC = round(GC / float(2000), 1)
        if args.notation == "Ensemble":
            cmd = str(ms[0]) + "\t" +str(ms[1]) + "\t" +str(ms[2])+ "\t"+str(ms[3]) + "\t"+str(ms[4])+ "\t" +str(ms[5]) +"\t"+str(GC)
        else:
            cmd = "chr"+str(ms[0]) + "\t" +str(ms[1]) + "\t" +str(ms[2])+ "\t"+str(ms[3]) + "\t"+str(ms[4])+ "\t" +str(ms[5]) + "\t"+ str(GC)
        f.write(cmd+"\n")

f.close()
