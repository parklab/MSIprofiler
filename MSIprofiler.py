rus={1,2,3,4,5} # repeat units considered ## take this from the CONF FILE
match_score = 1 # score for a match
mismatch_score = -6 # score penalty for a mismatch # it detects min length of abs(mismatch_score) + 1
fail_score = -1 # score value to stop search at a given current_pos
min_score = 4 ##minimum score value to pass
#--------------------------------------------------------------------------------------------------------------------------------------
import pysam
import os
import multiprocessing as mp
import argparse
import csv
import numpy as np
import bisect
from bisect_function import binary_search
from scipy import stats
from scipy.stats import mannwhitneyu
#from sputnik import find_repeats
from sputnik_target import find_repeats_target
#from parameters import *
# to make numpy divisions show decimal values by default:
# https://stackoverflow.com/questions/1799527/numpy-show-decimal-values-in-array-results
from __future__ import division
#--------------------------------------------------------------------------------------------------------------------------------------

def find_repeats(seq,flank_size):
    cigar = seq #XX change
    bases=len(seq) # number of bases in the input sequence
    flank_size = flank_size-1

    # save output as a list of lists
    out = []; exclude=set() # use sets: they are much faster with 'in'[]#np.array([],dtype='int')
    
    for ru in rus:
        positions_motif = range(0,ru)
        nb_positions_motif = len(positions_motif)
    # note that the flank is a range, whether python is zero-based
        not_found = True
        base = flank_size 

        while base < bases-flank_size: #and base not in exclude:
            #print base, "BASE"
            if base in exclude:
                base+=1
                continue
            elif not_found:# and base != flank_size:
               # base=test_pos+1
                test_pos=base+ru#+1flank_size
                current_pos=base
            else:
                #base=1+base+mm#test_pos#+1
                #print base
                current_pos=base
                not_found=True
                test_pos=test_pos+ru#+1#2

            pos_in_motif = 0 #update_current_pos(ru)
            score = 0; depth = 0; keep = 0
    
            max_observed_score = 0
            scores = []
            while ( (test_pos ) < (bases-flank_size) )  and  score > fail_score and test_pos not in exclude: #XX the minus one check
            #while score > fail_score and test_pos not in exclude: #XX the minus one check
                #print base, test_pos,"   ", score, max_observed_score, seq[base:test_pos],  exclude, bases-flank_size, "RU",ru
                #print base,test_pos, "   ", score, max_observed_score, seq[base:test_pos], depth, exclude, bases-flank_size, "RU",ru
                match = (seq[current_pos + pos_in_motif] == seq[test_pos])
        
                if match:
                    test_pos+=1
                    pos_in_motif = positions_motif[(pos_in_motif + 1) % nb_positions_motif]
                    score+=match_score
                    scores.append(score)
                    depth = 0
            
                else: # no mismatch: check for N, insertions, deletions and missense
                    score+=mismatch_score
                    scores.append(score)
                    pos_in_motif = positions_motif[(pos_in_motif + 1) % nb_positions_motif]

                    if score > fail_score and depth < 5:
                        depth+=1
                        test_pos +=1
                # keep track of the best observed score
                if score > max_observed_score:
                    max_observed_score = score
                #print "RU current_pos  pos_in_motif  test_pos   bases, score"
                #print ru, current_pos, pos_in_motif, test_pos, bases,score,"\n" #, current_pos + pos_in_motif, bases
            #if test_pos >= bases-flank_size: ## nos metemos en el flaking de la derecha
                    #print base, current_pos,test_pos, bases-flank_size, "test pos  bases-flank_size"
            #        max_observed_score=-100

#            test_pos = test_pos #- depth 
#            if test_pos >= bases-flank_size:
#                not_found = False
#                mm = scores.index(max(scores)) if len(scores)!=0 else 0
#                mm = mm +ru 
#                exclude.update(range(base,base+mm+1)) #np.unique(np.append(exclude,np.arange(base,test_pos)))
#                base =bases+2# just to make it anumber higher than the number of bases


            if max_observed_score >= min_score:# and test_pos <= bases-flank_size: 
                mm = scores.index(max(scores)) 
                mm = mm +ru
                if base+mm < (bases-flank_size): # repeat not overlapping flanking region
                   out.append( [ru, base, base+mm, seq[base:base+mm+1]])#test_pos]] )
                   not_found = False
                exclude.update(range(base,base+mm+1)) #np.unique(np.append(exclude,np.arange(base,test_pos)))
                test_pos = base + mm  ##X
                base = test_pos
            else:
                pass
            # increment base 
            base+=1
        else:
            base+=1
    return out
    
#--------------------------------------------------------------------------------------------------------------------------------------
# https://stackoverflow.com/questions/212358/binary-search-bisection-in-python
def binary_search(a, x, lo=0, hi=None):  # can't use a to specify default for hi
    hi = hi if hi is not None else len(a)  # hi defaults to len(a)   
    pos = bisect_left(a, x, lo, hi)  # find insertion position
    return (pos if pos != hi and a[pos] == x else -1) 
#--------------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='usage: samtools_view.py --bam reads.bam --bed bed_file.bed')
parser.add_argument('--tumor_bam', help='Tumor bam file name',required=True)
parser.add_argument('--normal_bam', help='Normal bam file name',required=True)
parser.add_argument('--bed', help='Input bedfile name',required=True)
parser.add_argument('--fasta', help='Input fasta reference file name',required=True)
parser.add_argument('--reference_set', help='Input reference set of microsatellites',required=True)
parser.add_argument('--output_prefix', help='Prefix for the output files',required=True)
parser.add_argument('--mode', help='Both, phased, unphased',required=True)
#parser.add_argument('--genomic_region', help='Analyse the exome or the whole genome. Accepted values: exome or genome',required=True)
parser.add_argument('--nprocs', help='Number of processes',required=True,type=int)
parser.add_argument('-ru', dest='rus',default=[], help='MS repeats units to be considered', type=int, action='append')

# Optional arguments
parser.add_argument('--min_MS_length', help='Minimum length of microsatellites to be considered. Minimum available is 6; default is 10.',required=False,default=10,type=int)
parser.add_argument('--max_MS_length', help='Maximum length of microsatellites to be considered. Maximum available is 60; default is 60.',required=False,default=60,type=int)
parser.add_argument('--mapping_quality', help='Minimum mapping quality. Default is 40.',required=False,default=40,type=int)
parser.add_argument('--flank_size', help='Minimum length of the flanking regions. Default is 10',required=False,default=10,type=int)
parser.add_argument('--min_coverage', help='Minimum coverage at each MS locus -both in the case and control bams-. Default is 10',required=False,default=10,type=int)
parser.add_argument('--tolerated_mismatches', help='Maximum number of tolerated mismatches in the flanking regions. Default is 0',required=False,default=0,type=int)

args = parser.parse_args()
args.rus = set(args.rus)

if not all(i <= 6 and i >0 for i in args.rus):
    raise Exception('The repeat units, i.e. ru, supported are in the 1-6 range. Exiting..')

if args.mode not in ['both','phased','unphased']:
    raise Exception('The mode argument needs to be one of the following: both, phased, unphased. Exiting..')

#if args.genomic_region not in ['exome','genome']:
#    raise Exception('The argument genomic_region needs to be one of the following: exome or genome. Exiting..')

if os.path.exists(args.tumor_bam) == False:
    raise "Tumor bam file does not exist. Exiting.."

if os.path.exists(args.normal_bam) == False:
    raise "Normal bam file does not exist. Exiting."

if os.path.exists(args.reference_set) == False and args.mode in ['both','unphased']:
    raise "Reference set file does not exist. Exiting.."

if os.path.exists(args.bed) == False and args.mode in ['both','phased']:
    raise "Bed file containing heterozygous SNPs does not exist. Exiting.."

if os.path.exists(args.fasta) == False:
    raise "Fasta file correspoding to the reference genome does not exist. Exiting.."
#--------------------------------------------------------------------------------------------------------------------------------------

fastafile = pysam.FastaFile(args.fasta)
with open(args.bed) as bed:
    reader = csv.reader(bed, delimiter="\t")
    sites = list(reader)

#----------------------------------------------------------------------------------------------------------------------------------------------------
# load reference set
#----------------------------------------------------------------------------------------------------------------------------------------------------
# remove those that are not of the sizes requested
#def loadcsv_exome(filename, criterion1,criterion2):
#    with open(filename, "rb") as csvfile:
#        datareader = csv.reader(csvfile, delimiter="\t")
#        for row in datareader:
#            if int(row[6]) >=criterion1 and int(row[6]) <= criterion2  and row[5] == "exonic" and int(row[4]) in rus:
#                yield row

def loadcsv(filename, criterion1,criterion2):
    with open(filename, "rb") as csvfile:
        datareader = csv.reader(csvfile, delimiter="\t")
        for row in datareader:
            if int(row[5]) >=criterion1 and int(row[5]) <= criterion2: #in ("column header", criterion):
            #if int(row[6]) >=criterion1 and int(row[6]) <= criterion2: #in ("column header", criterion):
                yield row

#if args.genomic_region == "genome":
refsetgen = loadcsv(args.reference_set,args.min_MS_length, args.max_MS_length)
refset = [x for x in refsetgen]
if args.mode in ['both','phased']:
    refset_ini_end = [x[1] for x in refset]

#if args.genomic_region == "exome":
#    refsetgen = loadcsv_exome(args.reference_set,args.min_MS_length, args.max_MS_length)
#    refset = [x for x in refsetgen]
#    if args.mode in ['both','phased']:
#        refset_ini_end = [x[1] for x in refset]
#
#------------------------------------------------------------------------------------------------------------
# Function to extract the MS lengths from the reads
#------------------------------------------------------------------------------------------------------------
def unphased(sites,bam_path):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    dict_out = {}; visited_reads = []
    for site in sites:
        start = int(site[1]);  end = int(site[2])+1;  chr=site[0]; ru=int(site[4])
        reads = [read for read in bamfile.fetch(str(chr), start,start+1,multiple_iterators=True)] 
        reads = [read for read in reads if read.is_proper_pair and read.is_duplicate == False and read.mapping_quality >= args.mapping_quality] 
        if  len(reads) > args.min_coverage:
            for read in reads:
                start_read = read.reference_start; end_read = read.reference_end
                #read_sequence = read.seq XXX
                # to remove soft-clipping, we can use
                read_sequence = read.query_alignment_sequence ##XX the docs are here: http://pysam.readthedocs.io/en/latest/api.html
                reps = find_repeats_target(read_sequence,args.flank_size,ru) #XX
                if len(reps) > 0: 
                    aligned_pos = read.get_reference_positions(full_length=True)
                    try:
                        idx = aligned_pos.index(start) 
                    except:
                        continue
                    for microsatellite in reps:
                        ru = microsatellite[0]; rs = microsatellite[1]; re = microsatellite[2]
                        if start != start_read + rs + 1:#XXX mal si hay ins/del antes del inicio del STR
                            continue 
                        difference = re - rs + 1

                        # get flinking sequence from reference
                        flank_left_ref = fastafile.fetch("chr"+chr,start_read+rs-args.flank_size, start_read+rs).upper()
                        flank_right_ref = fastafile.fetch("chr"+chr,int(site[2])-1,int(site[2])-1 +args.flank_size).upper() 
                        # get flinking sequence from the reads
                        posfl = (start_read+rs-args.flank_size)
                        if posfl >= start_read:
                            flank_left = read_sequence[rs-args.flank_size:rs];  mismatches_left = sum(a!=b for a, b in zip(flank_left,flank_left_ref))
                        else:
                            flank_left = ""; mismatches_left=10000

                        posflr = start_read+re+args.flank_size
                        if posflr <= end_read:
                            #flank_right = read_sequence[re:re+args.flank_size]; mismatches_right = sum(a!=b for a, b in zip(flank_right,flank_right_ref))
                            flank_right = read_sequence[re+1:re+args.flank_size+1]; mismatches_right = sum(a!=b for a, b in zip(flank_right,flank_right_ref))
                        else:
                            flank_right = ""; mismatches_right=10000

                        mismatches = mismatches_left + mismatches_right 

                        if mismatches <= args.tolerated_mismatches: 
                            key_now = site[0] + "\t"+site[1]+"\t"+site[2]+"\t"+site[3]+"\t"+site[4]+"\t"+site[5]+"\t"+site[6]
                            if dict_out.has_key(key_now):
                                dict_out[key_now] = np.append(dict_out[key_now], difference)
                            else:
                                dict_out[key_now] = difference
    bamfile.close()
    return dict_out

#------------------------------------------------------------------------------------------------------------
def phased(sites,bam_path,index):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    dict_out = {}; visited_reads = []
    for site in sites:
        start = int(site[1]);  end = int(site[2]);  chr=str(site[0]); base1 = site[3]; base2=site[4]; bases = [site[3], site[4] ]
        reads = [read for read in bamfile.fetch(chr, start,end ,multiple_iterators=True)] ## keep this as is, do not put conditions inside here
        reads = [read for read in reads if read.is_proper_pair and read.is_duplicate == False and read.mapping_quality >= args.mapping_quality] 
        if  len(reads) > args.min_coverage:
            for read in reads:
                #read_sequence = read.seq
                read_sequence = read.query_alignment_sequence
                reps = find_repeats(read_sequence,args.flank_size)
                if len(reps) > 0: 
                    # get the SNP in this read
                    start_read = read.reference_start; end_read = read.reference_end
                    aligned_pos = read.get_reference_positions(full_length=True)
                    try:
                        idx = aligned_pos.index(start) 
                    except:
                        continue
                    snp_read = read_sequence[idx]
                    if snp_read not in bases:
                        continue
                    for microsatellite in reps:
                        ru = microsatellite[0]; rs = microsatellite[1]; re = microsatellite[2]; difference = re - rs +1
                    
                        # use the reference set here to get the position on the right
                        ini = start_read + rs  #
                        #ini = aligned_pos[rs]  #using aligned_pos in case there are mismatches in the read before the MS. +1 to convert 1-based
                        #try: # check if this repeat is in the reference set
                            #idx2 = refset_ini_end.index(str(ini+1))
                            #idx2 = bisect.bisect(refset_ini_end,str(ini+1)) -1 # be careful, bisect gives the indices in 1-based format
                        idx2=binary_search(refset_ini_end,str(ini+1))
                        if idx2 == -1:
                        #except:
                            continue
                       
                        refset_now = refset[idx2]
                        diff_ref = int(refset_now[2])- int(refset_now[1])  + 1
                        #flank_right_ref = fastafile.fetch("chr"+str(site[0]),start_read+rs+diff_ref-1, start_read+rs+diff_ref-1+args.flank_size).upper()
                        #flank_left_ref = fastafile.fetch("chr"+str(site[0]),start_read+rs-args.flank_size, start_read+rs).upper()
                        # considering alignment mismatches
                        flank_right_ref = fastafile.fetch("chr"+str(site[0]), ini +diff_ref, ini +diff_ref+args.flank_size).upper()
                        flank_left_ref = fastafile.fetch("chr"+str(site[0]),ini-args.flank_size, ini).upper()

                        posfl = (start_read+rs-args.flank_size)
                        if posfl >= start_read:
                            flank_left = read_sequence[rs-args.flank_size:rs]; mismatches_left = sum(a!=b for a, b in zip(flank_left,flank_left_ref))
                        else:
                            flank_left = ""; mismatches_left = 10000
                        posflr = start_read+re+args.flank_size
                        if posflr <= end_read:
                            flank_right = read_sequence[re+1:re+1+args.flank_size]; mismatches_right = sum(a!=b for a, b in zip(flank_right,flank_right_ref))
                        else:
                            flank_right = ""; mismatches_right = 10000
                        mismatches = mismatches_left + mismatches_right 
                        #refse = fastafile.fetch("chr"+str(site[0]),start_read+rs-args.flank_size,  start_read+rs+diff_ref +args.flank_size).upper()
                        #readse = read_sequence[rs-args.flank_size: re+1+args.flank_size]
                        #if difference != diff_ref:
                        #print microsatellite,difference, site,snp_read,read_sequence[idx-2:idx+2],read_sequence[idx-10:idx+10], diff_ref,flank_right, flank_right_ref, "     ", flank_left, flank_left_ref,"  ",ini,start_read+rs,"\n\n"
                        #print read,"\n\n"
                            #print microsatellite,flank_right, flank_right_ref,flank_right_ref, "     ", flank_left, flank_left_ref,flank_left_ref,"  ",ini,start_read+rs,"\n\n"
                        if mismatches <= args.tolerated_mismatches: 
                            start_STR = ini #start_read +rs;
                            key_now = site[0] + "\t"+str(start_STR)+"\t"+snp_read+"\t"+str(site[1])
                            if dict_out.has_key(key_now):
                                dict_out[key_now] = np.append(dict_out[key_now], difference)
                            else:
                                dict_out[key_now] = difference
    bamfile.close()
    return dict_out

#--------------------------------------------------------------------------------------------------------------------------------------------------
# PHASED
#--------------------------------------------------------------------------------------------------------------------------------------------------
if args.mode in ['both','phased']:
    #------------------------------------------------------
    print "PHASED: Extracting MS repeats from tumor bam file..\n"
    #------------------------------------------------------
    # this list will contain the dictionaries returned by the different processes
    read_lengths_tumor = []

    def log_result(result):
        read_lengths_tumor.append(result)
    
    if args.nprocs == None:
        args.nprocs = mp.cpu_count() 
    if args.nprocs == 0:
        print "The value of the argument nprocs needs to be at least 1\n\n"
        raise
    if args.nprocs == 1:
        read_lengths_tumor = phased(sites,args.tumor_bam,1)
    else:
        pool = mp.Pool(args.nprocs)
        chunk_size= len(sites)/args.nprocs
        for index in np.arange(0,args.nprocs):
            #print index, index*chunk_size, (index+1)*chunk_size
            if index != (args.nprocs-1):
                pool.apply_async(phased, args = (sites[index*chunk_size:(index+1)*chunk_size], args.tumor_bam,index,), callback = log_result)
            else:
                pool.apply_async(phased, args = (sites[index*chunk_size: len(sites)], args.tumor_bam,index,), callback = log_result)
        # close the pool
        pool.close()
        pool.join()

    
    #------------------------------------------------------
    print "PHASED: tumor bam file processed correctly..\n"
    #------------------------------------------------------
    #------------------------------------------------------------------------------------------------------------
    #------------------------------------------------------
    print "PHASED: extracting MS repeats from normal bam file..\n"
    #------------------------------------------------------
    read_lengths_normal = []
    def log_result(result):
        read_lengths_normal.append(result)
    
    if args.nprocs == None:
        args.nprocs = mp.cpu_count() 
    
    if args.nprocs == 1:
        read_lengths_normal = phased(sites,args.normal_bam,1)
    else:
        pool = mp.Pool(args.nprocs)
        chunk_size= len(sites)/args.nprocs
        for index in np.arange(0,args.nprocs):
            if index != (args.nprocs-1):
                pool.apply_async(phased, args = (sites[index*chunk_size:(index+1)*chunk_size], args.normal_bam,index,), callback = log_result)
            else:
                pool.apply_async(phased, args = (sites[index*chunk_size:len(sites)], args.normal_bam,index,), callback = log_result)
        pool.close()
        pool.join()
    
    #------------------------------------------------------
    print "Normal bam file processed correctly..\n"
    #------------------------------------------------------
    f = open(args.output_prefix+'_phased.txt', 'w')
    if args.nprocs ==1:
        all_normal = read_lengths_normal[0]
        all_tumor = read_lengths_tumor[0]
    else:
        all_normal = read_lengths_normal[0]
        all_tumor = read_lengths_tumor[0]

    if args.nprocs >1:
        for i in range(1,args.nprocs):
            all_normal.update(read_lengths_normal[i])
            all_tumor.update(read_lengths_tumor[i])

    keys_normal =  set(all_normal);  keys_tumor =  set(all_tumor)
    common_keys= keys_tumor.intersection(keys_normal); counter = 0
        
    for name in common_keys:
        nor = all_normal[name]
        canc = all_tumor[name]
        if isinstance(nor,int) == False and isinstance(canc,int) == False:
            if len(nor) >= args.min_coverage and len(canc) >= args.min_coverage:
                pval = stats.ks_2samp(nor, canc)[1] ##read_lengths_normal[i][name], read_lengths_tumor[i][name])[1]
                f.write(name+"\t"+ ",".join([str(x) for x in nor])  +"\t"+ ",".join([str(x) for x in canc ]) +"\t"+str(pval)+"\n" )
    f.close()
    print "Phased microsatellites writen to: "+args.output_prefix+'_phased.txt'
    
    #------------------------------------------------------
    print "Calculation of the phased microsatellites successfully finished.."
    #------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------------
# UNPHASED
#-------------------------------------------------------------------------------------------------------------------------------------------
if args.mode in ['both','unphased']:
    
    #------------------------------------------------------
    print "Extracting MS repeats (UNPHASED) from tumor bam file..\n"
    #------------------------------------------------------
    read_lengths_tumor_unphased = []
    def log_result(result):
        read_lengths_tumor_unphased.append(result)
    pool = mp.Pool(args.nprocs)
    chunk_size= len(refset)/args.nprocs
    print chunk_size
    
    if args.nprocs == 0:
        print "The value of the argument nprocs needs to be at least 1\n\n"
        raise
    if args.nprocs == None:
        args.nprocs = mp.cpu_count() 
    if args.nprocs == 1:
        read_lengths_tumor_unphased = phased(sites,args.tumor_bam)
    else:
        for index in np.arange(0,args.nprocs):
            if index != (args.nprocs-1):
                pool.apply_async(unphased, args = (refset[index*chunk_size:(index+1)*chunk_size], args.tumor_bam,), callback = log_result)
            else:
                pool.apply_async(unphased, args = (refset[index*chunk_size:len(refset)], args.tumor_bam,), callback = log_result)
            #pool.apply_async(unphased, args = (refset[index*chunk_size:(index+1)*chunk_size], args.tumor_bam,), callback = log_result)
        pool.close()
        pool.join()
    
    #------------------------------------------------------
    print "UNPHASED: tumor bam file processed correctly..\n"
    #------------------------------------------------------
    #------------------------------------------------------
    print "Extracting MS repeats (UNPHASED) from normal bam file..\n"
    #------------------------------------------------------
    read_lengths_normal_unphased = []
    def log_result(result):
        read_lengths_normal_unphased.append(result)
    pool = mp.Pool(args.nprocs)
    
    if args.nprocs == 1:
        read_lengths_normal_unphased = phased(sites,args.normal_bam)
    else:
        for index in np.arange(0,args.nprocs):
            if index != (args.nprocs-1):
                pool.apply_async(unphased, args = (refset[index*chunk_size:(index+1)*chunk_size], args.normal_bam,), callback = log_result)
            else:
                pool.apply_async(unphased, args = (refset[index*chunk_size:len(refset)], args.normal_bam,), callback = log_result)
            #pool.apply_async(unphased, args = (refset[index*chunk_size:(index+1)*chunk_size], args.normal_bam,), callback = log_result)
        pool.close()
        pool.join()

    #------------------------------------------------------
    print "UNPHASED: normal bam file processed correctly..\n"
    #------------------------------------------------------
    f = open(args.output_prefix+'_unphased.txt', 'w')
    all_normal = read_lengths_normal_unphased[0]
    all_tumor = read_lengths_tumor_unphased[0]

    if args.nprocs >1:
        for i in range(1,args.nprocs):
            all_normal.update(read_lengths_normal_unphased[i])
            all_tumor.update(read_lengths_tumor_unphased[i])

    keys_normal =  set(all_normal)
    keys_tumor =  set(all_tumor)
    common_keys= keys_tumor.intersection(keys_normal)
    counter = 0
    #for name in common_keys:
    #    nor = all_normal[name]
    #    canc = all_tumor[name]
    #    if isinstance(nor,int) == False and isinstance(canc,int) == False:
    #        if len(nor) >= args.min_coverage and len(canc) >= args.min_coverage:
    #            counter+=1
    #            pval = stats.ks_2samp(nor, canc)[1] ##read_lengths_normal[i][name], read_lengths_tumor[i][name])[1]
    #            f.write(name+"\t"+ ",".join([str(x) for x in nor])  +"\t"+ ",".join([str(x) for x in canc ]) +"\t"+str(pval)+"\n" )
    #f.close()
# XXXX
    
    #for i in range(0,args.nprocs):
    #    keys_normal =  set(read_lengths_normal_unphased[i])
    #    keys_tumor =  set(read_lengths_tumor_unphased[i])
    #    common_keys= keys_tumor.intersection(keys_normal)
    #    print len(common_keys), len(keys_normal), len(keys_tumor)
        
    for name in common_keys:
        nor = all_normal[name]
        canc = all_tumor[name]
        if isinstance(nor,int) == False and isinstance(canc,int) == False:
            pval = stats.ks_2samp(nor,canc)[1] #read_lengths_normal_unphased[i][name], read_lengths_tumor_unphased[i][name])[1]
            mo = stats.mode(nor)
            percentage = (nor == mo).sum() / len(nor)
            confidence = "high" if percentage >=.7 else "low"
            f.write(name+"\t"+ ",".join([str(x) for x in nor])  +"\t"+ ",".join([str(x) for x in canc ]) +"\t"+str(pval)+"\t"+confidence+"\n" )
    f.close()

#------------------------------------------------------
print "All calculations finished successfully!\n"
#------------------------------------------------------
