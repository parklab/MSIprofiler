# Isidro Cortes-Ciriano
# Harvard Medical School
# isidrolauscher@gmail.com

# to make numpy divisions show decimal values by default:
# https://stackoverflow.com/questions/1799527/numpy-show-decimal-values-in-array-results
from __future__ import division

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
#--------------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='MSIprofiler serves to detect microsatellite instability from sequencing data. Type MSIprofiler.py --help for further instructions.')
parser.add_argument('--tumor_bam', help='Tumor bam file name',required=True)
parser.add_argument('--normal_bam', help='Normal bam file name',required=True)
parser.add_argument('--bed', help='Input bedfile name',required=True)
parser.add_argument('--fasta', help='Input fasta reference file name',required=True)
parser.add_argument('--reference_set', help='Input reference set of microsatellites',required=True)
parser.add_argument('--output_prefix', help='Prefix for the output files',required=True)
parser.add_argument('--mode', help='Both, phased, unphased',required=True)
parser.add_argument('--nprocs', help='Number of processes',required=True,type=int)
parser.add_argument('-ru', dest='rus',default=[], help='MS repeats units to be considered', type=int, action='append')
parser = argparse.ArgumentParser(
    description=(
        'MSIprofiler serves to detect microsatellite instability '
        'from sequencing data. Type MSIprofiler.py --help for '
        'further instructions.'
    )
)
parser.add_argument(
    '--tumor_bam',
    help='Tumor bam file name',
    required=True
)
parser.add_argument(
    '--normal_bam',
    help='Normal bam file name',
    required=True
)
parser.add_argument(
    '--bed',
    help='Input bedfile name',
    required=True
)
parser.add_argument(
    '--fasta',
    help='Input fasta reference file name',
    required=True
)
parser.add_argument(
    '--reference_set',
    help='Input reference set of microsatellites',
    required=True
)
parser.add_argument(
    '--output_prefix',
    help='Prefix for the output files',
    required=True
)
parser.add_argument(
    '--mode',
    help='Both, phased, unphased',
    required=True
)
parser.add_argument(
    '--nprocs',
    help='Number of processes',
    required=True,
    type=int
)
parser.add_argument(
    '-ru',
    dest='rus',
    default=[],
    help='MS repeats units to be considered',
    type=int,
    action='append'
)

# Optional arguments
parser.add_argument(
    '--min_MS_length',
    help=(
        'Minimum length of microsatellites to be considered. '
        'Minimum available is 6; default is 10.'
    ),
    required=False,
    default=10,
    type=int
)
parser.add_argument(
    '--max_MS_length',
    help=(
        'Maximum length of microsatellites to be considered. '
        'Maximum available is 60; default is 60.'
    ),
    required=False,
    default=60,
    type=int
)
parser.add_argument(
    '--mapping_quality',
    help='Minimum mapping quality. Default is 40.',
    required=False,
    default=40,
    type=int
)
parser.add_argument(
    '--flank_size',
    help='Minimum length of the flanking regions. Default is 10',
    required=False,
    default=10,
    type=int
)
parser.add_argument(
    '--min_coverage',
    help=(
        'Minimum coverage at each MS locus -both in the case and control bams. '
        'Default is 10'
    ),
    required=False,
    default=10,
    type=int
)
parser.add_argument(
    '--tolerated_mismatches',
    help='Maximum number of tolerated mismatches in the flanking regions. Default is 0',
    required=False,
    default=0,
    type=int
)

args = parser.parse_args()
args.rus = set(args.rus)
rus=args.rus

if not all(i <= 6 and i >0 for i in args.rus):
    raise Exception('The repeat units, i.e. ru, supported are in the 1-6 range. Exiting..')

if args.mode not in ['both','phased','unphased']:
    raise Exception('The mode argument needs to be one of the following: both, phased, unphased. Exiting..')

if os.path.exists(args.tumor_bam) == False:
    raise "Tumor/case bam file does not exist. Exiting.."

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

if args.mode in ['both','phased']:
    with open(args.bed) as bed:
        reader = csv.reader(bed, delimiter="\t")
        sites = list(reader)

#----------------------------------------------------------------------------------------------------------------------------------------------------
# load reference set
#----------------------------------------------------------------------------------------------------------------------------------------------------


refsetgen = loadcsv(args.reference_set,args.min_MS_length, args.max_MS_length)
refset = [x for x in refsetgen]
if args.mode in ['both','phased']:
    refset_ini_end = [x[1] for x in refset]

#------------------------------------------------------------------------------------------------------------
# Function to extract the MS lengths from the reads
#------------------------------------------------------------------------------------------------------------
def unphased(sites,bam_path):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    dict_out = {}; visited_reads = []
    for site in sites:
        start = int(site[1]); end = int(site[2])+1; chr=site[0]; ru=int(site[4])
        reads = [read for read in bamfile.fetch(str(chr), start,start+1,multiple_iterators=True)] 
        reads = [read for read in reads if read.is_proper_pair and read.is_duplicate == False and read.mapping_quality >= args.mapping_quality] 
        if  len(reads) > args.min_coverage:
            for read in reads:
                start_read = read.reference_start; end_read = read.reference_end
                read_sequence = read.seq
                # to remove soft-clipping, we can use
                #read_sequence = read.query_alignment_sequence ## the docs are here: http://pysam.readthedocs.io/en/latest/api.html
                reps = find_repeats_target(read_sequence,args.flank_size,ru)
                if len(reps) > 0: 
                    aligned_pos = read.get_reference_positions(full_length=True)
                    try:
                        idx = aligned_pos.index(start) 
                    except:
                        continue
                    for microsatellite in reps:
                        ru = microsatellite[0]; rs = microsatellite[1]; re = microsatellite[2]
                        if start != start_read + rs + 1:# do not consider if there are ins/del upstream of the repeat
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
                            flank_right = read_sequence[re:re+args.flank_size]; mismatches_right = sum(a!=b for a, b in zip(flank_right,flank_right_ref))
                        else:
                            flank_right = ""; mismatches_right=10000
                        mismatches = mismatches_left + mismatches_right 
                        if mismatches <= args.tolerated_mismatches: 
                            key_now = site[0] + "\t"+site[1]+"\t"+site[2]+"\t"+site[3]+"\t"+site[4]+"\t"+site[5]#+"\t"+site[6]
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
                read_sequence = read.seq
                #read_sequence = read.query_alignment_sequence
                reps = find_repeats(read_sequence,args.flank_size)
                if len(reps) > 0: 
                    # get the SNP allele in this read
                    start_read = read.reference_start; end_read = read.reference_end
                    aligned_pos = read.get_reference_positions(full_length=True) #True) reports none for soft-clipped positions
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
                        idx2=binary_search(refset_ini_end,str(ini+1))
                        #print "read, start_read, aligned_pos, idx, snp_read, microsatellite"
                        #print read_sequence,start_read, aligned_pos, idx, snp_read, microsatellite,read_sequence[idx-3:idx+3],"\n\n\n"
                        if idx2 == -1:
                            continue
                        refset_now = refset[idx2]
                        diff_ref = int(refset_now[2])- int(refset_now[1])  + 1
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
                        #print microsatellite,difference, site,snp_read,read_sequence[idx-2:idx+2],read_sequence[idx-10:idx+10], diff_ref,flank_right, flank_right_ref, "     ", flank_left, flank_left_ref,"  ",ini,start_read+rs,"\n\n"
                        #print read,"\n\n"
                            #print microsatellite,flank_right, flank_right_ref,flank_right_ref, "     ", flank_left, flank_left_ref,flank_left_ref,"  ",ini,start_read+rs,"\n\n"
                        if mismatches <= args.tolerated_mismatches: 
                            key_now = site[0] + "\t"+str(ini)+"\t"+snp_read+"\t"+str(site[1])
                            if dict_out.has_key(key_now):
                                dict_out[key_now] = np.append(dict_out[key_now], difference)
                            else:
                                dict_out[key_now] = difference
    bamfile.close()
    #print dict_out
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
        chunk_size= int(len(sites)/args.nprocs)
        for index in np.arange(0,args.nprocs):
            if index != (args.nprocs-1):
                pool.apply_async(phased, args = (sites[index*chunk_size:(index+1)*chunk_size], args.tumor_bam,index,), callback = log_result)
            else:
                pool.apply_async(phased, args = (sites[index*chunk_size: len(sites)], args.tumor_bam,index,), callback = log_result)
        # close the pool
        pool.close()
        pool.join()

    
    #------------------------------------------------------
    print "PHASED: tumor/case bam file processed correctly..\n"
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
        chunk_size= int(len(sites)/args.nprocs)
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
                pval = stats.ks_2samp(nor, canc)[1] 
                f.write(name+"\t"+ ",".join([str(x) for x in nor])  +"\t"+ ",".join([str(x) for x in canc ]) +"\t"+str(pval)+"\n" )
    f.close()
    print "Phased microsatellites writen to: "+args.output_prefix+'_phased.txt'
    
    #------------------------------------------------------
    print "Calculation of the phased microsatellites finished successfully.."
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
    chunk_size= int( len(refset)/args.nprocs )
    
    if args.nprocs == 0:
        print "The value of the argument nprocs needs to be at least 1\n\n"
        raise
    if args.nprocs == None:
        args.nprocs = mp.cpu_count() 
    if args.nprocs == 1:
        read_lengths_tumor_unphased = unphased(refset,args.tumor_bam)
    else:
        for index in np.arange(0,args.nprocs):
            if index != (args.nprocs-1):
                pool.apply_async(unphased, args = (refset[index*chunk_size:(index+1)*chunk_size], args.tumor_bam,), callback = log_result)
            else:
                pool.apply_async(unphased, args = (refset[index*chunk_size:len(refset)], args.tumor_bam,), callback = log_result)
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
        read_lengths_normal_unphased = unphased(refset,args.normal_bam)
    else:
        for index in np.arange(0,args.nprocs):
            if index != (args.nprocs-1):
                pool.apply_async(unphased, args = (refset[index*chunk_size:(index+1)*chunk_size], args.normal_bam,), callback = log_result)
            else:
                pool.apply_async(unphased, args = (refset[index*chunk_size:len(refset)], args.normal_bam,), callback = log_result)
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
        
    for name in common_keys:
        nor = all_normal[name]
        canc = all_tumor[name]
        if isinstance(nor,int) == False and isinstance(canc,int) == False:
            if len(nor) >= args.min_coverage and len(canc) >= args.min_coverage:
                pval = stats.ks_2samp(nor,canc)[1] #read_lengths_normal_unphased[i][name], read_lengths_tumor_unphased[i][name])[1]
                mo = stats.mode(nor)
                percentage = (nor == mo).sum() / len(nor)
                confidence = "high" if percentage >=.7 else "low"
                f.write(name+"\t"+ ",".join([str(x) for x in nor])  +"\t"+ ",".join([str(x) for x in canc ]) +"\t"+str(pval)+"\t"+confidence+"\n" )
    f.close()

#------------------------------------------------------
print "All calculations finished successfully!\n"
#------------------------------------------------------
