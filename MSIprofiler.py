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
from scipy import stats
from .utils import loadcsv, unphased, phased

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
repeat_units = set(args.rus)

if not all(6 >= i > 0 for i in repeat_units):
    raise RuntimeError(
        'The repeat units, i.e. ru, supported are in the 1-6 range. Exiting..'
    )

if args.mode not in ['both', 'phased', 'unphased']:
    raise RuntimeError(
        'The mode argument needs to be one of the following: '
        'both, phased, unphased. Exiting..'
    )

if not os.path.exists(args.tumor_bam):
    raise RuntimeError("Tumor/case bam file does not exist. Exiting..")

if not os.path.exists(args.normal_bam):
    raise RuntimeError("Normal bam file does not exist. Exiting.")

if (not os.path.exists(args.reference_set) and
            args.mode in ['both', 'unphased']):
    raise RuntimeError("Reference set file does not exist. Exiting..")

if not os.path.exists(args.bed) and args.mode in ['both', 'phased']:
    raise RuntimeError(
        "Bed file containing heterozygous SNPs does not exist. Exiting.."
    )

if not os.path.exists(args.fasta):
    raise RuntimeError(
        "Fasta file correspoding to the reference "
        "genome does not exist. Exiting.."
    )

fastafile = pysam.FastaFile(args.fasta)

if args.mode in ['both', 'phased']:
    with open(args.bed) as bed:
        reader = csv.reader(bed, delimiter="\t")
        sites = list(reader)

# -----------------------------------------------------------------------------
# load reference set
# -----------------------------------------------------------------------------

refsetgen = loadcsv(args.reference_set, args.min_MS_length, args.max_MS_length)
refset = [x for x in refsetgen]
if args.mode in ['both', 'phased']:
    refset_ini_end = [x[1] for x in refset]

# -----------------------------------------------------------------------------
# PHASED
# -----------------------------------------------------------------------------
if args.mode in ['both', 'phased']:
    # ------------------------------------------------------
    print "PHASED: Extracting MS repeats from tumor bam file..\n"
    # ------------------------------------------------------
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
        read_lengths_tumor = phased(
            args.tumor_bam,
            fastafile,
            args.flank_size,
            args.mapping_quality,
            args.min_coverage,
            refset,
            refset_ini_end,
            repeat_units,
            sites,
            args.tolerated_mismatches
        )
    else:
        pool = mp.Pool(args.nprocs)
        chunk_size = int(len(sites) / args.nprocs)
        for index in np.arange(0, args.nprocs):
            if index != (args.nprocs - 1):
                pool.apply_async(
                    phased,
                    args=(
                        args.tumor_bam,
                        fastafile,
                        args.flank_size,
                        args.mapping_quality,
                        args.min_coverage,
                        refset,
                        refset_ini_end,
                        repeat_units,
                        sites[index * chunk_size:(index + 1) * chunk_size],
                        args.tolerated_mismatches
                    ),
                    callback=log_result
                )
            else:
                pool.apply_async(
                    phased,
                    args=(
                        args.tumor_bam,
                        fastafile,
                        args.flank_size,
                        args.mapping_quality,
                        args.min_coverage,
                        refset,
                        refset_ini_end,
                        repeat_units,
                        sites[index * chunk_size: len(sites)],
                        args.tolerated_mismatches
                    ),
                    callback=log_result
                )
        # close the pool
        pool.close()
        pool.join()

    print "PHASED: tumor/case bam file processed correctly..\n"

    print "PHASED: extracting MS repeats from normal bam file..\n"

    read_lengths_normal = []


    def log_result(result):
        read_lengths_normal.append(result)


    if args.nprocs is None:
        args.nprocs = mp.cpu_count()

    if args.nprocs == 1:
        read_lengths_normal = phased(sites, args.normal_bam, 1)
    else:
        pool = mp.Pool(args.nprocs)
        chunk_size = int(len(sites) / args.nprocs)
        for index in np.arange(0, args.nprocs):
            if index != (args.nprocs - 1):
                pool.apply_async(
                    phased,
                    args=(
                        args.normal_bam,
                        fastafile,
                        args.flank_size,
                        args.mapping_quality,
                        args.min_coverage,
                        refset,
                        refset_ini_end,
                        repeat_units,
                        sites[index * chunk_size:(index + 1) * chunk_size],
                        args.tolerated_mismatches,
                    ),
                    callback=log_result
                )
            else:
                pool.apply_async(
                    phased,
                    args=(
                        args.normal_bam,
                        fastafile,
                        args.flank_size,
                        args.mapping_quality,
                        args.min_coverage,
                        refset,
                        refset_ini_end,
                        repeat_units,
                        sites[index * chunk_size: len(sites)],
                        args.tolerated_mismatches,
                    ),
                    callback=log_result
                )
        pool.close()
        pool.join()

    # ------------------------------------------------------
    print "Normal bam file processed correctly..\n"
    # ------------------------------------------------------
    f = open(args.output_prefix + '_phased.txt', 'w')
    if args.nprocs == 1:
        all_normal = read_lengths_normal[0]
        all_tumor = read_lengths_tumor[0]
    else:
        all_normal = read_lengths_normal[0]
        all_tumor = read_lengths_tumor[0]

    if args.nprocs > 1:
        for i in range(1, args.nprocs):
            all_normal.update(read_lengths_normal[i])
            all_tumor.update(read_lengths_tumor[i])

    keys_normal = set(all_normal)
    keys_tumor = set(all_tumor)
    common_keys = keys_tumor.intersection(keys_normal)
    counter = 0

    for name in common_keys:
        nor = all_normal[name]
        canc = all_tumor[name]
        if isinstance(nor, int) == False and isinstance(canc, int) == False:
            if len(nor) >= args.min_coverage and len(
                    canc) >= args.min_coverage:
                pval = stats.ks_2samp(nor, canc)[1]
                f.write(name + "\t" + ",".join(
                    [str(x) for x in nor]) + "\t" + ",".join(
                    [str(x) for x in canc]) + "\t" + str(pval) + "\n")
    f.close()
    print "Phased microsatellites writen to: " + args.output_prefix + '_phased.txt'

    # ------------------------------------------------------
    print "Calculation of the phased microsatellites finished successfully.."
    # ------------------------------------------------------

# -------------------------------------------------------------------------------------------------------------------------------------------
# UNPHASED
# -------------------------------------------------------------------------------------------------------------------------------------------
if args.mode in ['both', 'unphased']:
    # ------------------------------------------------------
    print "Extracting MS repeats (UNPHASED) from tumor bam file..\n"
    # ------------------------------------------------------
    read_lengths_tumor_unphased = []


    def log_result(result):
        read_lengths_tumor_unphased.append(result)


    pool = mp.Pool(args.nprocs)
    chunk_size = int(len(refset) / args.nprocs)

    if args.nprocs == 0:
        print "The value of the argument nprocs needs to be at least 1\n\n"
        raise
    if args.nprocs is None:
        args.nprocs = mp.cpu_count()
    if args.nprocs == 1:
        read_lengths_tumor_unphased = unphased(
            args.tumor_bam,
            fastafile,
            args.flank_size,
            args.mapping_quality,
            args.min_coverage,
            refset,
            args.tolerated_mismatches
        )
    else:
        for index in np.arange(0, args.nprocs):
            if index != (args.nprocs - 1):
                pool.apply_async(
                    unphased,
                    args=(
                        args.tumor_bam,
                        fastafile,
                        args.flank_size,
                        args.mapping_quality,
                        args.min_coverage,
                        refset[index * chunk_size:(index + 1) * chunk_size],
                        args.tolerated_mismatches,
                    ),
                    callback=log_result
                )
            else:
                pool.apply_async(
                    unphased,
                    args=(
                        args.tumor_bam,
                        fastafile,
                        args.flank_size,
                        args.mapping_quality,
                        args.min_coverage,
                        refset[index * chunk_size:len(refset)],
                        args.tolerated_mismatches,
                    ),
                    callback=log_result
                )
        pool.close()
        pool.join()

    # ------------------------------------------------------
    print "UNPHASED: tumor bam file processed correctly..\n"
    # ------------------------------------------------------
    # ------------------------------------------------------
    print "Extracting MS repeats (UNPHASED) from normal bam file..\n"
    # ------------------------------------------------------
    read_lengths_normal_unphased = []


    def log_result(result):
        read_lengths_normal_unphased.append(result)


    pool = mp.Pool(args.nprocs)

    if args.nprocs == 1:
        read_lengths_normal_unphased = unphased(
            args.normal_bam,
            fastafile,
            args.flank_size,
            args.mapping_quality,
            args.min_coverage,
            refset,
            args.tolerated_mismatches
        )
    else:
        for index in np.arange(0, args.nprocs):
            if index != (args.nprocs - 1):
                pool.apply_async(
                    unphased,
                    args=(
                        args.normal_bam,
                        fastafile,
                        args.flank_size,
                        args.mapping_quality,
                        args.min_coverage,
                        refset[index * chunk_size:(index + 1) * chunk_size],
                        args.tolerated_mismatches,
                    ),
                    callback=log_result
                )
            else:
                pool.apply_async(
                    unphased,
                    args=(
                        args.normal_bam,
                        fastafile,
                        args.flank_size,
                        args.mapping_quality,
                        args.min_coverage,
                        refset[index * chunk_size:len(refset)],
                        args.tolerated_mismatches,
                    ),
                    callback=log_result
                )
        pool.close()
        pool.join()

    # ------------------------------------------------------
    print "UNPHASED: normal bam file processed correctly..\n"
    # ------------------------------------------------------
    f = open(args.output_prefix + '_unphased.txt', 'w')
    all_normal = read_lengths_normal_unphased[0]
    all_tumor = read_lengths_tumor_unphased[0]

    if args.nprocs > 1:
        for i in range(1, args.nprocs):
            all_normal.update(read_lengths_normal_unphased[i])
            all_tumor.update(read_lengths_tumor_unphased[i])

    keys_normal = set(all_normal)
    keys_tumor = set(all_tumor)
    common_keys = keys_tumor.intersection(keys_normal)
    counter = 0

    for name in common_keys:
        nor = all_normal[name]
        canc = all_tumor[name]
        if isinstance(nor, int) == False and isinstance(canc, int) == False:
            if len(nor) >= args.min_coverage and len(
                    canc) >= args.min_coverage:
                pval = stats.ks_2samp(nor, canc)[
                    1]  # read_lengths_normal_unphased[i][name], read_lengths_tumor_unphased[i][name])[1]
                mo = stats.mode(nor)
                percentage = (nor == mo).sum() / len(nor)
                confidence = "high" if percentage >= .7 else "low"
                f.write(name + "\t" + ",".join(
                    [str(x) for x in nor]) + "\t" + ",".join(
                    [str(x) for x in canc]) + "\t" + str(
                    pval) + "\t" + confidence + "\n")
    f.close()

# ------------------------------------------------------
print "All calculations finished successfully!\n"
# ------------------------------------------------------
