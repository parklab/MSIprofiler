# to make numpy divisions show decimal values by default:
# https://stackoverflow.com/questions/1799527/numpy-show-decimal-values-in-array-results
from __future__ import division
import csv
import os

import multiprocessing as mp
import numpy as np
import pysam
from scipy import stats

from .utils import loadcsv, phased, unphased


class MicroSatelliteProfiler:
    """
    Class that aids in the detection of microsatellite instability (MSI) from
    sequencing data
    """

    def __init__(self, arguments):
        """
        Constructor for MicroSatelliteProfiler
        :param arguments: parsed argparse.ArgumentParser() object
        """
        self.bed_filename = arguments.bed
        self.fasta_filename = arguments.fasta
        self.flank_size = arguments.flank_size
        self.mapping_quality = arguments.mapping_quality
        self.max_microsatellite_length = arguments.max_MS_length
        self.min_coverage = arguments.min_coverage
        self.min_microsatellite_length = arguments.min_MS_length
        self.mode = arguments.mode
        self.normal_bam = arguments.normal_bam
        self.number_of_processors = arguments.nprocs
        self.output_prefix = arguments.output_prefix
        self.reference_set = arguments.reference_set
        self.repeat_units = set(arguments.repeat_units)
        self.tolerated_mismatches = arguments.tolerated_mismatches
        self.tumor_bam = arguments.tumor_bam

        if self.number_of_processors is None:
            self.number_of_processors = mp.cpu_count()

        if self.number_of_processors == 0:
            raise RuntimeError(
                "The value of the argument `nprocs` needs to be at least 1"
            )

        self.pool = mp.Pool(self.number_of_processors)
        self.chunk_size = int(
            len(self.reference_set) / self.number_of_processors
        )

        if not all(6 >= i > 0 for i in self.repeat_units):
            raise RuntimeError(
                'The repeat units, i.e. ru, '
                'supported are in the 1-6 range. Exiting..'
            )

        if self.mode not in ['both', 'phased', 'unphased']:
            raise RuntimeError(
                'The mode argument needs to be one of the following: '
                'both, phased, unphased.'
            )

        if not os.path.exists(self.tumor_bam):
            raise RuntimeError("Tumor/case bam file does not exist.")

        if not os.path.exists(self.normal_bam):
            raise RuntimeError("Normal bam file does not exist.")

        if (not os.path.exists(self.reference_set) and
                self.mode in ['both', 'unphased']):
            raise RuntimeError("Reference set file does not exist.")

        if not os.path.exists(self.bed_filename) \
                and self.mode in ['both', 'phased']:
            raise RuntimeError(
                "Bed file containing heterozygous SNPs does not exist."
            )

        if not os.path.exists(self.fasta_filename):
            raise RuntimeError(
                "Fasta file correspoding to the reference "
                "genome does not exist. Exiting.."
            )

        self.fasta_file = pysam.FastaFile(self.fasta_filename)

        self.reference_set = [
            x for x in loadcsv(
                self.reference_set,
                self.min_microsatellite_length,
                self.max_microsatellite_length
            )
            ]

        with open(self.bed_filename) as bed:
            reader = csv.reader(bed, delimiter="\t")
            self.sites = list(reader)

    def run_phased(self):
        print "PHASED: Extracting MS repeats from tumor bam file..\n"

        reference_set_ini_end = [x[1] for x in self.reference_set]

        # This list will contain the dictionaries
        # returned by the different processes
        read_lengths_tumor = []

        def log_result(result):
            read_lengths_tumor.append(result)

        if self.number_of_processors == 1:
            read_lengths_tumor = phased(
                self.tumor_bam,
                self.fasta_file,
                self.flank_size,
                self.mapping_quality,
                self.min_coverage,
                self.reference_set,
                reference_set_ini_end,
                self.repeat_units,
                self.sites,
                self.tolerated_mismatches
            )
        else:
            for index in np.arange(0, self.number_of_processors):
                if index != (self.number_of_processors - 1):
                    self.pool.apply_async(
                        phased,
                        args=(
                            self.tumor_bam,
                            self.fasta_file,
                            self.flank_size,
                            self.mapping_quality,
                            self.min_coverage,
                            self.reference_set,
                            reference_set_ini_end,
                            self.repeat_units,
                            self.sites[
                                index *
                                self.chunk_size:(index + 1) *
                                self.chunk_size
                            ],
                            self.tolerated_mismatches
                        ),
                        callback=log_result
                    )
                else:
                    self.pool.apply_async(
                        phased,
                        args=(
                            self.tumor_bam,
                            self.fasta_file,
                            self.flank_size,
                            self.mapping_quality,
                            self.min_coverage,
                            self.reference_set,
                            reference_set_ini_end,
                            self.repeat_units,
                            self.sites[
                                index * self.chunk_size: len(self.sites)
                            ],
                            self.tolerated_mismatches
                        ),
                        callback=log_result
                    )
            # close the pool
            self.pool.close()
            self.pool.join()

        print "PHASED: tumor/case bam file processed correctly..\n"

        print "PHASED: extracting MS repeats from normal bam file..\n"

        read_lengths_normal = []

        def log_result(result):
            read_lengths_normal.append(result)

        if self.number_of_processors == 1:
            read_lengths_normal = phased(
                self.normal_bam,
                self.fasta_file,
                self.flank_size,
                self.mapping_quality,
                self.min_coverage,
                self.reference_set,
                reference_set_ini_end,
                self.repeat_units,
                self.sites,
                self.tolerated_mismatches,
            )
        else:
            for index in np.arange(0, self.number_of_processors):
                if index != (self.number_of_processors - 1):
                   self.pool.apply_async(
                        phased,
                        args=(
                            self.normal_bam,
                            self.fasta_file,
                            self.flank_size,
                            self.mapping_quality,
                            self.min_coverage,
                            self.reference_set,
                            reference_set_ini_end,
                            self.repeat_units,
                            self.sites[
                                index *
                                self.chunk_size:(index + 1) *
                                self.chunk_size
                            ],
                            self.tolerated_mismatches,
                        ),
                        callback=log_result
                    )
                else:
                   self.pool.apply_async(
                        phased,
                        args=(
                            self.normal_bam,
                            self.fasta_file,
                            self.flank_size,
                            self.mapping_quality,
                            self.min_coverage,
                            self.reference_set,
                            reference_set_ini_end,
                            self.repeat_units,
                            self.sites[
                                index * self.chunk_size: len(self.sites)
                            ],
                            self.tolerated_mismatches,
                        ),
                        callback=log_result
                    )
            self.pool.close()
            self.pool.join()

        print "Normal bam file processed correctly..\n"

        with open(self.output_prefix + '_phased.txt', 'w') as f:
            if self.number_of_processors == 1:
                all_normal = read_lengths_normal[0]
                all_tumor = read_lengths_tumor[0]
            else:
                all_normal = read_lengths_normal[0]
                all_tumor = read_lengths_tumor[0]

            if self.number_of_processors > 1:
                for i in range(1, self.number_of_processors):
                    all_normal.update(read_lengths_normal[i])
                    all_tumor.update(read_lengths_tumor[i])

            keys_normal = set(all_normal)
            keys_tumor = set(all_tumor)
            common_keys = keys_tumor.intersection(keys_normal)

            for name in common_keys:
                nor = all_normal[name]
                canc = all_tumor[name]
                if isinstance(nor, int) == False and isinstance(canc,
                                                                int) == False:
                    if len(nor) >= self.min_coverage and len(
                            canc) >= self.min_coverage:
                        pval = stats.ks_2samp(nor, canc)[1]
                        f.write(name + "\t" + ",".join(
                            [str(x) for x in nor]) + "\t" + ",".join(
                            [str(x) for x in canc]) + "\t" + str(pval) + "\n")

        print "Phased microsatellites writen to: {}_phased.txt".format(
            self.output_prefix
        )
        print (
            "Calculation of the phased microsatellites finished successfully.."
        )

    def run_unphased(self):
        print "Extracting MS repeats (UNPHASED) from tumor bam file..\n"

        read_lengths_tumor_unphased = []

        def log_result(result):
            read_lengths_tumor_unphased.append(result)

        if self.number_of_processors == 1:
            read_lengths_tumor_unphased = unphased(
                self.tumor_bam,
                self.fasta_file,
                self.flank_size,
                self.mapping_quality,
                self.min_coverage,
                self.reference_set,
                self.tolerated_mismatches
            )
        else:
            for index in np.arange(0, self.number_of_processors):
                if index != (self.number_of_processors - 1):
                   self.pool.apply_async(
                        unphased,
                        args=(
                            self.tumor_bam,
                            self.fasta_file,
                            self.flank_size,
                            self.mapping_quality,
                            self.min_coverage,
                            self.reference_set[
                                index *
                                self.chunk_size:(index + 1) *
                                self.chunk_size
                            ],
                            self.tolerated_mismatches,
                        ),
                        callback=log_result
                    )
                else:
                   self.pool.apply_async(
                        unphased,
                        args=(
                            self.tumor_bam,
                            self.fasta_file,
                            self.flank_size,
                            self.mapping_quality,
                            self.min_coverage,
                            self.reference_set[
                                index * self.chunk_size:len(self.reference_set)
                            ],
                            self.tolerated_mismatches,
                        ),
                        callback=log_result
                    )
            self.pool.close()
            self.pool.join()

        print "UNPHASED: tumor bam file processed correctly..\n"
        print "Extracting MS repeats (UNPHASED) from normal bam file..\n"

        read_lengths_normal_unphased = []

        def log_result(result):
            read_lengths_normal_unphased.append(result)

        if self.number_of_processors == 1:
            read_lengths_normal_unphased = unphased(
                self.normal_bam,
                self.fasta_file,
                self.flank_size,
                self.mapping_quality,
                self.min_coverage,
                self.reference_set,
                self.tolerated_mismatches
            )
        else:
            for index in np.arange(0, self.number_of_processors):
                if index != (self.number_of_processors - 1):
                   self.pool.apply_async(
                        unphased,
                        args=(
                            self.normal_bam,
                            self.fasta_file,
                            self.flank_size,
                            self.mapping_quality,
                            self.min_coverage,
                            self.reference_set[
                                index *
                                self.chunk_size:(index + 1) *
                                self.chunk_size
                            ],
                            self.tolerated_mismatches,
                        ),
                        callback=log_result
                    )
                else:
                   self.pool.apply_async(
                        unphased,
                        args=(
                            self.normal_bam,
                            self.fasta_file,
                            self.flank_size,
                            self.mapping_quality,
                            self.min_coverage,
                            self.reference_set[
                                index * self.chunk_size:len(self.reference_set)
                            ],
                            self.tolerated_mismatches,
                        ),
                        callback=log_result
                    )
            self.pool.close()
            self.pool.join()

        print "UNPHASED: normal bam file processed correctly..\n"

        with open(self.output_prefix + '_unphased.txt', 'w') as f:
            all_normal = read_lengths_normal_unphased[0]
            all_tumor = read_lengths_tumor_unphased[0]

            if self.number_of_processors > 1:
                for i in range(1, self.number_of_processors):
                    all_normal.update(read_lengths_normal_unphased[i])
                    all_tumor.update(read_lengths_tumor_unphased[i])

            keys_normal = set(all_normal)
            keys_tumor = set(all_tumor)
            common_keys = keys_tumor.intersection(keys_normal)

            for name in common_keys:
                nor = all_normal[name]
                canc = all_tumor[name]
                if isinstance(nor, int) == False and isinstance(canc,
                                                                int) == False:
                    if (len(nor) >= self.min_coverage and
                            len(canc) >= self.min_coverage):
                        pval = stats.ks_2samp(nor, canc)[1]
                        mo = stats.mode(nor)
                        percentage = (nor == mo).sum() / len(nor)
                        confidence = "high" if percentage >= .7 else "low"
                        f.write(name + "\t" + ",".join(
                            [str(x) for x in nor]) + "\t" + ",".join(
                            [str(x) for x in canc]) + "\t" + str(
                            pval) + "\t" + confidence + "\n")
