# to make numpy divisions show decimal values by default:
# https://stackoverflow.com/questions/1799527/numpy-show-decimal-values-in-array-results
from __future__ import division
import csv
from os import path
import multiprocessing as mp
import numpy as np
from scipy import stats

import utils


class MicroSatelliteProfiler:
    """
    Class that aids in the detection of microsatellite instability (MSI) from
    sequencing data
    """
    BED_FILE_ERROR_MESSAGE = (
        "Bed file containing heterozygous SNPs does not exist."
    )
    VALID_REPEAT_UNITS = [1, 2, 3, 4, 5, 6]

    CHROMOSOMES = [str(i) for i in range(1, 23)]
    CHROMOSOMES.extend(["X", "Y"])
    CHROMOSOMES_ERROR_MESSAGE = (
        "Valid chromosomes are: {}".format(CHROMOSOMES)
    )
    FASTA_DIRECTORY_ERROR_MESSAGE = (
        "Fasta directory correspoding to the reference genome does not exist."
    )
    NORMAL = "normal"
    NORMAL_BAM_ERROR_MESSAGE = "Normal bam file does not exist."
    NUMBER_OF_PROCESSORS_ERROR_MESSAGE = (
        "The value of the argument `nprocs` needs to be at least 1"
    )
    PHASED = "phased"
    REFERENCE_SET_ERROR_MESSAGE = "Reference set file does not exist."
    REPEAT_UNITS_ERROR_MESSAGE = (
        "Valid repeat_units are {} or any "
        "combination thereof".format(VALID_REPEAT_UNITS)
    )
    TUMOR = "tumor"
    TUMOR_BAM_ERROR_MESSAGE = "Tumor/case bam file does not exist."
    UNPHASED = "unphased"
    VALID_MODES = [PHASED, UNPHASED]
    VALID_MODES_ERROR_MESSAGE = (
        'The mode argument needs to be one of the following: {}.'.format(
            VALID_MODES
        )
    )

    def __init__(self, arguments):
        """
        Constructor for MicroSatelliteProfiler
        :param arguments: parsed argparse.ArgumentParser() object
        """
        self.bed_filename = arguments.bed
        self.chromosomes = np.sort(arguments.chromosomes)
        self.fasta_dict = {}
        self.fasta_directory = arguments.fasta
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
        self.repeat_units = set(arguments.rus)
        self.tolerated_mismatches = arguments.tolerated_mismatches
        self.tumor_bam = arguments.tumor_bam
        self.reference_sets = None
        self.sites = None
        self.reference_set_dict = {}
        self.reference_set_ini_end_dict = {}
        self.read_lengths_normal = []
        self.read_lengths_tumor = []
        self.read_lengths_normal_unphased = []
        self.read_lengths_tumor_unphased = []

        self._validate_arguments()

        self._set_fasta_dict()
        self._populate_reference_sets()

        if self.is_phased:
            with open(self.bed_filename) as bed:
                reader = csv.reader(bed, delimiter="\t")
                self.sites = list(reader)
            self.chunk_size = int(len(self.sites) / self.number_of_processors)

        if self.is_unphased:
            self.chunk_size = int(
                len(self.reference_sets) / self.number_of_processors
            )

    @property
    def is_phased(self):
        return True if self.mode == self.PHASED else False

    @property
    def is_unphased(self):
        return True if self.mode == self.UNPHASED else False

    def _check_bams(self):
        if not path.exists(self.tumor_bam):
            raise RuntimeError(self.TUMOR_BAM_ERROR_MESSAGE)

        if not path.exists(self.normal_bam):
            raise RuntimeError(self.NORMAL_BAM_ERROR_MESSAGE)

    def _check_bed_filename(self):
        if not path.exists(self.bed_filename) and self.mode == self.PHASED:
            raise RuntimeError(self.BED_FILE_ERROR_MESSAGE)

    def _check_chromosomes(self):
        for chromosome in self.chromosomes:
            if chromosome not in self.CHROMOSOMES:
                raise RuntimeError(self.CHROMOSOMES_ERROR_MESSAGE)

    def _check_fasta_filename(self):
        if not path.exists(self.fasta_directory):
            raise RuntimeError(self.FASTA_DIRECTORY_ERROR_MESSAGE)

    def _check_mode(self):
        if self.mode not in self.VALID_MODES:
            raise RuntimeError(self.VALID_MODES_ERROR_MESSAGE)

    def _check_processors(self):
        if self.number_of_processors is None:
            self.number_of_processors = mp.cpu_count()

        if self.number_of_processors == 0:
            raise RuntimeError(self.NUMBER_OF_PROCESSORS_ERROR_MESSAGE)

    def _check_reference_set(self):
        if not path.exists(self.reference_set) and self.mode == self.UNPHASED:
            raise RuntimeError(self.REFERENCE_SET_ERROR_MESSAGE)

    def _check_repeat_units(self):
        for repeat_unit in self.repeat_units:
            if repeat_unit not in self.VALID_REPEAT_UNITS:
                raise RuntimeError(self.REPEAT_UNITS_ERROR_MESSAGE)

    def _conclude_run(self):
        if self.number_of_processors > 1:
            all_normal = {}
            all_tumor = {}
            for i in range(0, self.number_of_processors):
                if self.is_phased:
                    all_normal.update(self.read_lengths_normal[i])
                    all_tumor.update(self.read_lengths_tumor[i])
                else:
                    all_normal.update(self.read_lengths_normal_unphased[i])
                    all_tumor.update(self.read_lengths_tumor_unphased[i])
        else:
            if self.is_phased:
                all_normal = self.read_lengths_normal[0]
                all_tumor = self.read_lengths_tumor[0]
            else:
                all_normal = self.read_lengths_normal_unphased[0]
                all_tumor = self.read_lengths_tumor_unphased[0]

        keys_normal = set(all_normal)
        keys_tumor = set(all_tumor)
        common_keys = keys_tumor.intersection(keys_normal)

        with open('{}_{}.txt'.format(self.output_prefix, self.mode), 'w') as f:
            if self.is_phased:
                self._write_phased_output(
                    f,
                    common_keys,
                    all_normal,
                    all_tumor
                )
            else:
                self._write_unphased_output(
                    f,
                    common_keys,
                    all_normal,
                    all_tumor
                )

        print "{} microsatellites writen to: {}_{}.txt".format(
            self.mode.title(),
            self.output_prefix,
            self.mode
        )
        print (
            "Calculation of the {} microsatellites finished successfully."
            .format(self.mode)
        )
        print "All calculations finished successfully!\n"

    def _log_normal_result(self, result):
        self.read_lengths_normal.append(result)

    def _log_tumor_result(self, result):
        self.read_lengths_tumor.append(result)

    def _log_unphased_normal_result(self, result):
        self.read_lengths_normal_unphased.append(result)

    def _log_unphased_tumor_result(self, result):
        self.read_lengths_tumor_unphased.append(result)

    def _populate_reference_sets(self, refsets=None):
        for chromosome in self.chromosomes:
            refsetgen = utils.loadcsv(
                self.reference_set+"/reference_set_"+str(chromosome)+"_sorted.txt",
                self.min_microsatellite_length,
                self.max_microsatellite_length,
                self.repeat_units
            )
            refsets = [x for x in refsetgen]
            self.reference_set_dict[chromosome] = refsets
            # get the index positions
            refset_ini_end = [x[1] for x in refsets]
            self.reference_set_ini_end_dict[chromosome] = refset_ini_end
        self.reference_sets = refsets

    def run(self):
        """
        Public method that wraps phased/unphased runs and executes them in
        their multiprocess Pools
        """
        if self.is_phased:
            # Phased Normal run
            self._run_in_pool(
                utils.phased,
                self.normal_bam,
                self.sites,
                self._log_normal_result,
                self.NORMAL
            )
            # Phased Tumor run
            self._run_in_pool(
                utils.phased,
                self.tumor_bam,
                self.sites,
                self._log_tumor_result,
                self.TUMOR
            )

        if self.is_unphased:
            # Unphased Normal run
            self._run_in_pool(
                utils.unphased,
                self.normal_bam,
                self.reference_sets,
                self._log_unphased_normal_result,
                self.NORMAL
            )
            # Unphased Tumor run
            self._run_in_pool(
                utils.unphased,
                self.tumor_bam,
                self.reference_sets,
                self._log_unphased_tumor_result,
                self.TUMOR
            )

        self._conclude_run()

    def _run_in_pool(self,
                     func_to_run,
                     bam_file,
                     sites,
                     logging_method,
                     tumor_type):
        """
        Run phased/unphased MSI detection for normal/tumor bamfiles in a
        multiprocessing Pool
        """
        print "{}: Extracting MS repeats from {} bam file..\n".format(
            self.mode.upper(),
            tumor_type
        )

        if self.number_of_processors == 1:
            logging_method(func_to_run(self, sites, bam_file))
        else:
            pool = mp.Pool(
                self.number_of_processors,
                initializer=utils.multiprocessing_lock_init,
                initargs=(mp.Lock(),)
            )
            for index in np.arange(0, self.number_of_processors):
                if index != (self.number_of_processors - 1):
                    pool.apply_async(
                        func_to_run,
                        args=(
                            self,
                            sites[
                                index *
                                self.chunk_size:(index + 1) *
                                self.chunk_size
                            ],
                            bam_file,
                        ),
                        callback=logging_method
                    )
                else:
                    pool.apply_async(
                        func_to_run,
                        args=(
                            self,
                            sites[
                                index *
                                self.chunk_size: len(sites)
                            ],
                            bam_file,
                        ),
                        callback=logging_method
                    )
            pool.close()
            pool.join()

        print "{}: tumor/case bam file processed correctly..\n".format(
            self.PHASED.upper()
        )

    def _set_fasta_dict(self):
        for index in self.chromosomes:
            self.fasta_dict[index] = "{}chr{}.fa".format(self.fasta_directory, index)

    def _validate_arguments(self):
        """
        Validate the contents of our arg parser argument values
        :raises RuntimeError
        """
        self._check_processors()
        self._check_repeat_units()
        self._check_chromosomes()
        self._check_mode()
        self._check_bams()
        self._check_reference_set()
        self._check_bed_filename()
        self._check_fasta_filename()

    def _write_phased_output(self, outf, common_keys, all_normal, all_tumor):
        for name in common_keys:
            nor = all_normal[name]
            canc = all_tumor[name]
            if isinstance(nor, int) == False and isinstance(canc,
                                                            int) == False:
                if len(nor) >= self.min_coverage and len(
                        canc) >= self.min_coverage:
                    pval_ks = stats.ks_2samp(nor, canc)[1]
                    outf.write(name + "\t" + ",".join(
                        [str(x) for x in nor]) + "\t" + ",".join(
                        [str(x) for x in canc]) + "\t" + str(pval_ks) + "\n")

    def _write_unphased_output(self, outf, common_keys, all_normal, all_tumor):
        for name in common_keys:
            nor = all_normal[name]
            canc = all_tumor[name]
            if isinstance(nor, int) == False and \
                            isinstance(canc, int) == False:
                if len(nor) >= self.min_coverage and \
                                len(canc) >= self.min_coverage:
                    pval = stats.ks_2samp(nor, canc)[1]
                    mo = stats.mode(nor)
                    percentage = (nor == mo).sum() / len(nor)
                    confidence = "high" if percentage >= .7 else "low"
                    outf.write(name + "\t" + ",".join(
                        [str(x) for x in nor]) + "\t" + ",".join(
                        [str(x) for x in canc]) + "\t" + str(
                        pval) + "\t" + confidence + "\n")
