# Isidro Cortes-Ciriano
# Harvard Medical School
# isidrolauscher@gmail.com

import argparse
import models


def main():
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
        '--chromosomes',
        help='Chromosomes to be {}'.format(
            models.MicroSatelliteProfiler.PHASED
        ),
        required=True,
        nargs='+'
    )
    parser.add_argument(
        '--fasta',
        help=(
            'Path to the directory containing the fasta sequences (one per '
            'chromosome). The expected names for the fasta files are e.g. chr1.fa'
        ),
        required=True
    )
    parser.add_argument(
        '--reference_set',
        help=(
            'Path to the directory containing the reference '
            'sets of microsatellites'
        ),
        required=True
    )
    parser.add_argument(
        '--output_prefix',
        help=(
            'Path and prefix for the output files. E.g. '
            'path_to_out_dir/out_prefix'
        ),
        required=True
    )
    parser.add_argument(
        '--mode',
        help='{} or {}'.format(
            models.MicroSatelliteProfiler.PHASED.title(),
            models.MicroSatelliteProfiler.UNPHASED.title()
        ),
        required=True,
        choices=[
            models.MicroSatelliteProfiler.PHASED,
            models.MicroSatelliteProfiler.UNPHASED
        ]
    )
    parser.add_argument(
        '--nprocs',
        help='Number of processes',
        required=True,
        type=int
    )
    parser.add_argument(
        '--rus',
        help=(
            'MS repeat units. Supported from 1 (i.e. mono repeats) to 6 (i.e. '
            'hexarepeats)'
        ),
        required=True,
        nargs='+',
        type=int
    )

    # Optional arguments
    parser.add_argument(
        '--min_MS_length',
        help=(
            'Minimum length of microsatellites to be considered. Minimum '
              'available is 6; default is 10.'
        ),
        required=False,
        default=10,
        type=int
    )
    parser.add_argument(
        '--max_MS_length',
        help=(
            'Maximum length of microsatellites to be considered. Maximum '
              'available is 60; default is 60.'
        ),
        required=False,
        default=60,
        type=int
    )
    parser.add_argument(
        '--mapping_quality',
        help=(
            'Minimum mapping quality. Default is 40.'
        ),
        required=False,
        default=40,
        type=int
    )
    parser.add_argument(
        '--flank_size',
        help=(
            'Minimum length of the flanking regions. Default is 10'
        ),
        required=False,
        default=10,
        type=int
    )
    parser.add_argument(
        '--min_coverage',
        help=(
            'Minimum coverage at each MS locus -both in the case and control '
              'bams-. Default is 10'
        ),
        required=False,
        default=10,
        type=int
    )
    parser.add_argument(
        '--tolerated_mismatches',
        help=(
            'Maximum number of tolerated mismatches in the flanking regions. '
              'Default is 0'
        ),
        required=False,
        default=0,
        type=int
    )

    args = parser.parse_args()
    msi_profiler = models.MicroSatelliteProfiler(args)
    msi_profiler.run()

    print "All calculations finished successfully!\n"

if __name__ == '__main__':
    main()
