# Isidro Cortes-Ciriano
# Harvard Medical School
# isidrolauscher@gmail.com

import argparse

from .models import MicroSatelliteProfiler

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
        'Minimum coverage at each MS locus -both in the case and control bams.'
        'Default is 10'
    ),
    required=False,
    default=10,
    type=int
)
parser.add_argument(
    '--tolerated_mismatches',
    help=('Maximum number of tolerated mismatches in the flanking regions. '
          'Default is 0'),
    required=False,
    default=0,
    type=int
)

args = parser.parse_args()
msi_profiler = MicroSatelliteProfiler(args)


# -----------------------------------------------------------------------------
# PHASED
# -----------------------------------------------------------------------------
if msi_profiler.mode in ['both', 'phased']:
    msi_profiler.run_phased()

# -----------------------------------------------------------------------------
# UNPHASED
# -----------------------------------------------------------------------------
if args.mode in ['both', 'unphased']:
    msi_profiler.run_unphased()

print "All calculations finished successfully!\n"
