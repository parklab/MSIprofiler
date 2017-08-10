import argparse
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utils import find_repeats_reference


def initialize_parser():
    parser = argparse.ArgumentParser(
        description=(
            'usage: get_reference_set_from_fasta.py '
            '--chr 1 '
            '--notation Ensemble '
            '--min_score 5'
        )
    )
    parser.add_argument(
        '--chr',
        help='Chromosome',
        required=False,
        default=22
    )

    # Optional arguments
    parser.add_argument(
        '--notation',
        help=(
            'Chromosome notation. The two options are: Ensemble (e.g. 1,2,3,'
              '..) or UCSC (e.g. chr1, chr2, chr3,..)'
        ),
        required=False,
        default="Ensemble",
        choices=["Ensemble","UCSC"]
    )
    parser.add_argument(
        '--min_score',
        help=(
            'Minimum score to detect microsatellites in the flanking regions. '
              'Default is 4.'
        ),
        required=False,
        default=5,
        type=int
    )
    parser.add_argument(
        '--max_length',
        help='Maximum microsatellite length. Default is 60.',
        required=False,
        default=60,
        type=int
    )
    return parser


def main():
    parser = initialize_parser()
    args = parser.parse_args()

    fi = "./chrs_fa/chr" + str(args.chr) + ".fa"

    assert os.path.exists(fi), "File: {} doesn't exist.".format(fi)

    #load fasta file:
    # we skip first line, as it contains the chr name
    with open(fi, "rb") as myfile:
        sequence = ''.join(
            [line.strip("\n").upper() for line in myfile.readlines()[1:]]
        )
    chromosome = args.chr

    o = find_repeats_reference(sequence, 0, chromosome, min_score=args.min_score)

    with open("reference_set_"+str(chromosome)+".txt", 'w') as f:
        for ms in o:
            repeats_flanking_left = find_repeats_reference(
                sequence[ms[1]-10:ms[1]],
                0,
                chromosome,
                min_score=args.min_score
            )
            repeats_flanking_right = find_repeats_reference(
                sequence[ms[2]:ms[2]+10],
                0,
                chromosome,
                min_score=args.min_score
            )
            if len(repeats_flanking_left) == 0 and \
                    len(repeats_flanking_right) == 0 and \
					int(ms[5]) <= args.max_length:
                # get the GC contents
                region = sequence[ms[1]-1000:ms[2]+1000]
                GC = region.count("G") + region.count("C")
                GC = round(GC / float(2000), 1)

                base_cmd = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    ms[0], ms[1], ms[2], ms[3], ms[4], ms[5], GC
                )

                if args.notation == "Ensemble":
                    cmd = base_cmd
                else:
                    cmd = "chr{}".format(base_cmd)

                f.write(cmd+"\n")

if __name__ == "__main__":
    main()
