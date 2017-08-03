
import sys
#add option for notation: ensemble (1,2,3) o ucsc (chr1, chr2,..)
import utils

chr = sys.argv[1]
fi = "chrs_fa/chr"+str(chr)+".fa"


def fetch_reference_sets():
    #load fasta file:
    # we skip first line
    sequence = load_fasta_file()
    o = utils.find_repeats(sequence, 0, chr, min_score=5)
    write_reference_set_file(o, sequence)


def load_fasta_file():
    with open(fi, "rb") as myfile:
        sequence = ''.join(
            [line.strip("\n").upper() for line in myfile.readlines()[1:]]
        )
    return sequence

def write_reference_set_file(o, sequence):
    ## XX add to REFINE THE SET
    ## add to find location

    with open("reference_set_" + str(chr) + "_up_to_5rep_in_flank.txt",
              'w') as f:
        for ms in o:
            repeats_flanking_left = utils.find_repeats(
                sequence[ms[1] - 10:ms[1]],
                0,
                chr,
                min_score=4
            )
            repeats_flanking_right = utils.find_repeats(
                sequence[ms[2]:ms[2] + 10],
                0,
                chr,
                min_score=4
            )
            if len(repeats_flanking_left) == 0 and \
                            len(repeats_flanking_right) == 0:
                cmd = "{}\t{}\t{}\t{}\t{}\t{}".format(
                    ms[0], ms[1], ms[2], ms[3], ms[4], ms[5]
                )
                f.write(cmd + "\n")

if __name__ == '__main__':
    fetch_reference_sets()
