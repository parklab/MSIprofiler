import bisect
import csv
import string
import pysam
import numpy as np

MATCH_SCORE = 1  # score for a match
MISMATCH_SCORE = -6  # score penalty for a mismatch
FAIL_SCORE = -1  # score value to stop searching
MIN_SCORE = 5  # minimum score value to pass. The minimum length of MS repeats


# https://stackoverflow.com/questions/212358/binary-search-bisection-in-python
def binary_search(a, x, lo=0,
                  hi=None):  # can't use a to specify default for hi
    hi = hi if hi is not None else len(a)  # hi defaults to len(a)
    pos = bisect.bisect_left(a, x, lo, hi)  # find insertion position
    return (pos if pos != hi and a[pos] == x else -1)


def find_repeats(seq, flank_size, repeat_units, min_score=MIN_SCORE):
    bases = len(seq)
    flank_size = flank_size - 1
    # save output as a list of lists
    out = [];
    exclude = set()  # use sets: they are much faster to apply 'in'
    for ru in repeat_units:
        positions_motif = range(0, ru)
        nb_positions_motif = len(positions_motif)
        not_found = True
        base = flank_size
        while base < bases - flank_size:  # and base not in exclude:
            if base in exclude:
                base += 1
                continue
            elif not_found:
                test_pos = base + ru
                current_pos = base
            else:
                current_pos = base
                not_found = True
                test_pos = test_pos + ru
            pos_in_motif = 0
            score = 0;
            depth = 0;
            keep = 0
            max_observed_score = 0
            scores = []
            while (
                        (test_pos) < (
                        bases - flank_size)) and score > FAIL_SCORE and \
                            test_pos not in exclude:
                match = (seq[current_pos + pos_in_motif] == seq[test_pos])
                if match:
                    test_pos += 1
                    pos_in_motif = positions_motif[
                        (pos_in_motif + 1) % nb_positions_motif]
                    score += MATCH_SCORE
                    scores.append(score)
                    depth = 0
                else:
                    score += MISMATCH_SCORE
                    scores.append(score)
                    pos_in_motif = positions_motif[
                        (pos_in_motif + 1) % nb_positions_motif]
                    if score > FAIL_SCORE and depth < 5:
                        depth += 1
                        test_pos += 1
                # keep track of the best observed score
                if score > max_observed_score:
                    max_observed_score = score
                    # debugging
                    # print "RU current_pos  pos_in_motif  test_pos   bases, score"
                    # print ru, current_pos, pos_in_motif, test_pos, bases,score,"\n" #, current_pos + pos_in_motif, bases
            if max_observed_score >= min_score:  # and test_pos <= bases-flank_size:
                mm = scores.index(max(scores))
                mm = mm + ru
                if base + mm < (
                            bases - flank_size):  # repeat not overlapping flanking region
                    out.append([ru, base, base + mm, seq[base:base + mm + 1]])
                    not_found = False
                exclude.update(range(base, base + mm + 1))
                test_pos = base + mm
                base = test_pos
            else:
                pass
            base += 1
        else:
            base += 1
    return out


def find_repeats_target(seq, flank_size, repeat_units):
    # rus = {1,2,3,4} # repeat units considered ## take this from the CONF FILE
    MATCH_SCORE = 1  # score for a match
    MISMATCH_SCORE = -6  # score penalty for a mismatch
    FAIL_SCORE = -1  # score value to stop search at a given current_pos
    MIN_SCORE = 5  ##minimum score value to pass

    bases = len(seq)  # number of bases in the input sequence

    # save output as a list of lists
    out = []
    # use sets: they are much faster with 'in'[]#np.array([], dtype='int')
    exclude = set()

    for ru in [repeat_units]:
        positions_motif = range(0, ru)
        nb_positions_motif = len(positions_motif)

        # note that the flank is a range, whether python is zero-based
        not_found = True
        base = flank_size

        while base < bases - flank_size:  # and base not in exclude:
            if base in exclude:
                base += 1
                continue
            elif not_found:
                test_pos = base + ru
                current_pos = base
            else:
                current_pos = base
                not_found = True
                test_pos = test_pos + ru

            pos_in_motif = 0  # update_current_pos(ru)
            score = 0
            depth = 0

            max_observed_score = 0
            scores = []
            while ((test_pos) < (bases - flank_size)) and \
                            score > FAIL_SCORE and \
                            test_pos not in exclude:  # XX the minus one check
                match = (seq[current_pos + pos_in_motif] == seq[test_pos])

                if match:
                    test_pos += 1
                    pos_in_motif = positions_motif[
                        (pos_in_motif + 1) % nb_positions_motif]
                    score += MATCH_SCORE
                    scores.append(score)
                    depth = 0

                # no mismatch: check for N, insertions, deletions and missense
                else:
                    score += MISMATCH_SCORE
                    scores.append(score)
                    pos_in_motif = positions_motif[
                        (pos_in_motif + 1) % nb_positions_motif]

                    if score > FAIL_SCORE and depth < 5:
                        depth += 1
                        test_pos += 1

                # keep track of the best observed score
                if score > max_observed_score:
                    max_observed_score = score
            test_pos = test_pos  # - depth
            if max_observed_score >= MIN_SCORE:
                mm = scores.index(max(scores))
                mm = mm + ru
                out.append([ru, base, base + mm,
                            seq[base:base + mm + 1]])
                not_found = False
                exclude.update(range(base, base + mm + 1))
                test_pos = base + mm
                base = test_pos

            # increment base
            base += 1
        else:
            base += 1
    return out


def loadcsv(filename, criterion1, criterion2, repeat_units):
    with open(filename, "rb") as csvfile:
        datareader = csv.reader(csvfile, delimiter="\t")
        for row in datareader:
            if int(row[5]) >= criterion1 and int(row[5]) <= criterion2 and \
                            int(row[4]) in repeat_units:
                yield row


def phased(msi_obj, sites, bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        dict_out = {}
        for site in sites:
            start = int(site[1])
            end = int(site[2])
            chr = site[0]
            #chr = str(site[0])
            bases = [site[3], site[4]]
            reads = [
                read for read in bamfile.fetch(
                    chr,
                    start,
                    end,
                    multiple_iterators=True
                )
                ]  ## keep this as is, do not put conditions inside here

            reads = [
                read for read in reads if read.is_proper_pair and
                read.is_duplicate == False and
                read.mapping_quality >= msi_obj.mapping_quality
                ]
            if len(reads) > msi_obj.min_coverage:
                for read in reads:
                    read_sequence = read.seq
                    # read_sequence = read.query_alignment_sequence
                    reps = find_repeats(
                        read_sequence,
                        msi_obj.flank_size,
                        msi_obj.repeat_units
                    )
                    if len(reps) > 0:
                        # get the SNP allele in this read
                        start_read = read.reference_start
                        end_read = read.reference_end
                        aligned_pos = read.get_reference_positions(
                            full_length=True
                        )  # True) reports none for  soft-clipped positions
                        try:
                            idx = aligned_pos.index(start)
                        except:
                            continue
                        snp_read = read_sequence[idx]
                        if snp_read not in bases:
                            continue
                        for microsatellite in reps:
                            rs = microsatellite[1]
                            re = microsatellite[2]
                            difference = re - rs + 1
                            # use the reference set here to get the
                            # position on the right
                            ini = start_read + rs  #
                            chr=string.strip(chr,"chr")
                            idx2 = binary_search(
                                msi_obj.reference_set_ini_end_dict[chr],
                                str(ini + 1)
                            )

                            if idx2 == -1:
                                continue
                            refset_now = msi_obj.reference_set_dict[chr][idx2]
                            diff_ref = int(refset_now[2]) - int(
                                refset_now[1]) + 1
                            with pysam.FastaFile(
                                    filename=msi_obj.fasta_dict[chr]
                            ) as fasta_file:
                                flank_right_ref = fasta_file.fetch(
                                    "chr" + str(site[0]), ini + diff_ref, #XX
                                    ini + diff_ref + msi_obj.flank_size).upper()
                                flank_left_ref = fasta_file.fetch(
                                    "chr" + str(site[0]), ini - #XX
                                    msi_obj.flank_size,
                                    ini).upper()
                            posfl = (start_read + rs - msi_obj.flank_size)
                            if posfl >= start_read:
                                flank_left = read_sequence[
                                             rs - msi_obj.flank_size:rs]
                                mismatches_left = sum(a != b for a, b in
                                                      zip(flank_left,
                                                          flank_left_ref))
                            else:
                                mismatches_left = 10000
                            posflr = start_read + re + msi_obj.flank_size
                            if posflr <= end_read:
                                flank_right = read_sequence[
                                              re + 1:re + 1 +
                                                     msi_obj.flank_size]
                                mismatches_right = sum(a != b for a, b in
                                                       zip(flank_right,
                                                           flank_right_ref))
                            else:
                                mismatches_right = 10000
                            mismatches = mismatches_left + mismatches_right
                            if mismatches <= msi_obj.tolerated_mismatches:
                                key_now = site[0] + "\t" + str(ini) + "\t" + \
                                          refset_now[3] + "\t" + refset_now[
                                              4] + "\t" + refset_now[
                                              5] + "\t" + refset_now[
                                              6] + "\t" + snp_read + "\t" + \
                                          str(site[1])
                                if dict_out.has_key(key_now):
                                    dict_out[key_now] = np.append(
                                        dict_out[key_now], difference)
                                else:
                                    dict_out[key_now] = difference

    return dict_out


def unphased(msi_obj, sites, bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        dict_out = {}
        visited_reads = []
        for site in sites:
            start = int(site[1])
            #end = int(site[2]) + 1
            chr = site[0]
            ru = int(site[4])
            reads = [read for read in bamfile.fetch(chr, start, start + 1,
                                                    multiple_iterators=True)]
            reads = [read for read in reads if
                     read.is_proper_pair and read.is_duplicate == False and
                     read.mapping_quality >= msi_obj.mapping_quality]
            if len(reads) > msi_obj.min_coverage:
                for read in reads:
                    start_read = read.reference_start;
                    end_read = read.reference_end
                    read_sequence = read.seq
                    reps =  find_repeats_target(read_sequence,
                                               msi_obj.flank_size,
                                               ru)
                    if len(reps) > 0:
                        aligned_pos = read.get_reference_positions(
                            full_length=True)
                        try:
                            idx = aligned_pos.index(start)
                        except:
                            continue
                        for microsatellite in reps:
                            ru = microsatellite[0];
                            rs = microsatellite[1];
                            re = microsatellite[2]
                            if start != start_read + rs + 1:  # do not consider if there are ins/del upstream of the repeat
                                continue
                            difference = re - rs + 1
                            chr=string.strip(chr,"chr")
                            # get flinking sequence from reference
                            fasta_file = pysam.FastaFile(
                                filename=msi_obj.fasta_dict[chr]
                            )
                            flank_left_ref = fasta_file.fetch(
                                "chr" + chr, #string.strip(chr,"chr"), #X
                                start_read + rs -
                                msi_obj.flank_size,
                                start_read + rs).upper()
                            flank_right_ref = fasta_file.fetch("chr" + chr, #string.strip(chr,"chr"), #X
                                                               int(site[
                                                                       2]) - 1,
                                                               int(site[
                                                                       2]) - 1
                                                               + msi_obj.flank_size).upper()
                            # get flinking sequence from the reads
                            posfl = (start_read + rs - msi_obj.flank_size)
                            if posfl >= start_read:
                                flank_left = read_sequence[
                                             rs - msi_obj.flank_size:rs];
                                mismatches_left = sum(a != b for a, b in
                                                      zip(flank_left,
                                                          flank_left_ref))
                            else:
                                flank_left = "";
                                mismatches_left = 10000
                            posflr = start_read + re + msi_obj.flank_size
                            if posflr <= end_read:
                                flank_right = read_sequence[
                                              re:re + msi_obj.flank_size];
                                mismatches_right = sum(a != b for a, b in
                                                       zip(flank_right,
                                                           flank_right_ref))
                            else:
                                flank_right = "";
                                mismatches_right = 10000
                            mismatches = mismatches_left + mismatches_right
                            if mismatches <= msi_obj.tolerated_mismatches:
                                key_now = site[0] + "\t" + site[1] + "\t" + \
                                          site[
                                              2] + "\t" + site[3] + "\t" + \
                                          site[4] + "\t" + \
                                          site[5] + "\t" + site[6]
                                if dict_out.has_key(key_now):
                                    dict_out[key_now] = np.append(
                                        dict_out[key_now], difference)
                                else:
                                    dict_out[key_now] = difference
        return dict_out


def multiprocessing_lock_init(l):
    global lock
    lock = l

rus={1,2,3,4,5,6} # repeat units considered
match_score = 1 # score for a match
mismatch_score = -6 # score penalty for a mismatch # it detects min length of abs(mismatch_score) + 1
fail_score = -1 # score value to stop search at a given current_pos
#min_score = 5 ##minimum score value to pass; which corresponds to length 6 for mono


def find_repeats_reference(seq,flank_size,chromo,min_score):
    cigar = seq #XX change
    bases=len(seq) # number of bases in the input sequence

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
                #print base, "   ", score, max_observed_score, seq[base:test_pos], depth, exclude, bases-flank_size, "RU",ru
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

            test_pos = test_pos #- depth
            if max_observed_score >= min_score:# and test_pos <= bases-flank_size:
                mm = scores.index(max(scores))
                mm = mm +ru
                seq_now = seq[base:base+mm+1]
                len_seq = len(seq_now)
                if len_seq <= 60:
                    out.append( [chromo,base+1, base+mm+1, seq_now,ru, len_seq])#test_pos]] ) ## chr, start, end, seq, ru, diff
                not_found = False
                exclude.update(range(base,base+mm+1)) #np.unique(np.append(exclude,np.arange(base,test_pos)))
                test_pos = base + mm  ##X
                base = test_pos
            # increment base
            base+=1
        else:
            base+=1
    return out
