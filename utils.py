MATCH_SCORE = 1  # score for a match
MISMATCH_SCORE = -6  # score penalty for a mismatch
FAIL_SCORE = -1  # score value to stop searching
# detected is min_score + 4

# https://stackoverflow.com/questions/212358/binary-search-bisection-in-python
def binary_search(a, x, lo=0,
                  hi=None):  # can't use a to specify default for hi
    hi = hi if hi is not None else len(a)  # hi defaults to len(a)
    pos = bisect.bisect_left(a, x, lo, hi)  # find insertion position
    return (pos if pos != hi and a[pos] == x else -1)


def find_repeats(seq, flank_size, rus, min_score=4):
    bases = len(seq)
    flank_size -= 1
    # save output as a list of lists
    out = []
    exclude = set()  # use sets: they are much faster to apply 'in'
    for ru in rus:
        ru = int(ru)
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
            score = 0
            depth = 0
            max_observed_score = 0
            scores = []
            while ((test_pos) < (
                        bases - flank_size)) and score > FAIL_SCORE and test_pos not \
                    in exclude:  # XX the minus one check
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

            if max_observed_score >= min_score:  # and test_pos <=
                # bases-flank_size:
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


def loadcsv(filename, criterion1, criterion2):
    with open(filename, "rb") as csvfile:
        datareader = csv.reader(csvfile, delimiter="\t")
        for row in datareader:
            if int(row[5]) >= criterion1 and int(row[5]) <= criterion2:
                yield row


def phased(bam_path,
           fastafile,
           flank_size,
           mapping_quality,
           min_coverage,
           refset,
           refset_ini_end,
           rus,
           sites,
           tolerated_mismatches):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    dict_out = {}

    for site in sites:
        start = int(site[1])
        end = int(site[2])
        chr = str(site[0])
        bases = [site[3], site[4]]
        reads = [read for read in bamfile.fetch(chr, start, end,
                                                multiple_iterators=True)]  ## keep this as is, do not put conditions inside here
        reads = [read for read in reads
                 if read.is_proper_pair
                 and not read.is_duplicate
                 and read.mapping_quality >= mapping_quality]

        if len(reads) > min_coverage:
            for read in reads:
                read_sequence = read.seq
                # read_sequence = read.query_alignment_sequence
                reps = find_repeats(read_sequence, flank_size, rus)
                if len(reps) > 0:
                    # get the SNP allele in this read
                    start_read = read.reference_start
                    end_read = read.reference_end
                    aligned_pos = read.get_reference_positions(
                        full_length=True)  # True) reports none for soft-clipped positions
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

                        # use the reference set here to get the position on the right
                        ini = start_read + rs  #
                        idx2 = binary_search(refset_ini_end, str(ini + 1))

                        if idx2 == -1:
                            continue
                        refset_now = refset[idx2]
                        diff_ref = int(refset_now[2]) - int(refset_now[1]) + 1
                        flank_right_ref = fastafile.fetch("chr" + str(site[0]),
                                                          ini + diff_ref,
                                                          ini + diff_ref + flank_size).upper()
                        flank_left_ref = fastafile.fetch("chr" + str(site[0]),
                                                         ini - flank_size,
                                                         ini).upper()
                        posfl = (start_read + rs - flank_size)
                        if posfl >= start_read:
                            flank_left = read_sequence[
                                         rs - flank_size:rs]
                            mismatches_left = sum(a != b for a, b in
                                                  zip(flank_left,
                                                      flank_left_ref))
                        else:
                            mismatches_left = 10000
                        posflr = start_read + re + flank_size
                        if posflr <= end_read:
                            flank_right = read_sequence[
                                          re + 1:re + 1 + flank_size]
                            mismatches_right = sum(a != b for a, b in
                                                   zip(flank_right,
                                                       flank_right_ref))
                        else:
                            mismatches_right = 10000
                        mismatches = mismatches_left + mismatches_right

                        if mismatches <= tolerated_mismatches:
                            key_now = site[0] + "\t" + str(
                                ini) + "\t" + snp_read + "\t" + str(site[1])
                            if dict_out.has_key(key_now):
                                dict_out[key_now] = np.append(
                                    dict_out[key_now], difference)
                            else:
                                dict_out[key_now] = difference
    bamfile.close()
    return dict_out


