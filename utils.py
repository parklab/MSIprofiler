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

