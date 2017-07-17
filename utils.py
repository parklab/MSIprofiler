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

