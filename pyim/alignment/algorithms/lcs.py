
def longest_common_subsequence(seqA, seqB, matchFunc=None):
    ja, jb, n = -1, -1, 0

    if seqA == [] or seqB == []:
        return ja, jb, n

    if matchFunc is None:
        matchFunc = _match

    s = 1
    l = len(seqA) + len(seqB) - 1
    ia, ib = len(seqA) - 1, 0

    for k in range(l):
        nCurrent = 0

        for r in range(s):
            if matchFunc(seqA[ia + r], seqB[ib + r]):
                nCurrent += 1

                if nCurrent > n:
                    ja = ia + r - nCurrent + 1
                    jb = ib + r - nCurrent + 1
                    n = nCurrent
            else:
                nCurrent = 0

        if k < min(len(seqA), len(seqB)) - 1:
            ia -= 1
            s += 1
        elif k > l - min(len(seqA), len(seqB)) - 1:
            ib += 1
            s -= 1
        elif ia > 0:
            ia -= 1
        else:
            ib += 1

    return ja, jb, n

def _match(a, b):
    return a == b