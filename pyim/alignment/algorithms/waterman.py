
import numpy as np
from pyim.alignment.model import Alignment

TRACE_END = 0
TRACE_DIAGONAL = 1
TRACE_UP = 2
TRACE_LEFT = 3

TRACE_LOOKUP = [TRACE_END,       # Max score = 0:   means end of path
                TRACE_DIAGONAL,  # Max score match: means trace diagonal
                TRACE_LEFT,      # Max score up:    means trace left
                TRACE_UP]        # Max score left:  means trace up


def water(query, target, matchScore=5, mismatchScore=-12, gapScore=-4, verbose=False):
    seqA, seqB = query.seq, target.seq
    m, n = len(seqA), len(seqB)

    # Generate DP table and traceback path pointer matrix
    score = np.zeros((m + 1, n + 1))      # DP table
    pointer = np.zeros((m + 1, n + 1))    # Traceback

    # Calculate DP table and mark pointers
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            scoreMatch = score[i - 1][j - 1] + match_score(seqA[i - 1], seqB[j - 1], matchScore, mismatchScore, gapScore)
            scoreUp = score[i][j - 1] + gapScore
            scoreLeft = score[i - 1][j] + gapScore

            stepScores = np.array([0, scoreMatch, scoreUp, scoreLeft])
            maxInd = stepScores.argmax()
            score[i][j] = stepScores[maxInd]
            pointer[i][j] = TRACE_LOOKUP[int(maxInd)]

    # Determine alignment + some basic statistics
    alignA, alignB, startA, startB, endA, endB = traceback(seqA, seqB, score, pointer)
    identity, score, matchStr = alignment_stats(alignA, alignB, matchScore, mismatchScore, gapScore, verbose=verbose)

    return Alignment(query.seqId, startA, endA, query.seq,
                     target.seqId, startB, endB, target.seq,
                     score, identity, matchStr, 'inexact')


def match_score(alpha, beta, matchScore, mismatchScore, gapScore):
    if alpha == beta:
        return matchScore
    elif alpha == '-' or beta == '-':
        return gapScore
    else:
        return mismatchScore


def traceback(seqA, seqB, scoreMat, pointerMat):
    maxI, maxJ = np.unravel_index(scoreMat.argmax(), scoreMat.shape)

    alignmentA = ''
    alignmentB = ''

    # Trace back from maximum
    i, j = maxI, maxJ
    while pointerMat[i][j] != TRACE_END:
        if pointerMat[i][j] == TRACE_DIAGONAL:
            alignmentA += seqA[i - 1]
            alignmentB += seqB[j - 1]
            i -= 1
            j -= 1
        elif pointerMat[i][j] == TRACE_LEFT:
            alignmentA += '-'
            alignmentB += seqB[j - 1]
            j -= 1
        elif pointerMat[i][j] == TRACE_UP:
            alignmentA += seqA[i - 1]
            alignmentB += '-'
            i -= 1

    # Correct orientation of alignments
    alignmentA = alignmentA[::-1]
    alignmentB = alignmentB[::-1]

    return alignmentA, alignmentB, i, j, maxI, maxJ


def alignment_stats(alignmentA, alignmentB, matchScore, mismatchScore, gapScore, verbose=False):
    # Calculate identity, score and aligned sequences
    matchStr = ''
    score, identity = 0, 0

    for i in range(0, len(alignmentA)):
        # Match
        if alignmentA[i] == alignmentB[i]:
            matchStr += '|'
            identity += 1
            score += matchScore

        # Not identical, no gaps
        elif alignmentA[i] != alignmentB[i] and alignmentA[i] != '-' and alignmentB[i] != '-':
            score += mismatchScore
            matchStr += ' '

        # Not identical, one is a gap
        elif alignmentA[i] == '-' or alignmentB[i] == '-':
            matchStr += ' '
            score += gapScore


    identity = float(identity) / len(alignmentA) * 100 if len(alignmentA) > 0 else 0.0

    if verbose:
        print_alignment(alignmentA, alignmentB, matchStr, identity, score)

    return identity, score, matchStr


def print_alignment(alignmentA, alignmentB, matchStr, identity, score):
    print('Identity =', "%3.3f" % identity, 'percent')
    print('Score =', score)
    print(alignmentA)
    print(matchStr)
    print(alignmentB)
