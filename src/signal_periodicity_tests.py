#!/usr/env python

from signal_periodicity import *
import numpy as np


def test_getConsensus():
    seq1 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    seq2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    result = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0,
           0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, -1])
    assert all(getConsensus(seq1, seq2) == result)


def test_concatenateBinary():
    # one item
    binaryDict = OrderedDict({("chr1", 1, 10): np.zeros(10)})
    assert all(concatenateBinary(binaryDict, 1).values()[0] == binaryDict.values()[0])

    # two items
    binaryDict = OrderedDict(sorted({("chr1", 0, 10): np.zeros(9), ("chr1", 9, 19): np.ones(9)}.items()))
    assert all(concatenateBinary(binaryDict, 1).values()[0] == np.hstack([v for n, v in sorted(binaryDict.items())]))

    binaryDict = OrderedDict(sorted({("chr1", 0, 10): np.zeros(8), ("chr1", 8, 18): np.ones(8)}.items()))
    assert all(concatenateBinary(binaryDict, 2).values()[0] == np.hstack([v for n, v in sorted(binaryDict.items())]))

    # different window and pattern length
    patternLength = 25
    window = 1000
    binaryDict = OrderedDict(sorted({
        ("chr1", 1000, 2000): np.zeros(window - patternLength),
        ("chr1", 1975, 2975): np.ones(window - patternLength),
        ("chr1", 2950, 3950): np.ones(window - patternLength)
    }.items()))
    assert all(concatenateBinary(binaryDict, patternLength).values()[0] == np.hstack([v for n, v in sorted(binaryDict.items())]))

    # more windows
    window = 1000
    patternLength = 50
    binaryDict = OrderedDict(sorted({
        ("chr1", 1, 1000): np.zeros(window - patternLength),
        ("chr1", 951, 1950): np.zeros(window - patternLength),
        ("chr1", 1901, 2900): np.zeros(window - patternLength),
        ("chr1", 2851, 3850): np.ones(window - patternLength)
        }.items()))
    # bad test. Windows not sorted in right side
    # assert all(concatenateBinary(binaryDict, patternLength).values()[0] == np.hstack([v for n, v in sorted(binaryDict.items())]))

    # Check length of concatenation is end - start - pattern of all intervals
    assert len(concatenateBinary(binaryDict, patternLength).values()[0]) == 3850 - patternLength


def test_getDARNS():
    # no sequence
    sequence = ""
    DARNS = getDARNS(sequence)
    assert DARNS == list()

    # no darn
    sequence = "B" * 7
    DARNS = getDARNS(sequence)
    assert DARNS == list()

    # 1bp darn
    sequence = "BBBQBBBB"
    DARNS = getDARNS(sequence)
    assert len(DARNS) == 1
    assert DARNS[0][0] == 3 and DARNS[0][1] == 4

    # one DARN
    sequence = "BBBQ!Qqwe"
    DARNS = getDARNS(sequence)
    assert len(DARNS) == 1
    assert DARNS[0][0] == 3 and DARNS[0][1] == 8

    # DARN until the end
    sequence = "BBBQ!QqwB"
    DARNS = getDARNS(sequence)
    assert len(getDARNS(sequence)) == 1
    assert getDARNS(sequence)[0][0] == 3 and getDARNS(sequence)[0][-1] == 7

    # Two DARNS
    sequence = "BBBQ!QqweBBBQ!QqweB"
    DARNS = getDARNS(sequence)
    assert len(DARNS) == 2
    assert DARNS[0][0] == 3 and DARNS[0][1] == 8 and DARNS[1][0] == 12 and DARNS[1][-1] == 17
