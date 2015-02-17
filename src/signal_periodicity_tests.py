#!/usr/env python

from signal_periodicity import *
import numpy as np
from collections import Counter
import itertools


def test_binarize():
    X = np.array([j for i in [range(-10, 10), range(-11, 9)[::-1]] * 10 for j in i])
    observed = binarize(X)
    assert len(observed) == len(X)

    X = np.zeros(2000)
    observed = binarize(X)
    assert len(observed) == len(X)


def test_getConsensus():
    seq1 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    seq2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    result = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0,
           0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, -1])
    assert all(getConsensus(seq1, seq2) == result)


def test_reduceToOne():
    seq1 = []
    seq2 = [0] * 10
    output = reduceToOne(seq1, seq2)
    assert output is None

    seq1 = []
    seq2 = []
    output = reduceToOne(seq1, seq2)
    assert output == np.array([])

    seq1 = [0] * 10
    seq2 = [0] * 10
    output = reduceToOne(seq1, seq2)
    assert all([i == 4 for i in output])

    seq1 = [0] * 10
    seq2 = [1] * 10
    output = reduceToOne(seq1, seq2)
    assert all([i == 3 for i in output])

    seq1 = [0, 1, -1, 0, 1, -1]
    seq2 = [1, -1, 0, 1, -1, 0]
    output = reduceToOne(seq1, seq2)
    assert all(output == [3, 2, 7, 3, 2, 7])

    seq1 = [0, 1, -1, 0, 1, -1]
    seq2 = [-1, 0, 1, -1, 0, 1]
    output = reduceToOne(seq1, seq2)
    assert all(output == [5, 1, 6, 5, 1, 6])


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

def test_measureDARNS():
    DARNS = []
    d, m = measureDARNS(DARNS)
    assert d == Counter() and m == Counter()
    
    DARNS = [()]
    d, m = measureDARNS(DARNS)
    assert d == Counter() and m == Counter()

    DARNS = [(0, 0), (0, 0)]
    d, m = measureDARNS(DARNS)
    assert d == Counter({0 : 1}) and m == Counter({0: 1})

    DARNS = [(1, 1), (2, 2)]
    d, m = measureDARNS(DARNS)
    assert d == Counter({1 : 1}) and m == Counter({1: 1})

    DARNS = [(1, 10), (2, 11)]
    d, m = measureDARNS(DARNS)
    assert d == Counter({-10 : 1}) and m == Counter({1: 1})

    DARNS = [(0, 9), (10, 19), (20, 29), (30, 39)]
    d, m = measureDARNS(DARNS)
    assert d == Counter({1: 3, 11: 2, 21: 1}) and m == Counter({10: 3, 20: 2, 30: 1})
