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
    binaryDict = {"chr1_1_10": np.zeros(10)}
    assert all(concatenateBinary(binaryDict, 1).values()[0] == binaryDict.values()[0])

    # two items
    binaryDict = {"chr1_0_10": np.zeros(9), "chr1_9_19": np.ones(9)}
    assert all(concatenateBinary(binaryDict, 1).values()[0] == np.hstack([v for n, v in sorted(binaryDict.items())]))

    binaryDict = {"chr1_0_10": np.zeros(8), "chr1_8_18": np.ones(8)}
    assert all(concatenateBinary(binaryDict, 2).values()[0] == np.hstack([v for n, v in sorted(binaryDict.items())]))

    # different window and pattern length
    patternLength = 25
    window = 1000
    binaryDict = {"chr1_1000_2000": np.zeros(window - patternLength),
                  "chr1_1975_2975": np.ones(window - patternLength + 1),
                  "chr1_2950_3950": np.ones(window - patternLength + 1)
                  }
    assert all(concatenateBinary(binaryDict, patternLength).values()[0] == np.hstack([v for n, v in sorted(binaryDict.items())]))

    # more windows
    window = 1000
    patternLength = 50
    binaryDict = {"chr1_1_1000": np.zeros(window - patternLength),
                  "chr1_951_1950": np.zeros(window - patternLength),
                  "chr1_1901_2900": np.zeros(window - patternLength),
                  "chr1_2851_3850": np.ones(window - patternLength)
                  } 
    # bad test. Windows not sorted in right side
    # assert all(concatenateBinary(binaryDict, patternLength).values()[0] == np.hstack([v for n, v in sorted(binaryDict.items())]))

    # Check length of concatenation is end - start - pattern of all intervals
    assert len(concatenateBinary(binaryDict, patternLength).values()[0]) == 3850 - patternLength
