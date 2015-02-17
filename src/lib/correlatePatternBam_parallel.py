#!/usr/env python

from argparse import ArgumentParser
import cPickle as pickle
import multiprocessing
import parmap
import numpy as np
from scipy.stats.stats import pearsonr
import os
import HTSeq
import random
import itertools


def correlatePatternBam(feature, pattern, bam, fragmentsize, orientation, duplicates, strand_wise, permutate, step):
    """
    Gets read coverage in a single genomic interval. Returns 1D np.array if strand_wise=False, 2D if True.

    feature=HTSeq.GenomicInterval object.
    bam=HTSeq.BAM_Reader object - Must be sorted and indexed with .bai file!
    fragmentsize=int.
    stranded=bool.
    duplicates=bool.
    """
    # Loop through TSSs, get coverage, append to dict
    # FIXME:
    # select only features in chroms
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX']

    # Initialize empty array for this feature
    if not strand_wise:
        profile = np.zeros(feature.length, dtype=np.float64)
    else:
        profile = np.zeros((2, feature.length), dtype=np.float64)

    # Check if feature is in bam index
    if feature.chrom in chroms:
        # Fetch all alignments in feature window
        alnsInWindow = bam[feature]
        if permutate:
            # randomize each alignment's position in window
            alns = list()
            for aln in alnsInWindow:
                aln.iv.start_d = random.randrange(feature.start, feature.end)
                alns.append(aln)
            alnsInWindow = alns

        for aln in alnsInWindow:
            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue
            aln.iv.length = fragmentsize  # adjust to size

            # get position relative to window
            if orientation:
                if feature.strand == "+" or feature.strand == ".":
                    start_in_window = aln.iv.start - feature.start - 1
                    end_in_window = aln.iv.end - feature.start - 1
                else:
                    start_in_window = feature.length - abs(feature.start - aln.iv.end) - 1
                    end_in_window = feature.length - abs(feature.start - aln.iv.start) - 1
            else:
                start_in_window = aln.iv.start - feature.start - 1
                end_in_window = aln.iv.end - feature.start - 1

            # check fragment is within window; this is because of fragmentsize adjustment
            if start_in_window <= 0 or end_in_window > feature.length:
                continue

            # add +1 to all positions overlapped by read within window
            if not strand_wise:
                profile[start_in_window: end_in_window] += 1
            else:
                if aln.iv.strand == "+":
                    profile[0][start_in_window: end_in_window] += 1
                else:
                    profile[1][start_in_window: end_in_window] += 1

    # Correlate pattern and coverage
    if not strand_wise:
        return binarize(correlatePatternProfile(profile, pattern, step))
    else:
        peaksPos = binarize(correlatePatternProfile(profile[0], pattern, step))
        peaksNeg = binarize(correlatePatternProfile(profile[1], pattern, step))
        return reduceToOne(peaksPos, peaksNeg)


def correlatePatternProfile(profile, pattern, step):
    """
    Fits a sliding window of len(pattern) through a profile and calculates the Pearson
    correlation between the pattern and each window.
    Returns array with correlation values of dimensions (((profile - pattern) + 1) / step).

    profile=np.array.
    pattern=np.array.
    """
    if (((len(profile) - len(pattern)) + 1) / step) <= 0:
        return None
    else:
        R = list()
        i = 0
        while i + len(pattern) <= len(profile):
            R.append(pearsonr(pattern, profile[i:i + len(pattern)])[0])
            i += 1
        return np.nan_to_num(np.array(R))  # replace nan with 0


def binarize(X, debug=False):
    """
    Convert a numerical iterable (X) in a np.array with binary values by finding local maximas (1)
    and setting the rest to zeros (0). Local maximas are required to be higher than 0.05, and have negative
    values between them.

    X=iterable with numerical values.
    """
    X = np.array(X)

    maximas = np.r_[True, X[1:] > X[:-1]] & np.r_[X[:-1] > X[1:], True]
    minimas = np.r_[True, X[1:] < X[:-1]] & np.r_[X[:-1] < X[1:], True]
    prev_max, prev_min = 0, 0

    l = len(X)

    binary = list()

    for i in xrange(l):
        # initialization (to allow peaks on the first base)
        if i == 0:
            if maximas[i] == True and X[i] > 0.05:
                if debug: print("start base is maxima")
                binary.append(1)
            elif minimas[i] == True and X[i] < -0.05:
                if debug: print("start base is minima")
                binary.append(-1)
            else:
                if debug: print("start base is zero")
                binary.append(0)
        # initialization (to allow peaks on the first base)
        elif i == l - 1:
            if maximas[i] == True and X[i] > 0.05:
                if debug: print("end base is maxima")
                binary.append(1)
            elif minimas[i] == True and X[i] < -0.05:
                if debug: print("end base is minima")
                binary.append(-1)
            else:
                if debug: print("end base is zero")
                binary.append(0)

        elif maximas[i] == True and X[i] > 0.05 and any([n < 0 for n in X[prev_max: i]]):  # watch out for namespace pollution with np.any
            # find next minima
            for j in xrange(i + 1, len(X)):
                if maximas[j] == True and X[j] > 0.05 or j == l - 1:
                    next_max = j
                    break

            # if from i to next_max there is neg, then is max
            if any([n < 0 for n in X[i: next_max]]):
                if debug: print("no maxima after, so is maxima")
                binary.append(1)
                prev_max = i
            else:
                if debug: print("maxima after, so is not maxima")
                binary.append(0)
        elif minimas[i] == True and X[i] < -0.05 and any([n > 0 for n in X[prev_min: i]]):
            # find next minima
            for j in xrange(i + 1, len(X)):
                if minimas[j] == True and X[j] < -0.05 or j == l - 1:
                    next_min = j
                    break

            # if from i to next_neg there is pos, then is min
            if any([n > 0 for n in X[i: next_min]]):
                if debug: print("no minima after, so is minima")
                binary.append(-1)
                prev_min = i
            else:
                if debug: print("minima after, so is not minima")
                binary.append(0)
        else:
            if debug: print("no maxima or minima")
            binary.append(0)
    return np.array(binary)


def reduceToOne(seq1, seq2):
    """
    Given two trinary (1, 0, -1) sequences, return one sequence where pairs of values map to an integer (0-8).

    seq1,seq2=iterables with trinary (1, 0, -1) numerical values.
    """
    if not len(seq1) == len(seq2):
        return None
    else:
        mapping = {v: k for k, v in enumerate(itertools.product([1, 0, -1], repeat=2))}

        binary = list()
        for i in xrange(len(seq1)):
            binary.append(mapping[(seq1[i], seq2[i])])
        return np.array(binary)


if __name__ == '__main__':
    parser = ArgumentParser(description='correlatePatternBam_parallel.py',
                            usage='python correlatePatternBam_parallel.py [options] input_pickle output_pickle pattern_pickle bam_file'
                            )
    # Global options
    # positional arguments
    parser.add_argument(dest='input_pickle', type=str, help='Pickle file to load.')
    parser.add_argument(dest='output_pickle', type=str, help='Pickle file to save to.')
    parser.add_argument(dest='pattern_pickle', type=str, help='Pickle file with pattern.')
    parser.add_argument(dest='bam_file', type=str, help='Bam file.')

    # optional arguments
    parser.add_argument('--step', dest='step', type=int, default=1)
    parser.add_argument('--fragment-size', dest='fragment_size', type=int, default=1)
    parser.add_argument('--strand-wise', dest='strand_wise', action='store_true')
    parser.add_argument('--duplicates', dest='duplicates', action='store_true')
    parser.add_argument('--permute', dest='permute', action='store_true')
    parser.add_argument('--orientation', dest='orientation', action='store_true')

    args = parser.parse_args()

    # Read pickle with windows
    windows = pickle.load(open(args.input_pickle, 'r'))  # dict-like list of tuples

    # Load pattern
    pattern = pickle.load(open(args.pattern_pickle, 'r'))  # numpy.array with pattern

    # Make bam object
    bam = HTSeq.BAM_Reader(os.path.abspath(args.bam_file))

    # Process in parallel and serialize result
    # parallel process and reduce to Counter
    coverage = dict(
        zip(
            [tup[0] for tup in windows],
            parmap.map(
                correlatePatternBam,
                [tup[1] for tup in windows],  # similar to windows.values()
                pattern,
                bam,
                args.fragment_size,
                args.orientation,
                args.duplicates,
                args.strand_wise,
                args.permute,
                args.step
            )
        )
    )

    # Serialize
    pickle.dump(coverage, open(args.output_pickle, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
