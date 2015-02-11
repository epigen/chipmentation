#!/usr/env python

from argparse import ArgumentParser
import cPickle as pickle
import multiprocessing
import parmap
import numpy as np
from scipy.stats.stats import pearsonr


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
        return np.array(R)


if __name__ == '__main__':
    parser = ArgumentParser(description='correlation_parallel.py',
                            usage='python correlation_parallel.py [options] input_pickle output_pickle pattern_pickle'
                            )

    # Global options
    # positional arguments
    parser.add_argument(dest='input_pickle', type=str, help='Pickle file to load.')
    parser.add_argument(dest='output_pickle', type=str, help='Pickle file to save to.')
    parser.add_argument(dest='pattern_pickle', type=str, help='Pickle file with pattern.')

    # optional arguments
    parser.add_argument('--step', dest='step', type=int, default=1)

    args = parser.parse_args()

    # Read pickle with windows
    windows = pickle.load(open(args.input_pickle, 'r'))  # dict-like list of tuples

    # Load pattern
    pattern = pickle.load(open(args.pattern_pickle, 'r'))

    # Process in parallel and serialize result
    # parallel process and reduce to Counter
    coverage = dict(
        zip(
            [tup[0] for tup in windows],  # similar to windows.keys()
            parmap.map(
                correlatePatternProfile,
                [tup[1] for tup in windows],  # similar to windows.values()
                pattern,
                args.step
            )
        )
    )

    # Serialize
    pickle.dump(coverage, open(args.output_pickle, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
