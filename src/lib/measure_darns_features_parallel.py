#!/usr/env python

from argparse import ArgumentParser
import sys
import os
import cPickle as pickle
import pandas as pd
import multiprocessing
import parmap


def measureDarnsFeatures(darnsNumber, regions, hmm, chrom, regionType, features, coverage):
    """
    """
    start, end, center = regions[chrom][darnsNumber]
    name = "_".join([chrom, str(start), str(end), str(center)])

    series = pd.Series(index=["type"] + ["name"] + features)
    series["type"] = regionType
    series["name"] = name

    # measure length
    series["length"] = end - start

    # measure distance to neighbours
    if darnsNumber != 0:
        series["space_upstream"] = abs(start - regions[chrom][darnsNumber - 1][1])  # current start minus previous end
    else:
        series["space_upstream"] = None
    if darnsNumber != len(regions[chrom]) - 1:
        series["space_downstream"] = abs(regions[chrom][darnsNumber + 1][0] - end)  # current end minus next start
    else:
        series["space_downstream"] = None

    # get posterior prob
    # sequence = genome_binary[chrom][start: end]
    # series["post_prob"] = model.retrieveProbabilities(sequence)

    # get read count
    series["read_count"] = coverage[regionType][name].sum()
    # get density
    # I defined it as (sum / length). It is not clear this is the same as in Winter2013
    series["read_density"] = coverage[regionType][name].sum() / (end - start)

    # n. of positive peaks in nucleosome model
    hmm = eval("hmmOutput") if regionType == "DARNS" else eval("hmmOutputP")
    series["n_pospeaks"] = hmm[chrom][start: end].count("+")

    # n. of negative peaks in nucleosome model
    series["n_negpeaks"] = hmm[chrom][start: end].count("-")

    # append
    return series


def main():
    parser = ArgumentParser()

    # Global options
    # positional arguments
    parser.add_argument(dest='input_pickle', type=str, help='Pickle file to load.')
    parser.add_argument(dest='output_pickle', type=str, help='Pickle file to save to.')

    parser.add_argument(dest='regionType', type=str)
    parser.add_argument(dest='chrom', type=str)

    args = parser.parse_args()

    results_dir = "projects/chipmentation/results/periodicity"
    sampleName = "H3K4me3_K562_500k_CM"

    # start work
    features = ["length", "space_upstream", "space_downstream",
                "read_count", "read_density", "n_pospeaks", "n_negpeaks"]

    # load darns
    hmmOutput, DARNS = pickle.load(open(os.path.join(results_dir, sampleName + "_HMMResult.pickle"), "r"))

    # load darnsp
    hmmOutputP, DARNSP = pickle.load(open(os.path.join(results_dir, sampleName + "_HMMResultPermuted.pickle"), "r"))

    # load signal coverage
    cov = pickle.load(open(os.path.join(results_dir, sampleName + "_DARNS.coverage.pickle"), "r"))

    # load type, chrom and iterable with indexes of darns to be measured
    indexes = pickle.load(open(os.path.join(results_dir, sampleName + "_DARNS.coverage.pickle"), "r"))  # ("DARNS", "chr1", [1,2,3])

    if args.regionType == "DARNS":
        regions = DARNS[args.chrom]
        hmm = hmmOutput[args.chrom]
    elif args.regionType == "DARNSP":
        regions = DARNSP[args.chrom]
        hmm = hmmOutputP[args.chrom]

    # for each darn compute features and
    # "reduce" to dataframe
    features = parmap.map(
        measureDarnsFeatures,  # function
        indexes,  # iterator
        regions,  # regions
        hmm,  # hmmOutput
        args.chrom,  # chrom
        args.regionType,  # either DARNS or DARNSP
        features,  # list of features to measure and append to series
        cov  # signal coverage frame
    )

    pickle.dump(features, open(args.output_pickle, "wb"), protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(0)
