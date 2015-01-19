#!/usr/env python

from argparse import ArgumentParser
import os, re
from pybedtools import BedTool
import HTSeq
import numpy as np
import pandas as pd
import string

import rpy2.robjects as robj # for ggplot in R
import rpy2.robjects.pandas2ri # for R dataframe conversion

import itertools
import cPickle as pickle

def makeWindows(windowWidth, genome):
    """Generate 1kb windows genome-wide."""
    w = BedTool.window_maker(BedTool(), genome=genome, w=windowWidth)
    windows = dict()
    for interval in w:
        feature = HTSeq.GenomicInterval(
            interval.chrom,
            interval.start,
            interval.end
        )
        name = string.join(interval.fields, sep="_")
        windows[name] = feature

    return windows


def distances(bam, intervals, fragmentsize, duplicates=True):
    """ Gets read coverage in bed regions, returns dict with region:count.
    bam - Bam object from HTSeq.BAM_Reader.
    intervals - dict with HTSeq.GenomicInterval objects as values.
    fragmentsize - integer.
    duplicates - boolean.
    """
    dists = dict()
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    for name, feature in intervals.items():
        if feature.chrom not in chroms:
            continue
        # Fetch alignments in feature window
        for aln1, aln2 in itertools.combinations(bam[feature], 2):
            # check if duplicate
            if not duplicates and (aln1.pcr_or_optical_duplicate or aln2.pcr_or_optical_duplicate):
                continue
            # adjust fragment to size
            aln1.iv.length = fragmentsize
            aln2.iv.length = fragmentsize
            dist = abs(aln1.iv.start - aln2.iv.start)

            # add +1 to dict
            if dist not in dists.keys():
                dists[dist] = 1
            else:
                dists[dist] += 1
    return dists


def main(args):
    args.plots_dir = os.path.abspath(args.plots_dir)

    # Get sample names
    names = list()
    for bam in args.bamfiles:
        names.append(re.sub("\.bam", "", os.path.basename(bam)))

    # Get genome-wide windows
    print("Making %ibp windows genome-wide" % args.window_width)
    windows = makeWindows(args.window_width, args.genome)

    # Loop through all signals, compute distances, plot
    for bam in xrange(len(args.bamfiles)):
        print("Sample " + names[bam])
        # Load bam
        bamfile = HTSeq.BAM_Reader(os.path.abspath(args.bamfiles[bam]))
        # Get dataframe of bam coverage in bed regions, append to dict
        dists = distances(bamfile, windows, args.fragment_size, args.duplicates)

        pickle.dump(dists, open(names[bam] + ".counts.pickle", "wb"), protocol = pickle.HIGHEST_PROTOCOL)

        x = range(50, 200)
        y = [dists[i] for i in x]
        plt.scatter(x, y)
        p20 = np.poly1d(np.polyfit(x, y, 20))
        p50 = np.poly1d(np.polyfit(x, y, 50))
        p100 = np.poly1d(np.polyfit(x, y, 100))

        plt.plot(x, p20(x), '-', x, p50(x), '--', x, p100(x), '.')
        plt.savefig(os.path.join(args.plots_dir, names[bam] + ".fit.pdf"))
        plt.close()

        x = np.array(range(50, 200))
        y = np.array([dists[i] for i in x])

        from scipy.interpolate import spline
        x_smooth = np.linspace(x.min(), x.max(), 200)
        y_smooth = spline(x, y, x_smooth)

        plt.plot(x, y, 'o', x_smooth, y_smooth, '-')

        plt.savefig(os.path.join(args.plots_dir, names[bam] + ".pdf"))
        plt.close()


if __name__ == '__main__':
    ### Parse command-line arguments
    parser = ArgumentParser(
        description = 'correlations.py',
        usage       = 'python correlations.py <directory> file1, file2... '
    )

    ### Global options
    # positional arguments
    parser.add_argument(dest='plots_dir', type=str, help='Directory to save plots to.')
    parser.add_argument('bamfiles', nargs = '*', help = 'bamFiles')
    # optional arguments
    parser.add_argument('--duplicates', dest='duplicates', action='store_true')
    parser.add_argument('--window-width', dest='window_width', type=int, default=1000)
    parser.add_argument('--fragment-size', dest='fragment_size', type=int, default=1)
    parser.add_argument('--genome', dest='genome', type=str, default='hg19')

    args = parser.parse_args()
    # args = parser.parse_args([".", "data/human/chipmentation/mapped/merged/DNase_UWashington_K562_mergedReplicates.bam", "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam", "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_ChIP.bam"])
    
    main(args)

