#!/usr/env python

import os, re
from pybedtools import BedTool
import HTSeq
import numpy as np
import pandas as pd
import string

import rpy2.robjects as robj # for ggplot in R
import rpy2.robjects.pandas2ri # for R dataframe conversion

import itertools

def makeWindows(windowWidth):
    """Generate 1kb windows genome-wide."""
    w = BedTool.window_maker(BedTool(), genome="hg19", w=windowWidth)
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
    i = 0
    for name, feature in itertools.islice(intervals.items(), 0, None, 10000):
        if i % 1000 == 0:
            print(i)
        if feature.chrom not in chroms:
            i += 1
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
        i += 1
    return dists


def main(args):
    
    # Get sample names
    names = list()
    for bam in bamfiles:
        names.append(re.sub(os.path.basename(bam), "\.bam", ""))

    # Get genome-wide windows
    windows = makeWindows(args.windowWidth)

    # Loop through all signals, compute coverage in bed regions, append to dict
    rawSignals = dict()
    for bam in xrange(len(args.bamfiles)):
        # Load bam
        bamfile = HTSeq.BAM_Reader(args.bamfiles[bam])
        # Get dataframe of bam coverage in bed regions, append to dict
        dists = distances(bamfile, windows, args.fragment_size, args.duplicates)


x = range(50, 200)
y = [dists[i] for i in x]
plt.scatter(x, y)
p20 = np.poly1d(np.polyfit(x, y, 20))
p50 = np.poly1d(np.polyfit(x, y, 50))
p100 = np.poly1d(np.polyfit(x, y, 100))

plt.plot(x, p20(x), '-', x, p50(x), '--', x, p100(x), '.')

x = np.array(range(50, 200))
y = np.array([dists[i] for i in x])

from scipy.interpolate import spline
x_smooth = np.linspace(x.min(), x.max(), 200)
y_smooth = spline(x, y, x_smooth, 2)

plt.plot(x, y, 'o', x_smooth, y_smooth, '-')


if __name__ == '__main__':
    ### Parse command-line arguments
    parser = ArgumentParser(
        description = 'correlations.py',
        usage       = 'python correlations.py <directory> file1, file2... '
    )

    parser.add_argument('bamfiles', nargs = '*', help = 'bamFiles')

    ### Global options
    # positional arguments
    # optional arguments
    parser.add_argument(dest='plots_dir', type=str, help='Directory to save plots to.')
    parser.add_argument('--duplicates', dest='duplicates', type=bool, action='store_true')
    parser.add_argument('--window-width', dest='window_width', type=int, default=1000)
    parser.add_argument('--fragment-size', dest='fragment_size', type=int, default=1)
    parser.add_argument('--genome', dest='genome', type=str, default='hg19')
    parser.add_argument('bamfiles', nargs = '*', help = 'bamFiles')

    main(args)
