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

from matplotlib import pyplot as plt

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


def distances(bam, intervals, fragmentsize, duplicates=True, orientation=True):
    """ Gets read coverage in bed regions, returns dict with region:count.
    bam - Bam object from HTSeq.BAM_Reader.
    intervals - dict with HTSeq.GenomicInterval objects as values.
    fragmentsize - integer.
    duplicates - boolean.
    """
    if orientation:
        distsPos = dict()
        distsNeg = dict()
    else:
        dists = dict()
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    for name, feature in itertools.islice(intervals.items(), 0, None, 10000):
        if feature.chrom not in chroms:
            continue
        # Fetch alignments in feature window
        for aln1, aln2 in itertools.combinations(bam[feature], 2):
            # check if duplicate
            if not duplicates and (aln1.pcr_or_optical_duplicate or aln2.pcr_or_optical_duplicate):
                continue

            # check if in same strand
            if not orientation and aln1.iv.strand != aln2.iv.strand:
                continue
            # adjust fragment to size
            aln1.iv.length = fragmentsize
            aln2.iv.length = fragmentsize

            # get position relative
            dist = abs(aln1.iv.start_d - aln2.iv.start_d)
            # add +1 to dict
            if orientation:
                if aln1.iv.strand == "+":
                    if dist not in distsPos.keys():
                        distsPos[dist] = 1
                    else:
                        distsPos[dist] += 1

                if aln1.iv.strand == "-":
                    if dist not in distsNeg.keys():
                        distsNeg[dist] = 1
                    else:
                        distsNeg[dist] += 1
            else:
                if dist not in dists.keys():
                    dists[dist] = 1
                else:
                    dists[dist] += 1
    if orientation:
        return (distsPos, distsNeg)
    else:
        return dists


def main(args):
    args.plots_dir = os.path.abspath(args.plots_dir)

    # Get sample names
    names = list()
    for bam in args.bamfiles:
        names.append(re.sub("\.bam", "", os.path.basename(bam)))

        # Get genome-wide windows
        print("Making %ibp windows genome-wide" % args.window_width)
        #windows = makeWindows(args.window_width, args.genome)
        windows = makeWindows(args.window_width, {'chr1': (0, 249250621)})

        # Loop through all signals, compute distances, plot
        for bam in xrange(len(args.bamfiles)):
            print("Sample " + names[bam])
            # Load bam
            bamfile = HTSeq.BAM_Reader(os.path.abspath(args.bamfiles[bam]))
            # Get dataframe of bam coverage in bed regions, append to dict

            #dists = distances(bamfile, windows, args.fragment_size, args.duplicates, orientation=False)
            #pickle.dump(dists, open(names[bam] + ".counts.pickle", "wb"), protocol = pickle.HIGHEST_PROTOCOL)
            distsPos, distsNeg = distances(bamfile, windows, args.fragment_size, args.duplicates, orientation=True)
            pickle.dump((distsPos, distsNeg), open(names[bam] + ".countsStranded.pickle", "wb"), protocol = pickle.HIGHEST_PROTOCOL)

            x = range(1, 200)
            y = [dists[i] for i in x]
            plt.scatter(x, y)
            p20 = np.poly1d(np.polyfit(x, y, 20))
            p50 = np.poly1d(np.polyfit(x, y, 50))
            p100 = np.poly1d(np.polyfit(x, y, 100))

            plt.plot(x, p20(x), '-', x, p50(x), '--', x, p100(x), '.')
            plt.savefig(os.path.join(args.plots_dir, names[bam] + ".fit.pdf"))
            plt.close()

            x = np.array(range(1, 200))
            y = np.array([dists[i] for i in x])

            from scipy.interpolate import spline
            x_smooth = np.linspace(x.min(), x.max(), 200)
            y_smooth = spline(x, y, x_smooth)

            plt.plot(x, y, 'o', x_smooth, y_smooth, '-')

            plt.savefig(os.path.join(args.plots_dir, names[bam] + ".pdf"))
            plt.close()




# Frequency analysis
from spectrum import Periodogram
from scipy.interpolate import spline

### From raw data
#dists = pickle.load(open("DNase_UWashington_K562_mergedReplicates.counts.pickle", "r"))
dists = pickle.load(open("H3K4me3_K562_500k_CM.counts.pickle", "r"))

# sum strands count
from collections import Counter
inp = [dict(x) for x in (distsPos, distsNeg)]
count = Counter()
for y in inp:
  count += Counter(y)
dists = dict(count)

# fit linear regression
x = range(30, 130)
y = [dists[i] for i in x]
p1 = np.poly1d(np.polyfit(x, y, 1))
plt.plot(y, 'o', p1(x), "-")

# measure distance to regression in that point
m , b = p1.coeffs
distsReg = [y[i] - (m*i + b) for i in range(len(x))]

# subtract value of minimum to all points
distsReg -= min(distsReg)
#plt.plot(distsReg)

### Spectral analysis
p = Periodogram(distsReg, sampling = len(distsReg))
p.run()
p.plot(marker='o')


# fft
time = x
signal = distsReg

# get frequencies from decomposed fft
W =  np.fft.fftfreq(signal.size, d=time[1]-time[0])
f_signal = np.fft.fft(signal)

# signal is now in Hz    
cut_f_signal = f_signal.copy()
# filter noisy frequencies
cut_f_signal[(W < 0.1)] = 0
cut_f_signal[(W > 0.15)] = 0

# inverse fourier to get filtered frequency
cut_signal = np.fft.ifft(cut_f_signal)

# plot transformations
plt.subplot(221)
plt.plot(time,signal)
plt.subplot(222)
plt.plot(W, f_signal)
plt.subplot(223)
plt.plot(W, cut_f_signal)
plt.subplot(224)
plt.plot(time, cut_signal)
plt.show()


# Mirror pattern
plt.plot(np.concatenate((cut_signal, cut_signal[::-1])))


### With fitted data
# get best fit
x = range(60, 160)
y = [dists[i] for i in x]
p1 = np.poly1d(np.polyfit(x, y, 1))
x_smooth = np.linspace(min(x), max(x), len(x))
y_smooth = spline(x, y, x_smooth)

# fit
p100 = np.poly1d(np.polyfit(x, y, 100))
plt.plot(y, 'o', p1(x), "--", p100(x), "-", y_smooth, ".")

# fit linear regression
yPoli = p100(x)
p1 = np.poly1d(np.polyfit(x, yPoli, 1))
plt.plot(yPoli, 'o', p1(x), "-")

# measure distance to regression in that point
m , b = p1.coeffs
distsReg = [yPoli[i] - (m*i + b) for i in range(len(x))]

# subtract value of minimum to all points
distsReg -= min(distsReg)
#plt.plot(distsReg)

### Spectral analysis
# Periodogram
p = Periodogram(distsReg, sampling=len(distsReg))
p.run()
p.plot(marker='o')

# fft
yFFT = np.fft.fft(distsReg)
freqs = np.fft.fftfreq(len(distsReg))
plt.plot(2.0/len(distsReg) * np.abs(yFFT[0:len(distsReg)/2]), 'o')

# Mirror pattern
plt.plot(np.concatenate((distsReg, distsReg[::-1])))



####### With stranded data
distsPos, distsNeg = pickle.load(open("H3K4me3_K562_500k_CM.countsStranded.pickle", "r"))

# raw
for strand in [distsPos, distsNeg]:
    # extract region
    x = range(60, 200)
    y = [strand[i] for i in x]
    # fit linear
    p1 = np.poly1d(np.polyfit(x, y, 1))
    plt.plot(y, 'o', p1(x), "-")

    # measure distance to regression in that point
    m , b = p1.coeffs
    distsReg = [yPoli[i] - (m*i + b) for i in range(len(x))]

    # subtract value of minimum to all points
    distsReg -= min(distsReg)
    plt.plot(distsReg)

    ### Spectral analysis
    p = Periodogram(distsReg, sampling = len(distsReg))
    p.run()
    p.plot()


# fitted data
for strand in [distsPos, distsNeg]:
    # get best fit
    x = range(55, 155)
    y = [strand[i] for i in x]
    p1 = np.poly1d(np.polyfit(x, y, 1))
    x_smooth = np.linspace(min(x), max(x), len(x))
    y_smooth = spline(x, y, x_smooth)

    # fit
    p100 = np.poly1d(np.polyfit(x, y, 100))
    plt.plot(y, 'o', p1(x), "--", p100(x), "-", y_smooth, ".")

    # fit linear regression
    yPoli = p100(x)
    p1 = np.poly1d(np.polyfit(x, yPoli, 1))
    plt.plot(yPoli, 'o', p1(x), "-")

    # measure distance to regression in that point
    m , b = p1.coeffs
    distsReg = [yPoli[i] - (m*i + b) for i in range(len(x))]
    
    # subtract value of minimum to all points
    distsReg -= min(distsReg)
    plt.plot(distsReg)

    ### Spectral analysis
    p = Periodogram(distsReg, sampling=len(distsReg))
    p.run()
    p.plot(marker='o')






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
    
    args = parser.parse_args(
        ["projects/chipmentation/results/plots",
        "data/human/chipmentation/mapped/merged/DNase_UWashington_K562_mergedReplicates.bam",
        "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam",
        "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_ChIP.bam",
        "data/human/chipmentation/mapped/merged/PU1_K562_10mio_CM.bam"
        ]
    )

    main(args)
