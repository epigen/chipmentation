#!/usr/env python

from argparse import ArgumentParser
import os, re
from pybedtools import BedTool
import HTSeq
import numpy as np
import string
import random
import itertools
import cPickle as pickle

from matplotlib import pyplot as plt
from collections import Counter, OrderedDict

from scipy.stats.stats import pearsonr

def makeGenomeWindows(windowWidth, genome, step=None):
    """Generate windows genome-wide."""
    if step == None:
        step = windowWidth
    w = BedTool.window_maker(BedTool(), genome=genome, w=windowWidth, s=step)
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


def makeBedWindows(windowWidth, bedtool, step=None):
    """Generate windows within specified regions (e.g. from bed file)."""
    if step == None:
        step = windowWidth
    w = BedTool.window_maker(BedTool(), b=bedtool, w=windowWidth, s=step)
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


def bedToolsInterval2GenomicInterval(bedtool):
    """
    Given a pybedtools.BedTool object returns, dictionary of HTSeq.GenomicInterval objects.
    bedtool - a pybedtools.BedTool with intervals.
    """
    intervals = OrderedDict()
    for iv in bedtool:
        intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, iv.strand)
    return intervals


def distances(bam, intervals, fragmentsize, duplicates=True, orientation=True, permutate=False):
    """ Gets read coverage in bed regions, returns dict with region:count.
    bam=Bam object from HTSeq.BAM_Reader.
    intervals=dict with HTSeq.GenomicInterval objects as values.
    fragmentsize=integer.
    duplicates=boolean.
    orientation=boolean.
    permutate=boolean.
    """
    if orientation:
        distsPos = dict()
        distsNeg = dict()
    else:
        dists = dict()
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    
    for name, feature in itertools.islice(intervals.items(), 0, None, 10):
        if feature.chrom not in chroms:
            continue

        # Fetch all alignments in feature window
        alnsInWindow = bam[feature]

        if permutate:
            # randomize each alignment's position in window
            alns = list()
            for aln in alnsInWindow:
                aln.iv.start = random.randrange(feature.start, feature.end)
                aln.iv.length = fragmentsize
                alns.append(aln)
            alnsInWindow = alns

        # Fetch alignments in feature window
        for aln1, aln2 in itertools.combinations(alnsInWindow, 2):
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


def coverage(bam, intervals, fragmentsize, orientation=False, duplicates=True, strand_specific=False):
    """
    Gets read coverage in bed regions.
    Returns dict of regionName:numpy.array if strand_specific=False, A dict of "+" and "-" keys with regionName:numpy.array.
    bam=HTSeq.BAM_Reader object. Must be sorted and indexed with .bai file!
    intervals=dict with HTSeq.GenomicInterval objects as values.
    fragmentsize=integer.
    stranded=boolean.
    duplicates=boolean.
    """
    # Loop through TSSs, get coverage, append to dict
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX']
    cov = OrderedDict()
    n = len(intervals)
    i = 0
    for name, feature in intervals.iteritems():
        if i % 1000 == 0:
            print(n - i)
        # Initialize empty array for this feature
        if not strand_specific:
            profile = np.zeros(feature.length, dtype=np.float64)
        else:
            profile = np.zeros((2, feature.length), dtype=np.float64)

        # Check if feature is in bam index 
        if feature.chrom not in chroms or feature.chrom == "chrM":
            i+=1
            continue

        # Fetch alignments in feature window
        for aln in bam[feature]:
            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue
            aln.iv.length = fragmentsize # adjust to size

            # get position in relative to window
            if orientation:
                if feature.strand == "+" or feature.strand == ".":
                    start_in_window = aln.iv.start - feature.start - 1
                    end_in_window   = aln.iv.end   - feature.start - 1
                else:
                    start_in_window = feature.length - abs(feature.start - aln.iv.end) - 1
                    end_in_window   = feature.length - abs(feature.start - aln.iv.start) - 1
            else:
                start_in_window = aln.iv.start - feature.start - 1
                end_in_window   = aln.iv.end   - feature.start - 1
            
            # check fragment is within window; this is because of fragmentsize adjustment
            if start_in_window <= 0 or end_in_window > feature.length:
                continue

            # add +1 to all positions overlapped by read within window
            if not strand_specific:
                profile[start_in_window : end_in_window] += 1
            else:
                if aln.iv.strand == "+":
                    profile[0][start_in_window : end_in_window] += 1
                else:
                    profile[1][start_in_window : end_in_window] += 1
            
        # append feature profile to dict
        cov[name] = profile
        i+=1
    return cov


def extractPattern(dists):
    """
    Extracts most abundant frequency from raw data.
    Returns mirrored pattern.
    dists=dict of distances:counts.
    """
    #dists = pickle.load(open(os.path.join(args.results_dir, "DNase_UWashington_K562_mergedReplicates.counts.pickle"), "r"))
    #dists = pickle.load(open(os.path.join(args.results_dir, "H3K4me3_K562_500k_CM.counts.pickle"), "r"))

    # restrict to signal
    x = range(30, 130)
    y = [dists[i] for i in x]
    # fit linear regression
    p1 = np.poly1d(np.polyfit(x, y, 1))
    m , b = p1.coeffs
    #plt.plot(y, 'o', p1(x), "-")

    # measure distance to regression in that point
    distsReg = [y[i] - (m*i + b) for i in range(len(x))]

    # subtract value of minimum to all points
    distsReg -= min(distsReg)
    #plt.plot(distsReg)

    ### fourier transform data
    time = x
    signal = distsReg

    # get frequencies from decomposed fft
    W =  np.fft.fftfreq(signal.size, d=time[1]-time[0])
    f_signal = np.fft.fft(signal)

    # plot all frequencies, ask what is the amplitude of the signals
    freqs = list()
    for i in np.arange(0.01, len(x)/10., 0.01):
        cut_f_signal = f_signal.copy()
        cut_f_signal[((W != i) & (W != -i))] = 0
        cut_signal = np.fft.ifft(cut_f_signal)
        plt.plot(time, cut_signal)
        freqs.append(np.abs(cut_f_signal).max())

    # get frequency of signal with highest amplitude
    top = np.arange(0.01, len(x)/10., 0.01)[freqs.index(max(freqs))]

    # signal is now in Hz
    cut_f_signal = f_signal.copy()

    # plt.plot(W, abs(cut_f_signal), 'o')
    
    # select frequency of 1/10bp = 0.1Hz
    cut_f_signal[(W < top)] = 0
    cut_f_signal[(W > top)] = 0

    # plt.plot(W, abs(cut_f_signal), 'o')

    # inverse fourier to get filtered frequency
    cut_signal = np.fft.ifft(cut_f_signal)

    # plt.plot(W, abs(cut_signal), '-')

    # Extract pattern
    extracted = cut_signal.real

    return extracted

    ### Mirror pattern
    # add one more value to join phases
    #extractedPlus = np.append(extracted, extracted[10])
    #mirrored = np.concatenate((extracted, extracted[::-1]))
    #plt.plot(mirrored)

    # Set minimum to zero and scale down x100
    #mirrored += abs(min(mirrored))
    #mirrored *= 0.01
    #plt.plot(mirrored)
    #return mirrored


def correlatePatternProfile(pattern, profile, step=1):
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
        bp = 0
        while bp + len(pattern) <= len(profile):
            R.append(pearsonr(pattern, profile[bp:bp+len(pattern)])[0])
            i += 1
        return np.array(R)


def exportWigFile():
    """.
    """
    pass


def main(args):
    args.plots_dir = os.path.abspath(args.plots_dir)

    ### Get sample names
    names = list()
    for bamfile in args.bamfiles:
        names.append(re.sub("\.bam", "", os.path.basename(bamfile)))

    ### Loop through all signals, compute distances, plot
    # Get genome-wide windows
    print("Making %ibp windows genome-wide" % args.window_width)
    #windows = makeGenomeWindows(args.window_width, args.genome)
    windows = makeGenomeWindows(args.window_width, {'chr1': (0, 249250621)}, step=args.window_step)

    signals = dict()
    permutatedSignals = dict()
    for bam in xrange(len(args.bamfiles)):
        print("Sample " + names[bam])
        
        # Load bam
        bamfile = HTSeq.BAM_Reader(os.path.abspath(args.bamfiles[bam]))

        # Get dict of distances between reads genome-wide
        dists = distances(bamfile, windows, args.fragment_size, args.duplicates, orientation=False)
        pickle.dump(dists, open(os.path.join(args.results_dir, names[bam] + ".counts.pickle"), "wb"), protocol = pickle.HIGHEST_PROTOCOL)
        signals[bam] = dists
        #distsPos, distsNeg = distances(bamfile, windows, args.fragment_size, args.duplicates, orientation=True)
        #pickle.dump((distsPos, distsNeg), open(os.path.join(args.results_dir, names[bam] + ".countsStranded.pickle"), "wb"), protocol = pickle.HIGHEST_PROTOCOL)
        #signals[bam] = (distsPos, distsNeg)

        # Get dict of distances between permutated reads genome-wide
        permutedDists = distances(bamfile, windows, args.fragment_size, args.duplicates, orientation=False, permutate=True)
        pickle.dump(permutedDists, open(os.path.join(args.results_dir, names[bam] + ".countsPermuted.pickle"), "wb"), protocol = pickle.HIGHEST_PROTOCOL)
        permutatedSignals[bam] = dists


    ### For each signal extract most abundant periodic signal, correlate it with read coverage on each strand
    # generate windows genome-wide (or under 3K4 peaks)
    for bam in xrange(len(args.bamfiles)):
        distsPos, distsNeg = signals[bam]

        # sum strands count
        inp = [dict(x) for x in (distsPos, distsNeg)]
        count = Counter()
        for y in inp:
            count += Counter(y)
        dists = dict(count)
        #dists = pickle.load(open(os.path.join(args.results_dir, "H3K4me3_K562_500k_CM.counts.pickle"), "r"))

        # extract most abundant periodic pattern from signal
        pattern = extractPattern(dists)
        pickle.dump(pattern, open(os.path.join(args.results_dir, names[bam] + ".pattern.pickle"), "wb"), protocol = pickle.HIGHEST_PROTOCOL)

        permutedPattern = extractPattern(permutedDists)
        pickle.dump(permutedPattern, open(os.path.join(args.results_dir, names[bam] + ".patternPermuted.pickle"), "wb"), protocol = pickle.HIGHEST_PROTOCOL)

        # calculate read coverage in H3K4me3 peaks
        bamfile = HTSeq.BAM_Reader(os.path.abspath(args.bamfiles[bam]))
        peaks = BedTool("data/human/chipmentation/peaks/{sample}_peaks/{sample}_peaks.narrowPeak".format(sample=names[bam]))
        peaks = bedToolsInterval2GenomicInterval(peaks)
        #peaks = {name : peak for name, peak in peaks.items() if peak.length >= 1000}         # filter out peaks smaller than 1kb
        cov = coverage(bamfile, peaks, args.fragment_size)
        pickle.dump(cov, open(os.path.join(args.results_dir, names[bam] + ".peakCoverage.pickle"), "wb"), protocol = pickle.HIGHEST_PROTOCOL)
        #cov = pickle.load(open(os.path.join(args.results_dir, names[bam] + ".peakCoverage.pickle"), "r"))

        #cov = {name : reads for name, reads in cov.items() if len(reads) >= 1000}         # filter out peaks smaller than 1kb

        # Correlate coverage and signal pattern
        R = OrderedDict()
        #for peak, coverage in cov.iteritems():
        for peak, reads in cov.items()[:20]:
            R[peak] = correlatePatternProfile(pattern, reads)

        # plot pattern for 20 peaks
        for i in xrange(len(R)):
            plt.subplot(2,10,i)
            plt.plot(R.values()[i])

        # concatenate all peaks
        [i for l in R.values() for i in l]


        # export wig file with correlations
        for name, _ in R.items():
            exportWigFile(peaks[name], R[name], pattern)        


        # with paralelization
        import multiprocessing as mp
        pool = mp.Pool()
        R = [pool.apply(correlatePatternProfile, args=(pattern, reads)) for reads in cov.values()[:50]]

        for i in xrange(len(R)):
            plt.plot(R[i])

        # Get regions in extreme quantiles of correlation 

        # Check for enrichment in ...

        # mnase signal
        # nucleosomes
        # clusters of CM signal around TSSs



if __name__ == '__main__':
    ### Parse command-line arguments
    parser = ArgumentParser(
        description = 'read_distances.py',
        usage       = 'python read_distances.py <directory> file1, file2... '
    )

    ### Global options
    # positional arguments
    parser.add_argument(dest='results_dir', type=str, help='Directory to save data to.')
    parser.add_argument(dest='plots_dir', type=str, help='Directory to save plots to.')
    parser.add_argument('bamfiles', nargs = '*', help = 'bamFiles')
    # optional arguments
    parser.add_argument('--duplicates', dest='duplicates', action='store_true')
    parser.add_argument('--window-width', dest='window_width', type=int, default=1000)
    parser.add_argument('--window-step', dest='window_step', type=int, default=900)
    parser.add_argument('--fragment-size', dest='fragment_size', type=int, default=1)
    parser.add_argument('--genome', dest='genome', type=str, default='hg19')

    args = parser.parse_args()

    args = parser.parse_args(
        ["projects/chipmentation/results",
        "projects/chipmentation/results/plots",
        "data/human/chipmentation/mapped/merged/DNase_UWashington_K562_mergedReplicates.bam",
        "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam",
        "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_ChIP.bam",
        "data/human/chipmentation/mapped/merged/PU1_K562_10mio_CM.bam"
        ]
    )

    main(args)
