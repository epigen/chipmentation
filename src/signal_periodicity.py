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
import seaborn as sns
from collections import Counter, OrderedDict

from scipy.stats.stats import pearsonr
import statsmodels.nonparametric.kde as kde
from scipy import signal

import yahmm

class WinterHMM(object):
    """
    Class to model the Hidden Markov model described in the Winter 2013 paper.
    """
    def __init__(self):
        super(WinterHMM, self).__init__()
        background = yahmm.DiscreteDistribution({
            1: 0.33,
            0: 0.33,
            -1 : 0.33
        })
        pos = yahmm.DiscreteDistribution({
            1: 0.7,
            0: 0.2,
            -1 : 0.1
        })
        neg = yahmm.DiscreteDistribution({
            1: 0.1,
            0: 0.1,
            -1 : 0.8
        })
        zero = yahmm.DiscreteDistribution({
            1: 0,
            0: 1,
            -1 : 0
        })

        # Create states
        b = yahmm.State(background, name="B")
        in1 = yahmm.State(zero, name="in1")
        in2 = yahmm.State(zero, name="in2")
        in3 = yahmm.State(zero, name="in3")
        in4 = yahmm.State(zero, name="in4")
        minus = yahmm.State(neg, name="-")
        in5 = yahmm.State(zero, name="in5")
        in6 = yahmm.State(zero, name="in6")
        in7 = yahmm.State(zero, name="in7")
        in8 = yahmm.State(zero, name="in8")
        in9 = yahmm.State(zero, name="in9")
        in10 = yahmm.State(zero, name="in10")
        in11 = yahmm.State(zero, name="in11")
        plus = yahmm.State(pos, name="+")

        self.model = yahmm.Model(name="winter-2013")

        self.model.add_state(b)
        self.model.add_state(in1)
        self.model.add_state(in2)
        self.model.add_state(in3)
        self.model.add_state(in4)
        self.model.add_state(minus)
        self.model.add_state(in5)
        self.model.add_state(in6)
        self.model.add_state(in7)
        self.model.add_state(in8)
        self.model.add_state(in9)
        self.model.add_state(in10)
        self.model.add_state(in11)
        self.model.add_state(plus)

        self.model.add_transition(self.model.start, b, 1)

        self.model.add_transition(b, b, 0.5)
        self.model.add_transition(b, minus, 0.25)
        self.model.add_transition(b, plus, 0.25)

        self.model.add_transition(in1, in2, 1)
        self.model.add_transition(in2, in3, 1)
        self.model.add_transition(in3, in4, 1)
        self.model.add_transition(in4, minus, 1)

        self.model.add_transition(minus, b, 0.25)
        self.model.add_transition(minus, in5, 0.25)

        self.model.add_transition(in5, in6, 1)
        self.model.add_transition(in6, in7, 1)
        self.model.add_transition(in7, in8, 1)
        self.model.add_transition(in8, in9, 1)
        self.model.add_transition(in9, in10, 1)
        self.model.add_transition(in10, in11, 1)
        self.model.add_transition(in11, plus, 1)

        self.model.add_transition(plus, b, 0.50)
        self.model.add_transition(plus, in1, 0.50)

        self.model.bake()


    def draw(self):
        """
        Draws the Markov chain of the model.
        """
        return self.model.draw(node_size=400, labels={state.name : str(state.name) for state in self.model.states}, font_size=20)


    def train(self, sequences, **kwargs):
        """
        Train model with given data.
        """
        self.model.train(sequences)


    def predict(self, observations):
        """
        Predict hidden states from observations.
        """
        logp, chain = self.model.maximum_a_posteriori(observations)
        return [state.name for num, state in chain][1:-1] # strip begin and end states
        #model.forward(observations)
        #model.backward(observations)
        #model.forward_backward(observations)
        #-model.log_probability(observations)
        #logp, chain = self.model.maximum_a_posteriori(observations)
        #return (logp, [(num, name.name) for num, name in chain])


def makeGenomeWindows(windowWidth, genome, step=None):
    """
    Generate windows genome-wide for a given genome with width=windowWidth and
    return dictionary of HTSeq.GenomicInterval objects.

    windowWidth=int.
    genome=str.
    """
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
    """
    Generate windows with width=windowWidth given a pybedtools.BedTool object and
    return dictionary of HTSeq.GenomicInterval objects.

    windowWidth=int.
    bedtool=pybedtools.BedTool object.
    """
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


def bedToolsInterval2GenomicInterval(bedtool , strand=True):
    """
    Given a pybedtools.BedTool object, return dictionary of HTSeq.GenomicInterval objects.

    bedtool=pybedtools.BedTool object.
    """
    intervals = OrderedDict()
    if strand:
        for iv in bedtool:
            intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, iv.strand)
    else:
        for iv in bedtool:
            intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end)

    return intervals


def bedToolsInterval2GenomicInterval_noFormat(bedtool):
    """
    Given a pybedtools.BedTool object, return dictionary of HTSeq.GenomicInterval objects.

    bedtool=pybedtools.BedTool object.
    """
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


def distances(bam, intervals, fragmentsize, duplicates=True, orientation=True, permutate=False):
    """
    Gets pairwise distance between reads in intervals. Returns dict with distance:count.
    If permutate=True, it will randomize the reads in each interval along it.

    bam=HTSeq.BAM_Reader object.
    intervals=dict - HTSeq.GenomicInterval objects as values.
    fragmentsize=int.
    duplicates=bool.
    orientation=bool.
    permutate=bool.
    """
    if orientation:
        distsPos = dict()
        distsNeg = dict()
    else:
        dists = dict()
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    
    #for name, feature in intervals.iteritems():
    for name, feature in itertools.islice(intervals.items(), 0, None, 10):
        if feature.chrom not in chroms:
            continue
        # Fetch all alignments in feature window
        alnsInWindow = bam[feature]
        if permutate:
            # randomize each alignment's position in window
            alns = list()
            for aln in alnsInWindow:
                aln.iv.start_d = random.randrange(feature.start, feature.end)
                alns.append(aln)
            alnsInWindow = alns

        # Measure distance between reads in window, pairwisely
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
        #plt.plot(dists.keys(), dists.values(), 'o')
        return dists


def coverageInWindows(bam, intervals, fragmentsize, orientation=False, duplicates=True, strand_specific=False, permutate=False):
    """
    Gets read coverage in intervals. Returns dict of regionName:numpy.array if strand_specific=False,
    a dict of "+" and "-" keys with regionName:numpy.array.

    bam=HTSeq.BAM_Reader object - Must be sorted and indexed with .bai file!
    intervals=dict - HTSeq.GenomicInterval objects as values.
    fragmentsize=int.
    stranded=bool.
    duplicates=bool.
    """
    # Loop through TSSs, get coverage, append to dict
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX']
    coverage = OrderedDict()
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
        coverage[name] = profile
        i+=1
    return coverage


def extractPattern(dists, distRange, filePrefix):
    """
    Performs the fast Fourier transform (fft) on a dict of x:y values in several frequencies, and
    returns the signal of the most abundant (highest signal amplitude) frequency through inverse fft.
    Produces plots of the process.

    dists=dict - distances:counts, both int.
    distRange=int - indexes of subset of dists to extract pattern from.
    filePrefix=str - prefix of files to save plots to (should include paths).
    """
    plt.close()
    plt.figure(0)
    plt.subplot(211)
    x = dists.keys()
    y = [dists[i] for i in x]
    p1 = np.poly1d(np.polyfit(x, y, 1)) # fit linear regression
    m , b = p1.coeffs
    plt.plot(
        x, y, 'o',
        x, p1(x), "-"
    )    
    
    # restrict to signal
    x = distRange
    y = [dists[i] for i in x]
    p1 = np.poly1d(np.polyfit(x, y, 1)) # fit linear regression
    m , b = p1.coeffs
    plt.subplot(212)
    plt.plot(
        x, y, 'o',
        x, p1(x), "-"
    )
    plt.savefig(filePrefix + ".read_distances.pdf", bbox_inches='tight')
    plt.close()

    # measure distance to regression in that point
    distsReg = [y[i] - (m*i + b) for i in range(len(x))]

    # subtract value of minimum to all points
    distsReg -= min(distsReg)

    ### fourier transform data
    time = x
    signal = distsReg

    # get frequencies from decomposed fft
    W =  np.fft.fftfreq(signal.size, d=time[1]-time[0])
    f_signal = np.fft.fft(signal)

    # plot all frequencies, ask what is the amplitude of the signals
    freqs = dict()
    plt.figure(1)
    for i in np.arange(0, len(x)/10., 0.01):
        if i == 0:
            continue
        else:
            cut_f_signal = f_signal.copy()
            cut_f_signal[((W < i) | (W > i))] = 0
            cut_signal = np.fft.ifft(cut_f_signal)
            plt.plot(time, cut_signal)
            if 0.05 <= i <= 0.4:
                freqs[i] = np.abs(cut_f_signal).max()
    plt.savefig(filePrefix + ".fft_frequencies.pdf", bbox_inches='tight')
    plt.close()

    # get frequency of signal with highest amplitude
    top = max(freqs, key=freqs.get)

    # signal is now in Hz
    cut_f_signal = f_signal.copy()

    # select frequency of top/10bp
    cut_f_signal[((W < top) | (W > top))] = 0

    # inverse fourier to get filtered frequency
    cut_signal = np.fft.ifft(cut_f_signal)

    plt.figure(2)
    plt.subplot(221)
    plt.plot(time, signal, '-')
    plt.subplot(222)
    plt.plot(W, abs(f_signal), 'o')
    plt.subplot(223)
    plt.plot(W, abs(cut_f_signal), 'o')
    plt.subplot(224)
    plt.plot(time, cut_signal, '-')
    plt.savefig(filePrefix + ".fft_filter-{0}bp_ifft.pdf".format(int(1/top)), bbox_inches='tight')
    plt.close()

    return cut_signal.real


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
        i = 0
        while i + len(pattern) <= len(profile):
            R.append(pearsonr(pattern, profile[i:i+len(pattern)])[0])
            i += 1
        return np.array(R)


def exportWigFile(intervals, profiles, offset, filename, trackname):
    """
    Exports a wig file track with scores contained in profile, .

    intervals=iterable with HTSeq.GenomicInterval objects.
    profiles=iterable with scores to write.
    offset=int.
    filename=str.
    trackname=str.
    """
    with open(filename, 'w') as handle:
        track = 'track type=wiggle_0 name="{0}" description="{0}" visibility=dense autoScale=off\n'.format(trackname)
        handle.write(track)
        for i in xrange(len(profiles)):
            header = "fixedStep  chrom={0}  start={1}  step=1\n".format(intervals[i].chrom, intervals[i].start + offset)
            handle.write(header)    
            for j in xrange(len(profiles[i])):
                handle.write(str(abs(profiles[i][j])) + "\n")


def fitKDE(X, bandwidth):
    """
    Fit gaussian kernel density estimator to array, return

    X=numerical iterable.
    bandwidth=kernel bandwidth - float.
    """
    return kde.kdensity(range(len(x)), bw=bandwidth, weights=x)[0]


def fitRegression(X, degree):
    """
    Fit regression with degree.

    X=numerical iterable.
    degree=int.
    """
    coeffs = np.poly1d(np.polyfit(range(len(X)), X, degree))
    return np.polyval(coeffs, range(len(X)))


def findPeaks(X, peakWidth):
    """
    Find peaks of local maxima with width peakWidth and returns them with that width.

    X=dict with iterable with numerical values as dict values.
    peakWidth=integer.
    """
    peaks = signal.find_peaks_cwt(X, widths=np.array([peakWidth]))
    return {peak : (peak - peakWidth/2, peak + peakWidth/2) for peak in peaks}


def binarize(X):
    """
    Convert a numerical iterable (X) in a np.array with binary values by finding local maximas (1)
    and setting the rest to zeros (0). Local maximas are required to be higher than 0.05, and have negative
    values between them.

    X=iterable with numerical values.
    """
    maximas = np.r_[True, X[1:] > X[:-1]] & np.r_[X[:-1] > X[1:], True]
    binary = list()
    prev_max = 0
    for i in xrange(len(X)):
        if maximas[i] == True and X[i] > 0.05 and any(n < 0 for n in X[prev_max : i]): # watch out for namespace pollution with np.any
            # if from i to next_max there is no neg, then is max
            for j in xrange(i, len(X)):
                if maximas[j] == True and X[j] > 0.05:
                    next_max = j
                    break

            if not any(n > 0 for n in X[i : next_max]):
                binary.append(1)
                prev_max = i
            else:
                binary.append(0)
        else:
            binary.append(0)
    return np.array(binary)


def getConsensus(seq1, seq2):
    """
    Given two binary (1, 0) sequences, return one sequence with:
    1 if only one of them is 1, 0 otherwise.

    seq1,seq2=iterables with binary numerical values (0 or 1).
    """
    if not len(seq1) == len(seq2):
        return None
    else:
        binary = list()
        for i in xrange(len(seq1)):
            if seq1[i] == np.nan or seq1[i] == np.nan:
                binary.append(0)
            else:
                if seq1[i] == 1:
                    if seq2[i] == 1:
                        binary.append(0)
                    else:
                        binary.append(1)
                else:
                    if seq2[i] == 1:
                        binary.append(-1)
                    else:
                        binary.append(0)
        return np.array(binary)


def exportBedFile(intervals, peaks, offset, filename, trackname, strand="."):
    """
    Exports a bed file track.

    intervals=iterable with HTSeq.GenomicInterval objects.
    peaks=dict with tuple of start,end positions.
    offset=int.
    filename=str.
    trackname=str.
    strand=str.
    """
    if len(intervals) != len(peaks):
        raise TypeError("Intervals and peaks sequences have different lengths.")

    with open(filename, 'w') as handle:
        track = 'track name="{0}" description="{0}" visibility=pack autoScale=off colorByStrand="255,0,0 0,0,255"\n'.format(trackname)
        handle.write(track)
        for name, peak in peaks.items():
            for center, tup in peak.items():
                ### TODO: check peak boundaries
                entry = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    intervals[name].chrom,
                    intervals[name].start + offset + tup[0],
                    intervals[name].start + offset + tup[1],
                    name,
                    1,      # score
                    strand  # strand
                )
                handle.write(entry)


def main(args):
    args.plots_dir = os.path.abspath(args.plots_dir)

    ### Get sample names
    names = [re.sub("\.bam", "", os.path.basename(name)) for name in args.bamfiles]

    ### Loop through all signals, compute distances, plot
    # Get genome-wide windows
    print("Making %ibp windows genome-wide" % args.window_width)
    windows = makeGenomeWindows(args.window_width, args.genome)
    #windows = makeGenomeWindows(args.window_width, {'chr1': (0, 249250621)}, step=args.window_step) # only chromosome 1

    for index in xrange(len(args.bamfiles)):
        print("Sample " + names[index])
        
        # Load bam
        bam = HTSeq.BAM_Reader(os.path.abspath(args.bamfiles[index]))

        ### Get dict of distances between reads genome-wide
        distsPos, distsNeg = distances(bam, windows, args.fragment_size, args.duplicates, orientation=True)
        pickle.dump((distsPos, distsNeg), open(os.path.join(args.results_dir, names[index] + ".countsStranded.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        ### Get dict of distances between permutated reads genome-wide
        permutedDistsPos, permutedDistsNeg = distances(bam, windows, args.fragment_size, args.duplicates, orientation=True, permutate=True)
        pickle.dump((permutedDistsPos, permutedDistsNeg), open(os.path.join(args.results_dir, names[index] + ".countsPermutedStranded.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        # calculate distances in non-dhs regions for DNase-data
        if index == 0:
            # get windows not overlapping with DHS
            dhs = BedTool(os.path.join("/home", "arendeiro", "wgEncodeOpenChromDnaseK562Pk.narrowPeak"))
            w = BedTool.window_maker(BedTool(), genome={'chr1': (0, 249250621)}, w=args.window_width, s=args.window_width)
            w = bedToolsInterval2GenomicInterval_noFormat(w.intersect(b=dhs,v=True))
            
            # compute distances
            distsPos, distsNeg = distances(bam, w, args.fragment_size, args.duplicates, orientation=True)
            pickle.dump((distsPos, distsNeg), open(os.path.join(args.results_dir, names[index] + ".countsStranded-nonDHS.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
            permutedDistsPos, permutedDistsNeg = distances(bam, w, args.fragment_size, args.duplicates, orientation=True, permutate=True)
            pickle.dump((permutedDistsPos, permutedDistsNeg), open(os.path.join(args.results_dir, names[index] + ".countsPermutedStranded-nonDHS.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

            # dhs = bedToolsInterval2GenomicInterval(dhs)
            # { name : window for name, window in windows.items() for _, dh in dhs.items() if not window.overlaps(dh) }

    ### For each signal extract most abundant periodic signal, correlate it with read coverage on each strand
    # generate windows genome-wide (or under 3K4 peaks)
    for index in xrange(len(args.bamfiles)):
        distsPos, distsNeg = pickle.load(open(os.path.join(args.results_dir, names[index] + ".countsStranded.pickle"), "r"))
        permutedDistsPos, permutedDistsNeg = pickle.load(open(os.path.join(args.results_dir, names[index] + ".countsPermutedStranded.pickle"), "r"))

        ### Extract most abundant periodic pattern from signal
        # for DNase, extract from a different window (70-150bp)
        if index == 0:
            patternPos = extractPattern(distsPos, range(70, 110), os.path.join(args.plots_dir, names[index] + "_posStrand"))
            patternNeg = extractPattern(distsNeg, range(70, 110), os.path.join(args.plots_dir, names[index] + "-_negStrand"))
            pattern = extractPattern(Counter(distsPos) + Counter(distsNeg), range(70, 110), os.path.join(args.plots_dir, names[index] + "_bothStrands"))

            permutedPatternPos = extractPattern(permutedDistsPos, range(70, 110), os.path.join(args.plots_dir, names[index] + "_posStrand_permuted"))
            permutedPatternNeg = extractPattern(permutedDistsNeg, range(70, 110), os.path.join(args.plots_dir, names[index] + "_negStrand_permuted"))
            permutedPattern = extractPattern(Counter(permutedDistsPos) + Counter(permutedDistsNeg), range(70, 110), os.path.join(args.plots_dir, names[index] + "_bothStrands_permuted"))
        else:
            patternPos = extractPattern(distsPos, range(60, 100), os.path.join(args.plots_dir, names[index] + "_posStrand"))
            patternNeg = extractPattern(distsNeg, range(60, 100), os.path.join(args.plots_dir, names[index] + "_negStrand"))
            pattern = extractPattern(Counter(distsPos) + Counter(distsNeg), range(60, 100), os.path.join(args.plots_dir, names[index] + "_bothStrands"))
        
            permutedPatternPos = extractPattern(permutedDistsPos, range(60, 100), os.path.join(args.plots_dir, names[index] + "_posStrand_permuted"))
            permutedPatternNeg = extractPattern(permutedDistsNeg, range(60, 100), os.path.join(args.plots_dir, names[index] + "_negStrand_permuted"))
            permutedPattern = extractPattern(Counter(permutedDistsPos) + Counter(permutedDistsNeg), range(60, 100), os.path.join(args.plots_dir, names[index] + "_bothStrands_permuted"))

        ### calculate read coverage in H3K4me3 peaks
        bam = HTSeq.BAM_Reader(os.path.abspath(args.bamfiles[index]))
        peaks = BedTool("data/human/chipmentation/peaks/{sample}_peaks/{sample}_peaks.narrowPeak".format(sample=names[1]))
        peaks = bedToolsInterval2GenomicInterval(peaks)
        # peaks = {name : peak for name, peak in peaks.items() if peak.length >= 1000}         # filter out peaks smaller than 1kb
        
        coverage = coverageInWindows(bam, peaks, args.fragment_size, strand_specific=True)
        pickle.dump(coverage, open(os.path.join(args.results_dir, names[index] + ".peakCoverageStranded.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
        
        coverage = pickle.load(open(os.path.join(args.results_dir, names[index] + ".peakCoverageStranded.pickle"), "r"))

        # coverage = {name : reads for name, reads in coverage.items() if len(reads) >= 1000} # Filter out peaks smaller than 1kb

        ### Correlate coverage and signal pattern
        pattern = patternNeg # choose one pattern (e.g. negative strand)
        correlationsPos = {peak : correlatePatternProfile(pattern, reads[0]) for peak, reads in coverage.items()}
        correlationsNeg = {peak : correlatePatternProfile(pattern, reads[1]) for peak, reads in coverage.items()}
        pickle.dump((correlationsPos, correlationsNeg), open(os.path.join(args.results_dir, names[index] + ".peakCorrelationStranded.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
        
        correlationsPos, correlationsNeg = pickle.load(open(os.path.join(args.results_dir, names[index] + ".peakCorrelationStranded.pickle"), "r"))
        
        ### Binarize data for HMM
        # binarize data strand-wise
        pos = dict()
        for name, cor in correlationsPos.items():
            b = binarize(cor)
            pos[name] = b

        neg = dict()
        for name, cor in correlationsNeg.items():
            b = binarize(cor)
            neg[name] = b

        # join strands in one sequence
        binary = {peak : getConsensus(pos[peak], neg[peak]) for peak in pos.keys()}

        # get zeroes for rest of the genome

        ### HMM
        model = WinterHMM()
        ## Train
        model.train(binary.values()[:10000]) # subset data for training
        # see new probabilities
        [(s.name, s.distribution) for s in model.model.states]
        # save model with trained parameters
        pickle.dump(model, open(os.path.join(args.results_dir, "hmm_trained_with_" + names[index] + ".pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
        ## Predict
        hmmOutput = {peak : model.predict(sequence) for peak, sequence in binary.items()}

        ### TODO:
        ## Get DARNS from hmmOutput
        ### Output bed/wig


        ## Export wig files with raw correlations 
        exportWigFile(
            [peaks[i] for i in correlationsPos.keys()],
            correlationsPos.values(),
            len(pattern)/2,
            os.path.join(args.results_dir, names[index] + ".peakCorrelationPos.wig"),
            names[index] + " raw absolute correlation - positive strand"
        )
        exportWigFile(
            [peaks[i] for i in correlationsNeg.keys()],
            correlationsNeg.values(),
            len(pattern)/2,
            os.path.join(args.results_dir, names[index] + ".peakCorrelationNeg.wig"),
            names[index] + " raw absolute correlation - negative strand"
        )
        ## Export wig files with correlations peaks

        ## Export wig files with DARNS

        ## Export bed files with DARNS


        ### Genome arithmetics with DARNS:
        ## Plot :
        # distances between DARNS borders
        # distances between DARNS centers
        # DARNS length



        # ### Find correlation peaks (broadish regions around region of high pattern-reads correlation)
        # # i = 0
        # corPeaksPos = dict()
        # for name, cor in correlationsPos.items():
        #     corPeak = findPeaks(abs(cor), 50)
        #     corPeaksPos[name] = corPeak
        #     # plt.subplot(5,2,i)
        #     # plt.plot(abs(cor))
        #     # plt.plot(corPeak.keys(), [0.2] * len(corPeak.keys()), 'o')
        #     # for center, peak in corPeak.items():
        #     #     plt.plot(range(peak[0], peak[1]), [0.19] * (peak[1] - peak[0]), '-')
        #     # i += 1

        # exportBedFile(
        #     {name : peaks[name] for name, peak in corPeaksPos.items()},
        #     corPeaksPos,
        #     len(pattern)/2,
        #     os.path.join(args.results_dir, names[index] + ".correlationPeaksPos.bed"),
        #     names[index] + " correlation peaks - positive strand"
        # )

        # corPeaksNeg = dict()
        # for name, cor in correlationsNeg.items():
        #     corPeak = findPeaks(abs(cor), 50)
        #     corPeaksNeg[name] = corPeak

        # exportBedFile(
        #     {name : peaks[name] for name, peak in corPeaksNeg.items()},
        #     corPeaksNeg,
        #     len(pattern)/2,
        #     os.path.join(args.results_dir, names[index] + ".correlationPeaksNeg.bed"),
        #     names[index] + " correlation peaks - negative strand"
        # )

        ##### ALTERNATIVE WAYS TO HANDLE CORRELATIONS
        # density = {name : fitKDE(values) for name, values in correlationsPos.items()}
        
        # X = abs(correlationsPos.values()[1])
        # for bandwidth in xrange(1, 10):
        #     density = kde.kdensity(range(len(X)), bw=bandwidth / 10., weights=X)[0]
        #     plt.plot(density * 100)

        # for degree in xrange(10, 50):
        #     coeffs = np.poly1d(np.polyfit(range(len(X)), X, 50))
        #     plt.plot(np.polyval(coeffs, range(len(X))) * 2)

        # # smooth with interpolate
        # from scipy.interpolate import interp1d
        # for degree in xrange(2, 15):
        #     xn_ax = interp1d(range(len(X)), X, kind=degree)
        #     plt.plot(xn_ax(X))

        # # find local minima
        # mins = np.r_[True, X[1:] < X[:-1]] & np.r_[X[:-1] < X[1:], True]

        # min_index = [np.where(X==X[x]) for x in mins]

        # [np.mean(x - 1, x + 1) for x in mins if x == True]


        # # find peaks
        # X = [abs(i) for l in correlationsPos.values() for i in l]
        # corPeaks = findPeaks(np.array(X), 50)

        # plt.plot(X)
        # plt.plot(corPeaks.keys(), [0.2] * len(corPeaks.keys()), 'o')
        # for center, peak in corPeaks.items():
        #     plt.plot(range(peak[0], peak[1]), [0.19] * (peak[1] - peak[0]), '-')

        
        # plot correlations for 20 peaks
        # for i in xrange(len(correlations)):
        #     plt.subplot(len(correlations)/10, len(correlations)/2,i)
        #     plt.plot(correlations.values()[i])
        
        # plot concatenation of correlations for some peaks
        # plt.plot([abs(i) for l in correlations.values()[:5] for i in l])



        # Winter2013 fig.4 - chr19 comparison
        # peaksWinter = {"winter" : HTSeq.GenomicInterval("chr19", 37016000, 37022000)}
        # covWinter = coverageInWindows(bam, peaksWinter, args.fragment_size)
        # RWinter = OrderedDict()
        # RWinter["winter"] = correlatePatternProfile(pattern, covWinter["winter"])
        # exportWigFile(
        #     [peaksWinter[i] for i in RWinter.keys()],
        #     RWinter.values(),
        #     len(pattern)/2,
        #     os.path.join(args.results_dir, names[index] + ".peakCorrelation.Winter.wig"),
        #     names[index]
        # )

        # with paralelization
        # import multiprocessing as mp
        # pool = mp.Pool()
        # correlations = [pool.apply(correlatePatternProfile, args=(pattern, reads)) for reads in coverage.values()[:50]] 
        # for i in xrange(len(correlations)):
        #    plt.plot(correlations[i])

        ### TODO:
        # Get regions in extreme quantiles of correlation 
        # Check for enrichment in ...
        #     mnase signal
        #     nucleosomes
        #     clusters of CM signal around TSSs



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
    main(args)
