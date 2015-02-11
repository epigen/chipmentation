#!/usr/env python

"""
A more deep analysis of ChIPmentation data.

In brief, this does:
    * Count pairwise distances between reads from several samples (from several
    techniques - ChIP, ChIPmentation, DNase-seq, ATAC-seq) in several genomic
    locations (H3K4me3/H3K27me3 peaks, DHS regions, whole genome, etc...) in strand-
    -specific way and also for permuted reads.
    * Extract pattern from the read distances distribution.
    * Decompose pattern into signals with various frequencies (using FFT) and retrieve
    most abundant.
    * Calculate correlations between pattern and ChIPmentation read coverage along the genome.
    * Extract local maxima from correlation into a binary signal.
    * Feed to a HMM modeling a nucleossome, output predicted nucleossome-associated regions.

    * ... plus several plots along the way.

TODO:
    Repeat read counting in regions minus repeats and gaps - ongoing
    Repeat read counting in regions without duplicates

    Implement coverage and correlation at same time Task


"""


from argparse import ArgumentParser
from collections import Counter, OrderedDict
from divideAndSlurm import DivideAndSlurm, Task
from matplotlib import pyplot as plt
from pybedtools import BedTool
from scipy import signal
from scipy.stats.stats import pearsonr
import cPickle as pickle
import HTSeq
import itertools
import numpy as np
import os
import sys
import random
import re
import seaborn as sns
import statsmodels.nonparametric.kde as kde
import string
import subprocess
import textwrap
import time
import yahmm

np.set_printoptions(linewidth=200)


def main(args):
    # Parse command-line arguments
    parser = ArgumentParser(
        description='read_distances.py',
        usage='python read_distances.py <directory> file1, file2... '
    )

    # Global options
    # positional arguments
    parser.add_argument(dest='data_dir', type=str, help='Directory to save data to.')
    parser.add_argument(dest='results_dir', type=str, help='Directory to save data to.')
    parser.add_argument(dest='plots_dir', type=str, help='Directory to save plots to.')
    parser.add_argument('bam_files', nargs='*', help='Bam files')
    # optional arguments
    parser.add_argument('--duplicates', dest='duplicates', action='store_true')
    parser.add_argument('--window-width', dest='window_width', type=int, default=1000)
    parser.add_argument('--window-step', dest='window_step', type=int, default=900)
    parser.add_argument('--fragment-size', dest='fragment_size', type=int, default=1)
    parser.add_argument('--genome', dest='genome', type=str, default='hg19')
    args = parser.parse_args()
    args = parser.parse_args(["projects/chipmentation/results/periodicity-data",
                              "projects/chipmentation/results/periodicity",
                              "projects/chipmentation/results/plots/periodicity",
                              "data/human/chipmentation/mapped/merged/DNase_UWashington_K562_mergedReplicates.bam",
                              "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam",
                              "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_ChIP.bam",
                              "data/human/chipmentation/mapped/merged/H3K27me3_K562_10k500k_CM.bam",
                              "data/human/chipmentation/mapped/merged/H3K27me3_K562_500k_ChIP.bam",
                              "data/human/chipmentation/mapped/merged/PU1_K562_10mio_CM.bam",
                              "data/human/chipmentation/mapped/merged/CTCF_K562_10mio_CM.bam",
                              #  "/fhgfs/groups/lab_bock/arendeiro/projects/atac-seq/data/mapped/ASP14_50k_ATAC-seq_nan_nan_DOX_ATAC10-8_0_0.trimmed.bowtie2.shifted.dups.bam",
                              "/fhgfs/groups/lab_bock/arendeiro/projects/atac-seq/data/mapped/ASP14_50k_ATAC-seq_nan_nan_untreated_ATAC10-7_0_0.trimmed.bowtie2.shifted.dups.bam"
                              ]
                             )
    # TODO: mkdirs
    args.data_dir = os.path.abspath(args.data_dir)
    args.results_dir = os.path.abspath(args.results_dir)
    args.plots_dir = os.path.abspath(args.plots_dir)

    # Get regions of interest in the genome
    gapsRepeats = BedTool(os.path.join("/home", "arendeiro", "reference/Homo_sapiens/hg19_gapsRepeats.bed"))

    # Whole genome in 1kb-windows
    whole_genome = BedTool.window_maker(BedTool(), genome='hg19', w=args.window_width, s=args.window_width)
    whole_genome = whole_genome.intersect(b=gapsRepeats, v=True, wa=True)

    # DHS and non-DHS in 1kb-windows
    dhs = BedTool(os.path.join("/home", "arendeiro", "wgEncodeOpenChromDnaseK562Pk.narrowPeak"))
    dhs = dhs.intersect(b=gapsRepeats, v=True, wa=True)
    non_dhs = whole_genome.intersect(b=dhs, v=True)

    # Bed files
    H3K27me3 = BedTool(os.path.join("/home", "arendeiro", "wgEncodeSydhHistoneK562H3k27me3bUcdPk.narrowPeak"))
    H3K27me3 = H3K27me3.intersect(b=gapsRepeats, v=True, wa=True)
    H3K4me3 = BedTool(os.path.join("/home", "arendeiro", "wgEncodeSydhHistoneK562H3k4me3bUcdPk.narrowPeak"))
    H3K4me3 = H3K4me3.intersect(b=gapsRepeats, v=True, wa=True)

    # peaks not overlapping the other mark in 1kb-windows
    H3K27me3_only = H3K27me3.intersect(b=H3K4me3, v=True)
    H3K4me3_only = H3K4me3.intersect(b=H3K27me3, v=True)

    # peaks overlapping each other in 1kb-windows
    H3K27me3_H3K4me3 = H3K4me3.intersect(b=H3K27me3)

    # Make 1kb windows and convert to HTSeq
    whole_genome = bedToolsInterval2GenomicInterval(whole_genome, strand=False, name=False)
    non_dhs = bedToolsInterval2GenomicInterval(BedTool.window_maker(non_dhs, b=non_dhs, w=args.window_width, s=args.window_width), strand=False, name=False)
    dhs = bedToolsInterval2GenomicInterval(BedTool.window_maker(dhs, b=dhs, w=args.window_width, s=args.window_width), strand=False, name=False)
    H3K27me3_only = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K27me3_only, b=H3K27me3_only, w=args.window_width, s=args.window_width), strand=False, name=False)
    H3K4me3_only = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K4me3_only, b=H3K4me3_only, w=args.window_width, s=args.window_width), strand=False, name=False)
    H3K27me3_H3K4me3 = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K27me3_H3K4me3, b=H3K27me3_H3K4me3, w=args.window_width, s=args.window_width), strand=False, name=False)
    H3K27me3 = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K27me3, b=H3K27me3, w=args.window_width, s=args.window_width), strand=False, name=False)
    H3K4me3 = bedToolsInterval2GenomicInterval(BedTool.window_maker(H3K4me3, b=H3K4me3, w=args.window_width, s=args.window_width), strand=False, name=False)

    regions = {"whole_genome": whole_genome,
               "dhs": dhs,
               "non_dhs": non_dhs,
               "H3K27me3": H3K27me3,
               "H3K4me3": H3K4me3,
               "H3K27me3_only": H3K27me3_only,
               "H3K4me3_only": H3K4me3_only,
               "H3K27me3_H3K4me3": H3K27me3_H3K4me3
               }
    assert all([len(region) > 1 for region in regions.values()])
    pickle.dump(regions, open(os.path.join(args.data_dir, "genomic_regions.no_repeats.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    regions = pickle.load(open(os.path.join(args.data_dir, "genomic_regions.no_repeats.pickle"), "r"))
    # with repeats, load:
    # regions = pickle.load(open(os.path.join(args.data_dir, "genomic_regions.pickle"), "r"))

    # Initialize Slurm object
    slurm = DivideAndSlurm()
    tasks = dict()

    # Get sample names
    samples = {re.sub("\.bam", "", os.path.basename(sampleFile)): os.path.abspath(sampleFile) for sampleFile in args.bam_files}

    # Submit tasks for combinations of regions and bam files
    for regionName, region in regions.items():
        for sampleName, sampleFile in samples.items():
            exportName = os.path.join(args.data_dir, sampleName + "_" + regionName)
            if os.path.isfile(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle")):
                continue
            # Add new task
            task = CountDistances(region, 4, os.path.abspath(sampleFile), permute=False, queue="shortq", time="4:00:00", permissive=True, cpusPerTask=4)
            slurm.add_task(task)
            slurm.submit(task)  # Submit new task
            tasks[task] = (sampleName, regionName, False)  # Keep track
            print(tasks[task])
            if os.path.isfile(os.path.join(exportName + ".countsPermutedStranded-noRepeats-slurm.pickle")):
                continue
            # Add permuted
            task = CountDistances(region, 4, os.path.abspath(sampleFile), permute=True, queue="shortq", time="4:00:00", permissive=True, cpusPerTask=4)
            slurm.add_task(task)
            slurm.submit(task)  # Submit new task
            tasks[task] = (sampleName, regionName, True)  # Keep track
            print(tasks[task])

    stored = list()
    pickle.dump((slurm, tasks, stored), open("/home/arendeiro/slurm_20150210_1kb_norepeats.pickle", "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    slurm, tasks, stored = pickle.load(open("/home/arendeiro/slurm_20150210_1kb_norepeats.pickle", "r"))

    # Collect processed data
    for task, (sampleName, regionName, permuted) in tasks.items():          # loop through tasks, see if ready
        if task.is_ready():                                                 # if yes, collect output and save
            print(textwrap.dedent("""\
            Task {0} is now ready! {1}, {2}, {3}
            Time to completion was: {4} minutes.
            """.format(task, sampleName, regionName, permuted, int(time.time() - task.submission_time) / 60.)))
            exportName = os.path.join(args.data_dir, sampleName + "_" + regionName)
            dists = task.collect()
            if not permuted:
                pickle.dump(dists, open(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
            else:
                pickle.dump(dists, open(os.path.join(exportName + ".countsPermutedStranded-noRepeats-slurm.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
            stored.append(task)

    # For each signal extract most abundant periodic signal through FFT, IFFT
    for regionName, region in regions.items():
        for sampleName, sampleFile in samples.items():
            exportName = os.path.join(args.data_dir, sampleName + "_" + regionName)

            if not os.path.isfile(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle")):
                continue
            dists = pickle.load(open(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle"), "r"))

            distsPos = {dist: count for dist, count in dists.items() if dist >= 0}
            distsNeg = {abs(dist): count for dist, count in dists.items() if dist <= 0}

            if not os.path.isfile(os.path.join(exportName + ".countsPermutedStranded-noRepeats-slurm.pickle")):
                continue
            permutedDists = pickle.load(open(os.path.join(exportName + ".countsPermutedStranded-noRepeats-slurm.pickle"), "r"))

            permutedDistsPos = {dist: count for dist, count in permutedDists.items() if dist >= 0}
            permutedDistsNeg = {abs(dist): count for dist, count in permutedDists.items() if dist <= 0}

            # Extract most abundant periodic pattern from signal
            # for DNase, extract from a different window (70-150bp)
            patternPos = extractPattern(distsPos, range(60, 100), os.path.join(args.data_dir, exportName + "_posStrand"))
            patternNeg = extractPattern(distsNeg, range(60, 100), os.path.join(args.data_dir, exportName + "_negStrand"))
            pattern = extractPattern(Counter(distsPos) + Counter(distsNeg), range(60, 100), os.path.join(args.data_dir, exportName + "_bothStrands"))

            permutedPatternPos = extractPattern(permutedDistsPos, range(60, 100), os.path.join(args.data_dir, exportName + "_posStrand_permuted"))
            permutedPatternNeg = extractPattern(permutedDistsNeg, range(60, 100), os.path.join(args.data_dir, exportName + "_negStrand_permuted"))
            permutedPattern = extractPattern(Counter(permutedDistsPos) + Counter(permutedDistsNeg), range(60, 100), os.path.join(args.data_dir, exportName + "_bothStrands_permuted"))

    # Focus on H3K4me3 data and nucleosome positioning
    #
    #
    # calculate read coverage in H3K4me3 peaks
    sampleName = "H3K4me3_K562_500k_CM"
    sampleFile = samples[sampleName]
    regionName = "H3K4me3_only"

    peaks = BedTool("data/human/chipmentation/peaks/{sample}_peaks/{sample}_peaks.narrowPeak".format(sample=sampleName))
    peaks = bedToolsInterval2GenomicInterval(peaks)

    exportName = os.path.join(args.data_dir, sampleName + "_" + sampleName)

    task = Coverage(peaks, 4, os.path.abspath(sampleFile), fragment_size=args.fragment_size, strand_wise=True, queue="develop", cpusPerTask=2)
    slurm.add_task(task)
    slurm.submit(task)
    coverage = task.collect()

    pickle.dump(coverage, open(os.path.join(args.data_dir, exportName + ".peakCoverageStranded.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    coverage = pickle.load(open(os.path.join(args.data_dir, exportName + ".peakCoverageStranded.pickle"), "r"))

    # coverage = {name : reads for name, reads in coverage.items() if len(reads) >= 1000} # Filter out peaks smaller than 1kb

    # Correlate coverage and signal pattern
    # get pattern
    dists = pickle.load(open(os.path.join(args.data_dir, sampleName + "_" + regionName + ".countsStranded-slurm.pickle"), "r"))
    distsPos = {dist: count for dist, count in dists.items() if dist >= 0}
    distsNeg = {abs(dist): count for dist, count in dists.items() if dist <= 0}
    pattern = extractPattern(Counter(distsPos) + Counter(distsNeg), range(60, 100),
                             os.path.join(args.data_dir, sampleName + "_" + regionName + "_bothStrands")
                             )

    # positive strand
    coveragePos = {peak: reads[0] for peak, reads in coverage.items()}
    taskPos = Correlation(coveragePos, 2, pattern, step=1, queue="develop", cpusPerTask=8)
    slurm.add_task(taskPos)
    slurm.submit(taskPos)
    correlationsPos = taskPos.collect()

    # negative strand
    coverageNeg = {peak: reads[1] for peak, reads in coverage.items()}
    taskNeg = Correlation(coverageNeg, 2, pattern, step=1, queue="develop", cpusPerTask=8)
    slurm.add_task(taskNeg)
    slurm.submit(taskNeg)
    correlationsNeg = taskNeg.collect()

    pickle.dump((correlationsPos, correlationsNeg), open(os.path.join(args.data_dir, exportName + ".peakCorrelationStranded.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    correlationsPos, correlationsNeg = pickle.load(open(os.path.join(args.data_dir, exportName + ".peakCorrelationStranded.pickle"), "r"))

    # Binarize data for HMM
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
    binary = {peak: getConsensus(pos[peak], neg[peak]) for peak in pos.keys()}

    # get zeroes for rest of the genome

    # HMM
    model = WinterHMM()
    # Train
    model.train(binary.values()[:100])  # subset data for training
    # see new probabilities
    [(s.name, s.distribution) for s in model.model.states]
    model.model.dense_transition_matrix()
    # save model with trained parameters
    pickle.dump(model, open(os.path.join(args.results_dir, "hmm_trained_with_" + sampleName + ".pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    model = pickle.load(open(os.path.join(args.results_dir, "hmm_trained_with_" + sampleName + ".pickle"), "r"))
    # Predict
    hmmOutput = {peak: model.predict(sequence) for peak, sequence in binary.items()}

    # Get DARNS from hmmOutput
    DARNS = getDARNS(peaks, hmmOutput)

    # Plot attributes
    widths, distances, midDistances = measureDARNS(DARNS)

    plt.plot(widths.keys(), widths.values(), 'o')
    plt.plot(distances.keys(), distances.values(), 'o')
    plt.plot(midDistances.keys(), midDistances.values(), 'o')
    sns.violinplot(widths)

    # Get scores from DARNS
    probs = {peak: model.retrieveProbabilities(sequence) for peak, sequence in binary.items()}

    # plot predicted darns and post probs
    i = 3
    sequence = hmmOutput.values()[i]

    limits = getDARNSlimits(sequence)
    probs = model.retrieveProbabilities(binary.items()[i][1])

    from matplotlib.patches import Rectangle
    colors = sns.color_palette('deep', n_colors=6, desat=0.5)
    sns.set_context(rc={"figure.figsize": (14, 6)})
    sns.plt.axhline(y=1.1, c=colors[0], alpha=0.7)
    sns.plt.xlim([1, len(sequence) + 1])
    sns.plt.ylim([0, 1.2])
    sns.plt.ylabel(r'posterior probs, $\gamma_k$')
    sns.plt.xlabel(r'$k$')
    axis = sns.plt.gca()

    # Plot viterbi predicted TMDs
    for start, end in limits:
        axis.add_patch(Rectangle((start + 1, 1.075), end - start + 1, 0.05,
                                 facecolor=colors[0], alpha=0.7))

    # Get indices of states
    indices = {state.name: i for i, state in enumerate(model.model.states)}

    sns.plt.plot(range(1, len(sequence) + 1), probs,
                 c=colors[2], alpha=0.7)
    plt.show()

    ### Output bed/wig




    # Export wig files with raw correlations
    exportWigFile(
        [peaks[i] for i in correlationsPos.keys()],
        correlationsPos.values(),
        len(pattern) / 2,
        os.path.join(args.results_dir, sampleName + ".peakCorrelationPos.wig"),
        sampleName + " raw absolute correlation - positive strand"
    )
    exportWigFile(
        [peaks[i] for i in correlationsNeg.keys()],
        correlationsNeg.values(),
        len(pattern) / 2,
        os.path.join(args.results_dir, sampleName + ".peakCorrelationNeg.wig"),
        sampleName + " raw absolute correlation - negative strand"
    )
    # Export wig files with correlations peaks

    # Export wig files with DARNS

    # Export bed files with DARNS


    # Genome arithmetics with DARNS:
    # Plot :
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
    #     os.path.join(args.results_dir, sampleName + ".correlationPeaksPos.bed"),
    #     sampleName + " correlation peaks - positive strand"
    # )

    # corPeaksNeg = dict()
    # for name, cor in correlationsNeg.items():
    #     corPeak = findPeaks(abs(cor), 50)
    #     corPeaksNeg[name] = corPeak

    # exportBedFile(
    #     {name : peaks[name] for name, peak in corPeaksNeg.items()},
    #     corPeaksNeg,
    #     len(pattern)/2,
    #     os.path.join(args.results_dir, sampleName + ".correlationPeaksNeg.bed"),
    #     sampleName + " correlation peaks - negative strand"
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


    # find peaks
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
    #     os.path.join(args.results_dir, sampleName + ".peakCorrelation.Winter.wig"),
    #     sampleName
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


class CountDistances(Task):
    """
    Task to perform counting of distances between reads under regions.
    """
    def __init__(self, data, fractions, *args, **kwargs):
        super(CountDistances, self).__init__(data, fractions, *args, **kwargs)
        # Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1, 1000))], sep="_")
        self.name = "count_distances_{0}".format(now)

        # Parse
        # required argument
        if len(self.args) != 1:
            raise TypeError("Bam file argument is missing")
        self.bam_file = self.args[0]
        # additional arguments
        if "strand_wise" in kwargs.keys():
            self.strand_wise = kwargs["strand_wise"]
        else:
            self.strand_wise = True
        if "duplicates" in kwargs.keys():
            self.duplicates = kwargs["duplicates"]
        else:
            self.duplicates = True
        if "permute" in kwargs.keys():
            self.permute = kwargs["permute"]
        else:
            self.permute = False
        if "fragment_size" in kwargs.keys():
            self.fragment_size = kwargs["fragment_size"]
        else:
            self.fragment_size = 1

    def _prepare(self):
        """
        Add task to be performed with data. Is called when task is added to DivideAndSlurm object.
        """
        self.log = os.path.join(self.slurm.logDir, string.join([self.name, "log"], sep="."))  # add abspath

        # Split data in fractions
        ids, groups, files = self._split_data()

        ### Make jobs with groups of data
        self.jobs = list(); self.jobFiles = list(); self.inputPickles = list(); self.outputPickles = list()

        # for each group of data
        for i in xrange(len(ids)):
            jobFile = files[i] + "_count_distances.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            # assemble job file
            # header
            job = self._slurmHeader(ids[i])

            # command - add abspath!
            task = """\
                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python count_distances_parallel.py {0} {1} {2} """.format(inputPickle, outputPickle, self.bam_file)

            if self.strand_wise:
                task += "--strand-wise "
            if self.duplicates:
                task += "--duplicates "
            if self.permute:
                task += "--permute "
            task += "--fragment-size {0}".format(self.fragment_size)

            task += """

                # Deactivate virtual environment
                deactivate
                    """

            job += textwrap.dedent(task)

            # footer
            job += self._slurmFooter()

            # add to save attributes
            self.jobs.append(job)
            self.jobFiles.append(jobFile)
            self.inputPickles.append(inputPickle)
            self.outputPickles.append(outputPickle)

            # write job file to disk
            with open(jobFile, 'w') as handle:
                handle.write(textwrap.dedent(job))

        # Delete data if jobs are ready to submit and data is serialized
        if hasattr(self, "jobs") and hasattr(self, "jobFiles"):
            del self.data

    def collect(self):
        """
        If self.is_ready(), return joined reduced data.
        """
        if hasattr(self, "output"):  # if output is already stored, just return it
            return self.output

        if self.is_ready():
            # load all pickles into list
            if self.permissive:
                outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
            else:
                outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
            # if all are counters, and their elements are counters, sum them
            if all([type(outputs[i]) == Counter for i in range(len(outputs))]):
                output = reduce(lambda x, y: x + y, outputs)  # reduce
                if type(output) == Counter:
                    self.output = output    # store output in object
                    self._rm_temps()  # delete tmp files
                    return self.output
        else:
            raise TypeError("Task is not ready yet.")


class Coverage(Task):
    """
    Task to get read coverage under regions.
    """
    def __init__(self, data, fractions, *args, **kwargs):
        super(Coverage, self).__init__(data, fractions, *args, **kwargs)
        # Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1, 1000))], sep="_")
        self.name = "coverage_{0}".format(now)

        # Parse
        # required argument
        if len(self.args) != 1:
            raise TypeError("Bam file argument is missing")
        self.bam_file = self.args[0]
        # additional arguments
        if "strand_wise" in kwargs.keys():
            self.strand_wise = kwargs["strand_wise"]
        else:
            self.strand_wise = True
        if "duplicates" in kwargs.keys():
            self.duplicates = kwargs["duplicates"]
        else:
            self.duplicates = True
        if "orientation" in kwargs.keys():
            self.orientation = kwargs["orientation"]
        else:
            self.orientation = False
        if "permute" in kwargs.keys():
            self.permute = kwargs["permute"]
        else:
            self.permute = False
        if "fragment_size" in kwargs.keys():
            self.fragment_size = kwargs["fragment_size"]
        else:
            self.fragment_size = 1

    def _prepare(self):
        """
        Add task to be performed with data. Is called when task is added to DivideAndSlurm object.
        """
        self.log = os.path.join(self.slurm.logDir, string.join([self.name, "log"], sep="."))  # add abspath

        # Split data in fractions
        ids, groups, files = self._split_data()

        # Make jobs with groups of data
        self.jobs = list(); self.jobFiles = list(); self.inputPickles = list(); self.outputPickles = list()

        # for each group of data
        for i in xrange(len(ids)):
            jobFile = files[i] + "_coverage.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            # assemble job file
            # header
            job = self._slurmHeader(ids[i])

            # command - add abspath!
            task = """\
                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python coverage_parallel.py {0} {1} {2} """.format(inputPickle, outputPickle, self.bam_file)

            if self.strand_wise:
                task += "--strand-wise "
            if self.duplicates:
                task += "--duplicates "
            if self.orientation:
                task += "--orientation "
            if self.permute:
                task += "--permute "
            task += "--fragment-size {0}".format(self.fragment_size)

            task += """

                # Deactivate virtual environment
                deactivate
                    """

            job += textwrap.dedent(task)

            # footer
            job += self._slurmFooter()

            # add to save attributes
            self.jobs.append(job)
            self.jobFiles.append(jobFile)
            self.inputPickles.append(inputPickle)
            self.outputPickles.append(outputPickle)

            # write job file to disk
            with open(jobFile, 'w') as handle:
                handle.write(textwrap.dedent(job))

        # Delete data if jobs are ready to submit and data is serialized
        if hasattr(self, "jobs") and hasattr(self, "jobFiles"):
            del self.data

    def collect(self):
        """
        If self.is_ready(), return joined reduced data.
        """
        if hasattr(self, "output"):  # if output is already stored, just return it
            return self.output

        if self.is_ready():
            # load all pickles into list
            if self.permissive:
                outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
            else:
                outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
            # if all are counters, and their elements are counters, sum them
            if all([type(outputs[i]) == dict for i in range(len(outputs))]):
                output = reduce(lambda x, y: dict(x, **y), outputs)
                if type(output) == dict:
                    self.output = output  # store output in object
                    self._rm_temps()  # delete tmp files
                    return self.output
        else:
            raise TypeError("Task is not ready yet.")


class Correlation(Task):
    """
    Task to get read coverage under regions.
    """
    def __init__(self, data, fractions, *args, **kwargs):
        super(Correlation, self).__init__(data, fractions, *args, **kwargs)
        # Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1, 1000))], sep="_")
        self.name = "correlation_{0}".format(now)

        # Parse
        # required argument
        if len(self.args) != 1:
            raise TypeError("Pattern argument is missing")
        self.pattern = self.args[0]
        # additional arguments
        if "step" in kwargs.keys():
            self.step = kwargs["step"]
        else:
            self.step = 1

    def _prepare(self):
        """
        Add task to be performed with data. Is called when task is added to DivideAndSlurm object.
        """
        self.log = os.path.join(self.slurm.logDir, string.join([self.name, "log"], sep="."))  # add abspath

        # Pickle pattern
        self.patternPickle = os.path.join(self.slurm.tmpDir, self.name + "_pattern.pickle")
        pickle.dump(self.pattern, open(self.patternPickle, "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        # Split data in fractions
        ids, groups, files = self._split_data()

        # Make jobs with groups of data
        self.jobs = list(); self.jobFiles = list(); self.inputPickles = list(); self.outputPickles = list()

        # for each group of data
        for i in xrange(len(ids)):
            jobFile = files[i] + "_correlation.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            # assemble job file
            # header
            job = self._slurmHeader(ids[i])

            # command - add abspath to python script!
            task = """\
                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python correlation_parallel.py {0} {1} {2} --step {3}

                # Deactivate virtual environment
                deactivate
                """.format(inputPickle, outputPickle, self.patternPickle, self.step)

            job += textwrap.dedent(task)

            # footer
            job += self._slurmFooter()

            # add to save attributes
            self.jobs.append(job)
            self.jobFiles.append(jobFile)
            self.inputPickles.append(inputPickle)
            self.outputPickles.append(outputPickle)

            # write job file to disk
            with open(jobFile, 'w') as handle:
                handle.write(textwrap.dedent(job))

        # Delete data if jobs are ready to submit and data is serialized
        if hasattr(self, "jobs") and hasattr(self, "jobFiles"):
            del self.data

    def collect(self):
        """
        If self.is_ready(), return joined reduced data.
        """
        if hasattr(self, "output"):  # if output is already stored, just return it
            return self.output

        if self.is_ready():
            # load all pickles into list
            if self.permissive:
                outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
            else:
                outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
            # if all are counters, and their elements are counters, sum them
            if all([type(outputs[i]) == dict for i in range(len(outputs))]):
                output = reduce(lambda x, y: dict(x, **y), outputs)
                if type(output) == dict:
                    self.output = output  # store output in object
                    self._rm_temps()  # delete tmp files
                    return self.output
        else:
            raise TypeError("Task is not ready yet.")


class WinterHMM(object):
    """
    Class to model the Hidden Markov model described in the Winter 2013 paper.
    """
    def __init__(self):
        super(WinterHMM, self).__init__()
        background = yahmm.DiscreteDistribution({
            1: 0.33,
            0: 0.33,
            -1: 0.33
        })
        pos = yahmm.DiscreteDistribution({
            1: 0.7,
            0: 0.2,
            -1: 0.1
        })
        neg = yahmm.DiscreteDistribution({
            1: 0.1,
            0: 0.2,
            -1: 0.7
        })
        zero = yahmm.DiscreteDistribution({
            1: 0,
            0: 1,
            -1: 0
        })

        # Create states
        b = yahmm.State(background, name="B")      
        minus = yahmm.State(neg, name="-")
        in1 = yahmm.State(zero, name="in1")
        in2 = yahmm.State(zero, name="in2")
        in3 = yahmm.State(zero, name="in3")
        plus = yahmm.State(pos, name="+")
        in4 = yahmm.State(zero, name="in4")
        in5 = yahmm.State(zero, name="in5")
        in6 = yahmm.State(zero, name="in6")
        in7 = yahmm.State(zero, name="in7")
        in8 = yahmm.State(zero, name="in8")
        in9 = yahmm.State(zero, name="in9")
        in10 = yahmm.State(zero, name="in10")
        in11 = yahmm.State(zero, name="in11")

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

        self.model.add_transition(self.model.start, b, 0.5)       # start in background
        self.model.add_transition(self.model.start, minus, 0.25)  # start in minus
        self.model.add_transition(self.model.start, plus, 0.25)   # start in plus

        self.model.add_transition(b, b, 0.8)                      # stay in background
        self.model.add_transition(b, minus, 0.1)                  # enter in minus
        self.model.add_transition(b, plus, 0.1)                   # enter in plus

        self.model.add_transition(minus, b, 0.25)                 # can exit from plus
        self.model.add_transition(minus, in1, 0.25)
        self.model.add_transition(minus, plus, 0.25)

        self.model.add_transition(in1, in2, 0.5)                  # states 1-3 can go directly to plus or follow a string 1 to 4
        self.model.add_transition(in1, plus, 0.5)
        self.model.add_transition(in2, in3, 0.5)
        self.model.add_transition(in2, plus, 0.5)
        self.model.add_transition(in3, in4, 0.5)
        self.model.add_transition(in3, plus, 0.5)

        self.model.add_transition(plus, b, 0.50)                  # can exit from plus
        self.model.add_transition(plus, in4, 0.50)
        self.model.add_transition(in4, in5, 1)                    # string of 4 to 8 at least
        self.model.add_transition(in5, in6, 1)
        self.model.add_transition(in6, in7, 1)
        self.model.add_transition(in7, in8, 1)
        self.model.add_transition(in8, in9, 0.5)                  # states 8-11 can go directly to minus
        self.model.add_transition(in8, minus, 0.5)
        self.model.add_transition(in9, in10, 0.5)
        self.model.add_transition(in9, minus, 0.5)
        self.model.add_transition(in10, in11, 0.5)
        self.model.add_transition(in10, minus, 0.5)
        self.model.add_transition(in11, minus, 1)

        self.model.bake()

    def draw(self):
        """
        Draws the Markov chain of the model.
        """
        return self.model.draw(node_size=400, labels={state.name: str(state.name) for state in self.model.states}, font_size=20)

    def train(self, sequences):
        """
        Train model with given data.
        """
        self.model.train(sequences)

    def predict(self, observations):
        """
        Predict hidden states from observations.
        """
        logp, chain = self.model.maximum_a_posteriori(observations)
        return [state.name for num, state in chain][1:-1]  # strip begin and end states

    def retrieveProbabilities(self, observations):
        """
        Retrieve the posterior probabilities of being in a "nucleossome" state for a sequence of observations.
        """
        trans, ems = self.model.forward_backward(observations)
        ems = np.exp(ems)
        self.probs = ems / np.sum(ems, axis=1)[:, np.newaxis]  # probs of all states
        background_state = 0
        prob_nucleossome = 1 - self.probs[:, background_state]  # prob of all states but background
        return prob_nucleossome


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


def bedToolsInterval2GenomicInterval(bedtool, strand=True, name=True):
    """
    Given a pybedtools.BedTool object, return dictionary of HTSeq.GenomicInterval objects.

    bedtool=pybedtools.BedTool object.
    """
    intervals = OrderedDict()
    if strand:
        for iv in bedtool:
            if name:
                intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, iv.strand)
            else:
                intervals[string.join(iv.fields[:3], sep="_")] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, iv.strand)
    else:
        for iv in bedtool:
            if name:
                intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end)
            else:
                intervals[string.join(iv.fields[:3], sep="_")] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end)

    return intervals


def distances(bam, intervals, fragmentsize, duplicates=True, strand_wise=True, permutate=False):
    """
    Gets pairwise distance between reads in intervals. Returns dict with distance:count.
    If permutate=True, it will randomize the reads in each interval along it.

    bam=HTSeq.BAM_Reader object.
    intervals=dict - HTSeq.GenomicInterval objects as values.
    fragmentsize=int.
    duplicates=bool.
    strand_wise=bool.
    permutate=bool.
    """
    dists = Counter()
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

    # for name, feature in intervals.iteritems():
    for name, feature in intervals.items():
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
            if not strand_wise and aln1.iv.strand != aln2.iv.strand:
                continue
            # adjust fragment to size
            aln1.iv.length = fragmentsize
            aln2.iv.length = fragmentsize

            # get position relative
            dist = abs(aln1.iv.start_d - aln2.iv.start_d)
            # add +1 to dict
            if strand_wise:
                if aln1.iv.strand == "+":
                    dists[dist] += 1
                if aln1.iv.strand == "-":
                    dists[-dist] += 1
            else:
                dists[dist] += 1
    return dists


def coverageInWindows(bam, intervals, fragmentsize, orientation=False, duplicates=True, strand_wise=False, permutate=False):
    """
    Gets read coverage in intervals. Returns dict of regionName:numpy.array if strand_wise=False,
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
        if not strand_wise:
            profile = np.zeros(feature.length, dtype=np.float64)
        else:
            profile = np.zeros((2, feature.length), dtype=np.float64)

        # Check if feature is in bam index
        if feature.chrom not in chroms or feature.chrom == "chrM":
            i += 1
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
            aln.iv.length = fragmentsize  # adjust to size

            # get position in relative to window
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

        # append feature profile to dict
        coverage[name] = profile
        i += 1
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
    p1 = np.poly1d(np.polyfit(x, y, 1))  # fit linear regression
    m, b = p1.coeffs
    plt.plot(
        x, y, 'o',
        x, p1(x), "-"
    )

    # restrict to signal
    x = distRange
    y = [dists[i] for i in x]
    p1 = np.poly1d(np.polyfit(x, y, 1))  # fit linear regression
    m, b = p1.coeffs
    plt.subplot(212)
    plt.plot(
        x, y, 'o',
        x, p1(x), "-"
    )
    plt.savefig(filePrefix + ".read_distances.pdf", bbox_inches='tight')
    plt.close()

    # measure distance to regression in that point
    distsReg = [y[i] - (m * i + b) for i in range(len(x))]

    # subtract value of minimum to all points
    distsReg -= min(distsReg)

    # fourier transform data
    time = x
    signal = distsReg

    # get frequencies from decomposed fft
    W = np.fft.fftfreq(signal.size, d=time[1] - time[0])
    f_signal = np.fft.fft(signal)

    # plot all frequencies, ask what is the amplitude of the signals
    freqs = dict()
    plt.figure(1)
    for i in np.arange(0, len(x) / 10., 0.01):
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
    plt.savefig(filePrefix + ".fft_filter-{0}bp_ifft.pdf".format(int(1 / top)), bbox_inches='tight')
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
            R.append(pearsonr(pattern, profile[i:i + len(pattern)])[0])
            i += 1
        return np.array(R)


def fitKDE(X, bandwidth):
    """
    Fit gaussian kernel density estimator to array, return

    X=numerical iterable.
    bandwidth=kernel bandwidth - float.
    """
    return kde.kdensity(range(len(X)), bw=bandwidth, weights=X)[0]


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
    return {peak: (peak - peakWidth / 2, peak + peakWidth / 2) for peak in peaks}


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
        if maximas[i] == True and X[i] > 0.05 and any(n < 0 for n in X[prev_max: i]):  # watch out for namespace pollution with np.any
            # if from i to next_max there is no neg, then is max
            for j in xrange(i, len(X)):
                if maximas[j] == True and X[j] > 0.05:
                    next_max = j
                    break

            if not any(n > 0 for n in X[i: next_max]):
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


def getDARNS(intervals, hmmOutput):
    """
    Returns list of HTSeq.GenomicInterval objects for all DARNS in hmmOutput (a dict of name:sequence).
    DARNS are stretches of sequences which belong to the "nucleosome" states (see WinterHMM object definition).

    DARNS cannot:
        - start at the last element of the sequence.
        - be of length <= 1

    intervals=dict - name:HTSeq.GenomicInterval objects.
    hmmOutput=dict - name:iterable with sequence from HMM.
    """
    outside = ["B"]
    DARNS = list()
    for name, sequence in hmmOutput.items():
        DARN = False
        for i in xrange(len(sequence)):
            if not DARN and sequence[i] not in outside and i != len(sequence):
                DARN = True
                DARN_start = i
            if DARN and sequence[i] in outside and i - DARN_start > 1:
                DARN = False
                DARN_end = i
                DARNS.append([
                    intervals[name].chrom,
                    intervals[name].start + DARN_start,
                    intervals[name].start + DARN_end,
                    name + "_DARN{0}".format(i)
                ])
    return DARNS


def getDARNSlimits(sequence):
    """
    sequence=iterable with sequence from HMM.
    """
    outside = ["B"]
    DARN = False
    starts = list()
    ends = list()
    for i in xrange(len(sequence)):
        if not DARN and sequence[i] not in outside:
            DARN = True
            starts.append(i)
        if DARN and sequence[i] in outside:
            DARN = False
            ends.append(i)
    return zip(starts, ends)


def measureDARNS(DARNS):
    """
    Computes distributions of DARNS attributes: width, interdistance, distance between midpoints.

    DARNS=iterable of list with [chrom, start, end, name] of DARN.
    """
    chrom, start, end, name = range(4)
    width = lambda x: x[end] - x[start]
    widths = Counter()
    distances = Counter()
    midDistances = Counter()

    for d1, d2 in itertools.permutations(DARNS, 2):
        if d1[chrom] == d2[chrom]:
            # widths
            widths[width(d1)] += 1
            widths[width(d2)] += 1

            # distance end-to-end
            if d1[end] <= d2[start]:
                distances[abs(d2[start] - d1[end])] += 1
            else:
                distances[d1[start] - d2[end]] += 1

            # distance midpoint-midpoint
            midDistances[abs((d1[end] - d1[start]) - (d2[end] - d2[start]))] += 1
    return (widths, distances, midDistances)


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
                # TODO: check peak boundaries
                entry = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    intervals[name].chrom,
                    intervals[name].start + offset + tup[0],
                    intervals[name].start + offset + tup[1],
                    name,
                    1,      # score
                    strand  # strand
                )
                handle.write(entry)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(0)
