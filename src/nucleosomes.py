#!/usr/env python

"""
Predicting nucleosome associated-regions from ChIPmentation data.

In brief, this does:
    * Pairwise distance count between reads from several samples (from several
    techniques - ChIP, ChIPmentation, DNase-seq, ATAC-seq) in several genomic
    locations (H3K4me3/H3K27me3 peaks, DHS regions, whole genome, etc...) in strand-
    -specific way and also for permuted reads.
    * Get pattern from the read distances distribution (correlogram).
    * Decompose pattern into signals with various frequencies (using FFT) and retrieve
    most abundant (10bp).
    * Calculate correlations between pattern and ChIPmentation read coverage/permuted reads along the genome.
    * Extract local maxima from correlation into a binary signal.
    * Feed to a HMM modeling a nucleosome, output predicted nucleosome-associated regions.
This will do:
    * Measure several features of predicted regions.
    * Plot and calculate p-values for each feature between real data and permuted.
    * Train multivariate linear regression (logistic) classifier with features from real and permuted data.
    * Classify the rest of data 10 times independently and calculate FDR for each nucleosome-associated region.
    * Export only nucleosome-associated regions with FDR<0.5.

    * ... plus several plots along the way.

TODO:
    * Plot CM, DNase, MNase in DARNS:
        * From the middle of the DARN (mid-peak)
        * From the 5' and 3' end of nucleosome Dyads (get external data)
        * DARNS frequency around TSSs (models and CAGE) and TTSs
        * DARNS frequency around CpGs islands

    * Try:
        * Do the same on concatenated histone ChIPmentation data.
        * Do the same using IGG data as background.
"""


from argparse import ArgumentParser
from collections import Counter, OrderedDict
from divideAndSlurm import DivideAndSlurm, Task
from matplotlib import pyplot as plt
from pybedtools import BedTool
from scipy import signal
import cPickle as pickle
import HTSeq
import itertools
import numpy as np
import pandas as pd
import os
import sys
import random
import re
import seaborn as sns
import string
import textwrap
import time
import yahmm
from sklearn.linear_model import LinearRegression

np.set_printoptions(linewidth=200)


def main(args):
    # Parse command-line arguments
    parser = ArgumentParser()

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

    # Skip creating regions if exist (dedicated to Nathan)
    regionsPickle = os.path.join(args.data_dir, "genomic_regions.no_repeats.pickle")
    if os.path.isfile(regionsPickle):
        regions = pickle.load(open(regionsPickle), "r")
    else:
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
        pickle.dump(regions, open(os.path.join(args.data_dir, "genomic_regions.no_repeats.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    assert all([len(region) > 1 for region in regions.values()])

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
            """.format(task, sampleName, regionName, permuted, int(time.time() - task.submissiontime) / 60.)))
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

            try:
                dists = pickle.load(open(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle"), "r"))
            except IOError("Can't open file."):
                continue

            distsPos = {dist: count for dist, count in dists.items() if dist >= 0}
            distsNeg = {abs(dist): count for dist, count in dists.items() if dist <= 0}

            # Extract most abundant periodic pattern from signal
            # for DNase, extract from a different window (70-150bp)
            # patternPos = extractPattern(distsPos, range(60, 100), os.path.join(args.data_dir, exportName + "_posStrand-noRepeats"))
            # patternNeg = extractPattern(distsNeg, range(60, 100), os.path.join(args.data_dir, exportName + "_negStrand-noRepeats"))
            pattern = extractPattern(Counter(distsPos) + Counter(distsNeg), range(60, 100), os.path.join(args.data_dir, exportName + "_bothStrands-noRepeats"))

            try:
                permutedDists = pickle.load(open(os.path.join(exportName + ".countsPermutedStranded-noRepeats-slurm.pickle"), "r"))
            except IOError("Can't open file."):
                continue

            permutedDistsPos = {dist: count for dist, count in permutedDists.items() if dist >= 0}
            permutedDistsNeg = {abs(dist): count for dist, count in permutedDists.items() if dist <= 0}

            # Extract most abundant periodic pattern from signal
            # for DNase, extract from a different window (70-150bp)
            permutedPatternPos = extractPattern(permutedDistsPos, range(60, 100), os.path.join(args.data_dir, exportName + "_posStrand_permuted-noRepeats"))
            permutedPatternNeg = extractPattern(permutedDistsNeg, range(60, 100), os.path.join(args.data_dir, exportName + "_negStrand_permuted-noRepeats"))
            permutedPattern = extractPattern(Counter(permutedDistsPos) + Counter(permutedDistsNeg), range(60, 100), os.path.join(args.data_dir, exportName + "_bothStrands_permuted-noRepeats"))

    # Focus on H3K4me3 data and nucleosome positioning
    #
    #
    # calculate read coverage in H3K4me3 peaks
    samples = {re.sub("\.bam", "", os.path.basename(sampleFile)): os.path.abspath(sampleFile) for sampleFile in args.bam_files}
    sampleName = "H3K4me3_K562_500k_CM"
    sampleFile = samples[sampleName]
    regionName = "H3K4me3_only"
    exportName = os.path.join(args.data_dir, sampleName + "_" + regionName)

    # Correlate coverage and signal pattern
    # get pattern
    dists = pickle.load(open(os.path.join(exportName + ".countsStranded-noRepeats-slurm.pickle"), "r"))
    distsPos = {dist: count for dist, count in dists.items() if dist >= 0}
    distsNeg = {abs(dist): count for dist, count in dists.items() if dist <= 0}
    pattern = extractPattern(Counter(distsPos) + Counter(distsNeg), range(60, 100),
                             os.path.join(args.data_dir, exportName + "_bothStrands-noRepeats"))

    # Task get coverage and correlate with pattern, separate strands
    slurm = DivideAndSlurm()
    # make windows genome-wide with overlapping size of pattern
    width = 1000
    step = width - len(pattern)
    genome_windows = BedTool.window_maker(BedTool(), g={"chr1": (1, 249250621)}, w=width, s=step)  # genome_windows = BedTool.window_maker(BedTool(), genome='hg19', w=width, s=step)
    genome_windows = bedToolsInterval2GenomicInterval(genome_windows, strand=False, name=False)

    task = CorrelatePatternBam(genome_windows.items(), 40, pattern, os.path.abspath(sampleFile), cpusPerTask=8)
    slurm.add_task(task)
    slurm.submit(task)
    taskP = CorrelatePatternBam(genome_windows.items(), 40, pattern, os.path.abspath(sampleFile), permute=True, cpusPerTask=8)
    slurm.add_task(taskP)
    slurm.submit(taskP)

    # collect and serialize
    binary = task.collect()
    assert len(genome_windows) == len(binary)  # compare sizes
    assert all([len(window) == width - len(pattern) + 1 for window in binary.values()])  # likely to fail - patched downstream
    pickle.dump(binary, open(os.path.join(args.data_dir, exportName + ".peakCorrelationBinary.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    binaryP = taskP.collect()
    assert len(genome_windows) == len(binaryP)  # compare sizes
    assert all([len(window) == width - len(pattern) + 1 for window in binaryP.values()])  # likely to fail - patched downstream
    pickle.dump(binaryP, open(os.path.join(args.data_dir, exportName + ".peakCorrelationBinaryPermuted.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    binary = pickle.load(open(os.path.join(args.data_dir, exportName + ".peakCorrelationBinary.pickle"), "r"))
    binaryP = pickle.load(open(os.path.join(args.data_dir, exportName + ".peakCorrelationBinaryPermuted.pickle"), "r"))

    # to subset, pass: OrderedDict(sorted((binary.items()[9000:10000]))) <- 1Mb of sequence
    genome_binary = concatenateBinary(binary, len(pattern))  # 40
    pickle.dump(genome_binary, open(os.path.join(args.data_dir, exportName + ".peakCorrelationBinaryConcatenated.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    genome_binaryP = concatenateBinary(binaryP, len(pattern))
    pickle.dump(genome_binaryP, open(os.path.join(args.data_dir, exportName + ".peakCorrelationBinaryConcatenatedPermuted.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    genome_binary = pickle.load(open(os.path.join(args.data_dir, exportName + ".peakCorrelationBinaryConcatenated.pickle"), "r"))
    genome_binaryP = pickle.load(open(os.path.join(args.data_dir, exportName + ".peakCorrelationBinaryConcatenatedPermuted.pickle"), "r"))

    # HMM
    model = WinterHMM()

    # Model training
    #
    # Train on subset of data, see probabilities
    # i = len(genome_binary)
    # model.train([genome_binary.values()[0][:10000]])  # subset data for training
    # # see new probabilities
    # [(s.name, s.distribution) for s in model.model.states]  # emission
    # print(model.model.dense_transition_matrix())  # transition

    # # save model with trained parameters
    # pickle.dump(model, open(os.path.join(args.results_dir, sampleName + "_hmModel_trained_%i.pickle" % i), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    # model = pickle.load(open(os.path.join(args.results_dir, sampleName + "_hmModel_trained_%i.pickle" % i), "r"))

    # Predict and get DARNS
    DARNS = dict()
    ite = range(0, len(genome_binary.values()[0]), 5000000)
    prev = 0
    last = -1
    for cur in ite[1:]:
        print(cur)
        if not cur == last:
            hmmOutput = {chrom: model.predict(sequence[prev:cur]) for chrom, sequence in genome_binary.items()}
        else:
            hmmOutput = {chrom: model.predict(sequence[prev:last]) for chrom, sequence in genome_binary.items()}
        # add darns to dict
        for chrom, sequence in hmmOutput.items():
            DARNS[chrom] = getDARNS(sequence)
        prev = cur
    pickle.dump(DARNS, open(os.path.join(args.results_dir, sampleName + "_DARNS.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # predict from permuted data
    DARNSP = dict()
    ite = range(0, len(genome_binaryP.values()[0]), 5000000)
    prev = 0
    last = -1
    for cur in ite[1:]:
        print(cur)
        if not cur == last:
            hmmOutputP = {chrom: model.predict(sequence[prev:cur]) for chrom, sequence in genome_binaryP.items()}
        else:
            hmmOutputP = {chrom: model.predict(sequence[prev:last]) for chrom, sequence in genome_binaryP.items()}
        # add darns to dict
        for chrom, sequence in hmmOutputP.items():
            DARNSP[chrom] = getDARNS(sequence)
        prev = cur
    pickle.dump(DARNSP, open(os.path.join(args.results_dir, sampleName + "_DARNSPermuted.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # Measure attributes:
    # value of maximum positive strand correlation peak
    # value of maximum negative strand correlation peak
    # mean of positive strand correlation peak values
    # mean of negative strand correlation peak values
    # sum of positive strand correlation peak values
    # sum of negative strand correlation peak values

    # number of cycles through HMM “nucleosome” states (0.5 for each peak state passed)
    # the HMM posterior probability of the path through the background state
    # and the HMM posterior probability of the path through the nucleosome set of states

    # length (bp)
    # sum of reads within the DARNS
    # read density
    # spacing to nearest upstream DARNS (bp)
    # spacing to nearest downstream DARNS (bp)

    # Measure read count and density in DARNS
    slurm = DivideAndSlurm()
    # make windows genome-wide with overlapping size of pattern
    intervals = DARNS2GenomicInterval(DARNS)

    task = Coverage(intervals, 40, os.path.abspath(samples["H3K4me3_K562_500k_CM"]), cpusPerTask=8)
    slurm.add_task(task)
    slurm.submit(task)

    while not task.is_ready():
        time.sleep(10)

    # collect and serialize
    cov = task.collect()

    # Add all features
    features = ["length", "space_upstream", "space_downstream", "read_count", "read_density", "n_pospeaks", "n_negpeaks"]
    DARNS_features = pd.DataFrame()
    for chrm, _ in DARNS.items():
        for i in range(len(DARNS[chrom])):
            start, end, center = DARNS[chrm][i]
            name = "_".join([chrm, str(start), str(end), str(center)])

            s = pd.Series(index=["name"] + features)
            s["name"] = name

            # measure length
            s["length"] = end - start

            # measure distance to neighbours
            if i != 0:
                s["space_upstream"] = abs(start - DARNS[chrom][i - 1][1])  # current start minus previous end
            else:
                s["space_upstream"] = None
            if i != len(DARNS[chrom]) - 1:
                s["space_downstream"] = abs(DARNS[chrom][i + 1][0] - end)  # current end minus next start
            else:
                s["space_downstream"] = None

            # get posterior prob
            sequence = genome_binary[chrom][start: end]
            s["post_prob"] = model.retrieveProbabilities(sequence)

            # get read count and density
            s["read_count"] = cov[name].sum()
            s["read_density"] = cov[name].sum()  # this should be the whole array

            # n. of positive peaks in nucleosome model
            s["n_pospeaks"] = hmmOutput[chrom][start: end].count("+")

            # n. of negative peaks in nucleosome model
            s["n_negpeaks"] = hmmOutput[chrom][start: end].count("-")

            # append
            DARNS_features = DARNS_features.append(s, ignore_index=True)

    # serialize
    pickle.dump(
        DARNS_features,
        open(os.path.join(args.results_dir, sampleName + "_DARNS.features.pickle"), "wb"),
        protocol=pickle.HIGHEST_PROTOCOL
    )

    # Train classifier on features
    # get random tenth of data to train on
    randomRows = [random.randrange(0, len(DARNS[chrom])) for _ in range(len(DARNS[chrom]) / 10)]
    train_features = DARNS_features.loc[randomRows, features]

    # get class labels for data
    labels = ["real" for _ in range(len(DARNS[chrom]) / 10)]

    lr = LinearRegression()
    lr.fit(np.array(train_features), labels, n_jobs=-1)

    # Predict for all data
    pred = lr.predict(DARNS_features.loc[: features])








    ##### OLD #####

    # Plot attributes
    assert len(DARNS) == len(DARNSP)
    for data in (DARNS, DARNSP):
        widths = Counter([darn[1] - darn[0] for darn in data.values()[0]])
        distances, midDistances = measureDARNS(data)

        plt.plot(widths.keys(), widths.values(), '-', color='orange')
        plt.plot(distances.keys(), distances.values(), '-', color='purple')
        plt.plot(midDistances.keys(), midDistances.values(), '-', color='green')
        # sns.violinplot(widths.values())
        # decide how to save

    # Get scores from DARNS
    probs = {peak: model.retrieveProbabilities(sequence) for peak, sequence in genome_binary.items()}

    # Plot predicted darns and post probs in an example region
    #
    chrom = "chr1"
    start = 5990
    end = 6000
    width = (end - start) * 1000  # length in bp

    # retrieve post. probs for selected region
    probs = model.retrieveProbabilities(OrderedDict(sorted((genome_binary[chrom][start:end]))))  # random 10kb region
    probsP = model.retrieveProbabilities(OrderedDict(sorted((genome_binaryP[chrom][start:end]))))  # random 10kb region

    # start plotting
    from matplotlib.patches import Rectangle
    colors = sns.color_palette('deep', n_colors=6, desat=0.5)
    sns.set_context(rc={"figure.figsize": (14, 6)})
    sns.plt.axhline(y=1.1, c=colors[0], alpha=0.7)
    sns.plt.xlim([1, width + 1])
    sns.plt.ylim([0, 1.2])
    sns.plt.ylabel(r'posterior probs, $\gamma_k$')
    sns.plt.xlabel(r'$k$')
    axis = sns.plt.gca()

    # viterbi predicted DARNS
    # get darns in window
    DARNSinWindow = [darn for darn in DARNS[chrom] if darn[0] >= start and darn[1] <= end]
    # for each darn in window draw a box
    for start, end in DARNSinWindow:
        axis.add_patch(Rectangle((start + 1, 1.075), end - start + 1, 0.05,
                                 facecolor=colors[0], alpha=0.7))

    # line plot of post. probs
    sns.plt.plot(range(1, width + 1), probs,  # post. probs of real data
                 c=colors[2], alpha=0.7)
    sns.plt.plot(range(1, width + 1), probsP,  # post. probs of permuted data
                 c=colors[3], alpha=0.5)
    plt.show()  # decide how to save

    # TO IMPLEMENT:
    # Export bed files of predicted DARNS for both real data and permuted
    exportBedFile(
        DARNS,
        os.path.join(args.results_dir, sampleName + ".DARNS.bed"),
        "DARNS predicted from %s" % sampleName
    )
    exportBedFile(
        DARNS,
        os.path.join(args.results_dir, sampleName + ".DARNS.permuted.bed"),
        "DARNS predicted from %s permuted reads" % sampleName
    )
    # Get overal score for darn from post. prob (or something else)

    # Get overal score for darn over permuted

    # Output bed/wig with scores!
    # Export wig files with raw correlations
    # exportWigFile(
    #     [peaks[i] for i in correlationsPos.keys()],
    #     correlationsPos.values(),
    #     len(pattern) / 2,
    #     os.path.join(args.results_dir, sampleName + ".peakCorrelationPos.wig"),
    #     sampleName + " raw absolute correlation - positive strand"
    # )
    # exportWigFile(
    #     [peaks[i] for i in correlationsNeg.keys()],
    #     correlationsNeg.values(),
    #     len(pattern) / 2,
    #     os.path.join(args.results_dir, sampleName + ".peakCorrelationNeg.wig"),
    #     sampleName + " raw absolute correlation - negative strand"
    # )

    # TODO:
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

        # Make jobs with groups of data
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
        if not hasattr(self, "output"):  # if output is already stored, just return it
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
        else:
            return self.output


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
            self.strand_wise = False
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
        if not hasattr(self, "output"):  # if output is already stored, just return it
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
        else:
            return self.output


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
        if not hasattr(self, "output"):  # if output is already stored, just return it
            if self.is_ready():
                # load all pickles into list
                if self.permissive:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
                else:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
                # if all are dicts, and their elements are dicts, sum them
                if all([type(outputs[i]) == dict for i in range(len(outputs))]):
                    output = reduce(lambda x, y: dict(x, **y), outputs)
                    if type(output) == dict:
                        self.output = output  # store output in object
                        self._rm_temps()  # delete tmp files
                        return self.output
            else:
                raise TypeError("Task is not ready yet.")
        else:
            return self.output


class CorrelatePatternBam(Task):
    """
    Task to get Correlation peaks from reads and a given pattern.
    """
    def __init__(self, data, fractions, *args, **kwargs):
        super(CorrelatePatternBam, self).__init__(data, fractions, *args, **kwargs)
        # Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1, 1000))], sep="_")
        self.name = "CorrelatePatternBam_{0}".format(now)

        # Parse
        # required argument
        if len(self.args) != 2:
            raise TypeError("You must specify two arguments: Pattern, Bam")
        self.pattern = self.args[0]
        self.bam_file = self.args[1]
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
            jobFile = files[i] + "_correlatePatternBam.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            # assemble job file
            # header
            job = self._slurmHeader(ids[i])

            # command - add abspath to python script!
            task = """\
                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python correlatePatternBam_parallel.py {0} {1} {2} {3} """.format(inputPickle, outputPickle, self.patternPickle, self.bam_file)

            task += "--step {0} ".format(self.step)
            task += "--fragment-size {0} ".format(self.fragment_size)
            if self.strand_wise:
                task += "--strand-wise "
            if self.duplicates:
                task += "--duplicates "
            if self.orientation:
                task += "--orientation "
            if self.permute:
                task += "--permute "

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
        if not hasattr(self, "output"):  # if output is already stored, just return it
            if self.is_ready():
                # load all pickles into list
                if self.permissive:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles if os.path.isfile(outputPickle)]
                else:
                    outputs = [pickle.load(open(outputPickle, 'r')) for outputPickle in self.outputPickles]
                # if all are dicts, and their elements are dicts, sum them
                if all([type(outputs[i]) == dict for i in range(len(outputs))]):
                    output = reduce(lambda x, y: dict(x, **y), outputs)
                    if type(output) == dict:
                        # Split name into chrom start end
                        output = {(key.split("_")[0], int(key.split("_")[1]), int(key.split("_")[2])): value for key, value in output.items()}
                        # Sort, make ordered dict, store output in object
                        self.output = OrderedDict(sorted(output.items(), key=lambda x: (x[0][0], x[0][1])))
                        # Delete tmp files
                        self._rm_temps()
                        return self.output
            else:
                raise TypeError("Task is not ready yet.")
        else:
            return self.output


class WinterHMM(object):
    """
    Class to model the Hidden Markov model described in the Winter 2013 paper.
    """
    def __init__(self):
        super(WinterHMM, self).__init__()

        self.mapping = {v: k for k, v in enumerate(itertools.product([1, 0, -1], repeat=2))}
        # self.mapping = {
        #     (-1, -1):   8,
        #     (-1, 0):    7,
        #     (-1, 1):    6,
        #     (0, -1):    5,
        #     (0, 0):     4,
        #     (0, 1):     3,
        #     (1, -1):    2,
        #     (1, 0):     1,
        #     (1, 1):     0
        # }

        # rewrite like this:
        # matrix = array([
        #     [-inf,-inf, -0.69314718,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf, -0.69314718,-inf,-inf, -0.69314718,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -0.10536052,-inf,-inf,-inf,-inf, -2.30258509, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -0.69314718, -0.69314718,-inf,-inf, -inf,-inf],
        #     [-inf,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,  0.,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,  0.,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf,-inf, -0.69314718,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,  0.,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,  0.,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf, -0.7985077 ,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -0.7985077 , -2.30258509, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,  0.,-inf,-inf,-inf, -inf,-inf],
        #     [-0.69314718,-inf, -0.69314718,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf],
        #     [-inf,-inf, -2.30258509,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -2.30258509,-inf,-inf, -0.22314355, -inf,-inf],
        #     [-inf,-inf, -2.30258509,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -2.30258509,-inf,-inf, -0.22314355, -inf,-inf],
        #     [-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf, -inf,-inf]
        # ])

        # distributions = [NormalDistribution(1, .5), NormalDistribution(5, 2)]
        # starts = [ 1., 0. ]
        # ends = [ .1., .1 ]
        # state_names= [ "A", "B" ]

        # model = Model.from_matrix( matrix, distributions, starts, ends,
        #         state_names, name="test_model" )

        background = yahmm.DiscreteDistribution({
            0: 0.001, 1: 0.001, 2: 0.001,
            3: 0.001, 4: 0.99, 5: 0.002,
            6: 0.001, 7: 0.002, 8: 0.001
        })
        pos = yahmm.DiscreteDistribution({
            0: 0.33, 1: 0.33, 2: 0.33,
            3: 0.001, 4: 0.001, 5: 0.001,
            6: 0.001, 7: 0.001, 8: 0.001
        })
        neg = yahmm.DiscreteDistribution({
            0: 0.33, 1: 0.001, 2: 0.001,
            3: 0.33, 4: 0.001, 5: 0.001,
            6: 0.33, 7: 0.001, 8: 0.001
        })
        zero = yahmm.DiscreteDistribution({
            0: 0.001, 1: 0.2, 2: 0.001,
            3: 0.2, 4: 0.2, 5: 0.2,
            6: 0.001, 7: 0.2, 8: 0.001
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

        for state in [b, minus, in1, in2, in3, plus, in4, in5, in6, in7, in8, in9, in10, in11]:
            self.model.add_state(state)

        self.model.add_transition(self.model.start, b, 0.8)       # start in background
        self.model.add_transition(self.model.start, minus, 0.1)   # start in minus
        self.model.add_transition(self.model.start, plus, 0.1)    # start in plus

        self.model.add_transition(b, b, 0.8)                      # stay in background
        self.model.add_transition(b, minus, 0.1)                 # enter in minus
        self.model.add_transition(b, plus, 0.1)                  # enter in plus

        self.model.add_transition(minus, b, 0.1)                  # can exit from plus
        self.model.add_transition(minus, in1, 0.45)
        self.model.add_transition(minus, plus, 0.45)

        self.model.add_transition(in1, in2, 0.5)                  # states 1-3 can go directly to plus or follow a string 1 to 4
        self.model.add_transition(in1, plus, 0.5)
        self.model.add_transition(in2, in3, 0.5)
        self.model.add_transition(in2, plus, 0.5)
        self.model.add_transition(in3, in4, 0.5)
        self.model.add_transition(in3, plus, 0.5)

        self.model.add_transition(plus, b, 0.1)                   # can exit from plus
        self.model.add_transition(plus, in4, 0.9)
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
        Retrieve the posterior probabilities of being in a "nucleosome" state for a sequence of observations.
        """
        trans, ems = self.model.forward_backward(observations)
        ems = np.exp(ems)
        self.probs = ems / np.sum(ems, axis=1)[:, np.newaxis]  # probs of all states
        background_state = 0
        prob_nucleosome = 1 - self.probs[:, background_state]  # prob of all states but background
        return prob_nucleosome

    def retrieveBackgroundProbabilities(self, observations):
        """
        Retrieve the posterior probabilities of being in a "nucleosome" state for a sequence of observations.
        """
        trans, ems = self.model.forward_backward(observations)
        ems = np.exp(ems)
        self.probs = ems / np.sum(ems, axis=1)[:, np.newaxis]  # probs of all states
        background_state = 0
        prob_background = self.probs[:, background_state]  # prob of all states but background
        return prob_background


def makeGenomeWindows(windowWidth, genome, step=None):
    """
    Generate windows genome-wide for a given genome with width=windowWidth and
    return dictionary of HTSeq.GenomicInterval objects.

    windowWidth=int.
    genome=str.
    """
    if step is None:
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
    if step is None:
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


def DARNS2GenomicInterval(DARNS):
    """
    """
    intervals = OrderedDict()
    for chrom, values in DARNS.items():
        for start, end, center in values:
            name = "_".join([chrom, str(start), str(end), str(center)])
            intervals[name] = HTSeq.GenomicInterval(chrom, start, end)

    return intervals


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


def concatenateBinary(binaryDict, patternLength, debug=False):
    """
    Concatenates string of binary signals from several consecutive windows.
    Returns dict of chr:binary.

    binaryDict=dict - name: 1D numpy.array of ~binary (1,0,-1) signals.
    patternLength=int - length of pattern used to create binary correlation peaks.
    """
    # TODO: make sure windows are sorted alphanumerically
    items = binaryDict.items()
    binaries = OrderedDict()
    prev_chrom = ""
    prev_end = 1

    for i in xrange(len(items)):
        window, sequence = items[i]
        chrom, start, end = (str(window[0]), int(window[1]), int(window[2]))
        if debug:
            name = "_".join([str(j) for j in window])
            print name, (chrom, start, end)
        if not chrom == prev_chrom:  # new chromossome
            binaries[chrom] = sequence
        else:  # same chromossome as before
            if window == items[-1][0]:  # if last window just append remaining
                if debug:
                    print("Last of all")
                binaries[chrom] = np.hstack((binaries[chrom], sequence))
            elif items[i + 1][0][0] != chrom:  # if last window in chromossome, just append remaining
                if debug:
                    print("Last of chromossome")
                binaries[chrom] = np.hstack((binaries[chrom], sequence))
            else:  # not last window
                if prev_end - patternLength == start:  # windows are continuous
                    if len(sequence) == (end - start) - patternLength + 1:
                        binaries[chrom] = np.hstack((binaries[chrom], sequence))
                    elif len(sequence) > (end - start) - patternLength + 1:
                        if debug:
                            print(name, len(sequence), (end - start) - patternLength + 1)
                        raise ValueError("Sequence is bigger than its coordinates.")
                    elif len(sequence) < (end - start) - patternLength + 1:
                        if debug:
                            print(name, len(sequence), (end - start) - patternLength + 1)
                        raise ValueError("Sequence is shorter than its coordinates.")
                else:
                    if debug:
                        print(name, prev_end, start)
                    raise ValueError("Windows are not continuous or some are missing.")

        prev_chrom = chrom
        prev_end = end
    return binaries


def splitAndSortDict(binary):
    """
    Splits dictionary keys in chrom, start, end; sorts dict based on that and returns collections.OrderedDict object.

    binary=dict - keys must be strings that when separated give a chrom (string), start and end (both ints).
    """
    b = {(key.split("_")[0], int(key.split("_")[1]), int(key.split("_")[2])): value for key, value in binary.items()}
    return OrderedDict(sorted(b.items(), key=lambda x: (x[0][0], x[0][1])))


def getDARNS(sequence, debug=False):
    """
    Returns list of tuples with start and end positions (0-based) of DARNS.
    DARNS are stretches of sequences which belong to the "nucleosome" states (see WinterHMM object definition).

    sequence=str.
    """
    outside = ["B"]
    DARNS = list()
    DARN = False
    l = len(sequence)

    # Detect DARNS start and end
    for i in xrange(l):
        if not i == l - 1:  # not last element
            if debug: print(i, "Not last element")
            # Entering DARN
            if not DARN and (sequence[i] not in outside) and i != len(sequence):
                if debug: print(i, "Entered DARN")
                DARN = True
                DARN_start = i
                # Handle 1-element DARNS
                if sequence[i + 1] in outside:
                    if debug: print(i, "1bp DARN")
                    #DARNS.append((DARN_start, i + 1))
                    DARN = False
            # DARN ends at next element
            elif DARN and (sequence[i] not in outside) and (sequence[i + 1] in outside):
                if debug: print(i, "Exited DARN")
                DARN = False
                DARN_end = i
                DARNS.append((DARN_start, DARN_end))
        elif i == l - 1 and DARN:  # last element
            if debug: print(i, "Last element")
            if DARN:  # finish DARN at end
                if debug: print(i, "Finishing DARN in last element")
                DARNS.append((DARN_start, i))

    # Find DARNS middle point
    DARNS_middle = list()
    for start, end in DARNS:
        seq = sequence[start: end + 1]
        middle = int(((end + 1) - start) / 2)
        negative_peaks = [i for i in xrange(len(seq)) if seq[i] == '-']

        if len(negative_peaks) is not 0:
            # find middle most negative peak
            peak = start + min(negative_peaks, key=lambda x: abs(x - middle))
        else:
            peak = start + middle

        DARNS_middle.append((start, end, peak))

    return DARNS_middle


def measureDARNS(DARNS):
    """
    Computes distributions of DARNS attributes: width, interdistance, distance between midpoints.

    tuples=list - list of tuples: (start,end).
    """
    start, end = range(2)
    distances = Counter()
    midDistances = Counter()

    try:
        for chrom, darns in DARNS.items():
            for d1, d2 in itertools.combinations(darns, 2):
                # distance end-to-end
                if d1[end] <= d2[start]:
                    distances[abs(d2[start] - d1[end])] += 1
                else:
                    distances[d1[start] - d2[end]] += 1

                # distance midpoint-midpoint
                midDistances[abs(((d1[end] + d1[start]) / 2) - ((d2[end] + d2[start]) / 2))] += 1
        return (distances, midDistances)
    except KeyboardInterrupt:
        return (distances, midDistances)


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


def exportBedFile(intervals, filename, trackname):
    """
    Exports a bed file track from dict with genomic positions.

    intervals=dict - chrom:(start, end) dict.
    filename=str.
    trackname=str.
    strand=str.
    """
    with open(filename, 'w') as handle:
        header = 'track name="{0}" description="{0}" visibility=pack autoScale=off colorByStrand="255,0,0 0,0,255"\n'.format(trackname)
        handle.write(header)
        for chrom, (start, end) in intervals.items():
            name = "DARN_{0}_{1}_{2}".format(chrom, start, end)
            entry = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                chrom, start, end,
                name,  # name
                1,  # score
                "."  # strand
            )
            handle.write(entry)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(0)
