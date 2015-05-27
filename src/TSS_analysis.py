#!/usr/env python
#############################################################################################
#
# This code was used to produce the plots in the ChIPmentation paper (Schmidl, et al. 2015).
#
# It produces plots of average signal profiles around TSSs,
# clusters of TSSs based on various signals,
# and heatmaps of the same
#
#############################################################################################

import os
import re
from collections import OrderedDict
import HTSeq
import pybedtools
import numpy as np
import pandas as pd
import cPickle as pickle
from sklearn.cluster import k_means
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns  # changes plt style (+ full plotting library)


# Force matplotlib to not use any Xwindows backend.

# Set seaborn style
sns.set_style("whitegrid")


def loadBed(filename):
    """
    Parses bed file and returns dict of featureName:HTSeq.GenomicInterval objects.
    filename - string.
    """
    from warnings import warn
    warn("Function is deprecated!")
    from HTSeq import GenomicInterval

    features = OrderedDict()
    for line in open(filename):
        fields = line.split("\t")
        feature = GenomicInterval(
            fields[0],                      # chrom
            int(fields[1]),                 # start
            int(fields[2]),                 # end
            fields[5]                       # strand
        )
        features[fields[4]] = feature       # append with name
    return features


def bedToolsInterval2GenomicInterval(bedtool):
    """
    Given a pybedtools.BedTool object returns, dictionary of HTSeq.GenomicInterval objects.
    bedtool - a pybedtools.BedTool with intervals.
    """
    intervals = OrderedDict()
    for iv in bedtool:
        intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, iv.strand)
    return intervals


def coverage(bam, intervals, fragmentsize, orientation=True, duplicates=True, strand_specific=False, switchChromsNames=False):
    """
    Gets read coverage in bed regions.
    Returns dict of regionName:numpy.array if strand_specific=False, A dict of "+" and "-" keys with regionName:numpy.array.
    bam - HTSeq.BAM_Reader object. Must be sorted and indexed with .bai file!
    intervals - dict with HTSeq.GenomicInterval objects as values.
    fragmentsize - integer.
    stranded - boolean.
    duplicates - boolean.
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
            i += 1
            continue

        # Replace chromosome reference 1 -> chr1 if not chr
        if switchChromsNames:
            feature.chrom = re.sub("chr", "", feature.chrom)

        # Fetch alignments in feature window
        for aln in bam[feature]:
            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue
            # check it's aligned
            if not aln.aligned:
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
            if not strand_specific:
                profile[start_in_window: end_in_window] += 1
            else:
                if aln.iv.strand == "+":
                    profile[0][start_in_window: end_in_window] += 1
                else:
                    profile[1][start_in_window: end_in_window] += 1

        # append feature profile to dict
        cov[name] = profile
        i += 1
    return cov


def coverageSum(bam, intervals, fragmentsize, orientation=True, duplicates=True, switchChromsNames=False):
    """
    Gets read coverage in bed regions.
    Returns dict of regionName:numpy.array if strand_specific=False, A dict of "+" and "-" keys with regionName:numpy.array.
    bam - HTSeq.BAM_Reader object. Must be sorted and indexed with .bai file!
    intervals - dict with HTSeq.GenomicInterval objects as values.
    fragmentsize - integer.
    stranded - boolean.
    duplicates - boolean.
    """
    # Loop through TSSs, get coverage, append to dict
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX']

    # Initialize empty array for this feature
    profile = np.zeros(intervals.items()[0][1].length, dtype=np.float64)

    n = len(intervals)
    i = 0
    for name, feature in intervals.iteritems():
        if i % 1000 == 0:
            print(n - i)

        # Check if feature is in bam index
        if feature.chrom not in chroms or feature.chrom == "chrM":
            i += 1
            continue

        # Replace chromosome reference 1 -> chr1 if not chr
        if switchChromsNames:
            feature.chrom = re.sub("chr", "", feature.chrom)

        # Fetch alignments in feature window
        for aln in bam[feature]:
            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue
            # check it's aligned
            if not aln.aligned:
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
            profile[start_in_window: end_in_window] += 1
        i += 1
    return profile


def plotHeatmap(df, filename):
    """
    Plot heatmap for data in dataframe using matplotlib.pyplot.imshow.
    df - pandas.DataFrame with numerical values.
    filename - string with path and filename to save pdf to.
    """
    plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    ax = plt.imshow(df, interpolation='nearest', aspect='auto', vmin=0, vmax=0.5).get_axes()
    ax.grid('off')
    plt.colorbar(orientation="vertical")
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.close()


def exportToJavaTreeView(df, filename):
    """
    Export cdt file of cluster to view in JavaTreeView.
    df - pandas.DataFrame object with numeric data.
    filename - string.
    """
    cols = ["X" + str(x) for x in df.columns]
    df.columns = cols
    df["X"] = df.index
    df["NAME"] = df.index
    df["GWEIGHT"] = 1
    df = df[["X", "NAME", "GWEIGHT"] + cols]
    df.to_csv(filename, sep="\t", index=False)


def normalize(x):
    """
    Implements standardization of data.
    x - numerical iterable (e.g. list, numpy.array).
    """
    return (x - min(x)) / (max(x) - min(x))


def smooth(x, window_len=8, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y


# Define variables
projectRoot = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/"
bamsDir = os.path.join(projectRoot, "data", "mapped/")
resultsDir = os.path.join(projectRoot, "results")
plotsDir = os.path.join(resultsDir, "footprints")
DNase = os.path.join(bamsDir, "wgEncodeUwDnaseK562Aln.merged.bam")
MNase = os.path.join(bamsDir, "wgEncodeSydhNsomeK562Aln.merged.bam")
NexteraBackground = "/home/arendeiro/PGA1Nextera/PGA_0001_Nextera-2.bam"

# Get samples
samples = pd.read_csv(os.path.abspath(projectRoot + "chipmentation.replicates.annotation_sheet.csv"))

# replace bam path with local
fhPath = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/"
samples.loc[:, "filePath"] = samples["filePath"].apply(lambda x: re.sub(fhPath, bamsDir, x))

# subset samples
sampleSubset = samples[
    samples["sampleName"].str.contains(
        "K562_500K_CM_H3K4ME3_nan_nan_0_0_hg19|K562_50K_ATAC_nan_nan_nan_0_0_hg19"
    )
].reset_index(drop=True)

sampleSubset = sampleSubset.append(pd.Series(data=["DNase", DNase], index=["sampleName", "filePath"]), ignore_index=True)
sampleSubset = sampleSubset.append(pd.Series(data=["NexteraBackground", NexteraBackground], index=["sampleName", "filePath"]), ignore_index=True)

sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

bedFilePath = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/annotation/hg19.cage_peak_coord_robust.TATA_Annotated.bed"
genome = "hg19"
windowRange = (-80, 95)
fragmentsize = 1
duplicates = True
n_clusters = 5

windowWidth = abs(windowRange[0]) + abs(windowRange[1])

# Load TSSs from bed file, transform by window width
tsss = pybedtools.BedTool(bedFilePath).slop(genome=genome, b=windowWidth / 2)
tsss = bedToolsInterval2GenomicInterval(tsss)

# Filter tsss near chrm borders
for name, interval in tsss.iteritems():
    if interval.length < windowWidth:
        tsss.pop(name)

# Loop through all signals, compute coverage in bed regions,
# save dicts with coverage and average profiles
rawSignals = OrderedDict()
aveSignals = OrderedDict()
for i in range(len(sampleSubset)):
    name = sampleSubset['sampleName'][i]

    if name not in aveSignals.keys():
        print(name)

        # Load bam
        bamfile = HTSeq.BAM_Reader(sampleSubset['filePath'][i])

        # Get dataframe of signal coverage in bed regions, append to dict
        cov = coverageSum(bamfile, tsss, fragmentsize)

        # Make multiindex dataframe
        levels = [cov.keys(), ["+", "-"]]
        labels = [[y for x in range(len(cov)) for y in [x, x]], [y for x in range(len(cov.keys())) for y in (0, 1)]]
        index = pd.MultiIndex(labels=labels, levels=levels, names=["tss", "strand"])
        df = pd.DataFrame(np.vstack(cov.values()), index=index, columns=range(windowRange[0], windowRange[1]))

        # For strand_specific=False
        # cov = coverage(bamfile, tsss, fragmentsize)
        # df = pd.DataFrame(cov).T
        # df.columns = range(windowRange[0], windowRange[1])

        # append to dict
        rawSignals[name] = df
        pickle.dump(rawSignals, open(os.path.join(resultsDir, "TSS_rawSignals.2.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        # Get average profiles and append to dict
        ave = {
            "name": name,
            "average": df.apply(np.mean, axis=0),                              # both strands
            "positive": df.ix[range(0, len(df), 2)].apply(np.mean, axis=0),    # positive strand
            "negative": df.ix[range(1, len(df), 2)].apply(np.mean, axis=0)     # negative strand
        }
        aveSignals[name] = pd.DataFrame(ave)
        pickle.dump(aveSignals, open(os.path.join(plotsDir, "TSS_aveSignals.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

# Add coverage of dyads in TSSs
name = "dyads"

tsss = pybedtools.BedTool(bedFilePath).slop(genome=genome, b=windowWidth / 2)
dyads = pybedtools.BedTool(os.path.join(projectRoot, "data", "dyads", "dyads.bed"))

cov = dyads.coverage(tsss, d=True)
cov = cov.to_dataframe(
    names=['chrom', 'start', 'end', 'name', '1', '2', '3', '4', '5', '6', '7', 'bp', "count"]
)
cov = cov[["bp", "count"]]

# Get average profiles and append to dict
cov = cov.groupby(["bp"]).apply(np.mean)
ave = {
    "name": name,
    "average": cov["count"],  # both strands
}
aveSignals[name] = pd.DataFrame(ave)
pickle.dump(aveSignals, open(os.path.join(resultsDir, "TSS_aveSignals.2.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

# convert dict to dataframe
aveSignals = pickle.load(open(os.path.join(resultsDir, "TSS_aveSignals.2.pickle"), "r"))
aveSignals = pd.concat(aveSignals, axis=0)
aveSignals["x"] = [i[1] for i in aveSignals.index]
aveSignals = aveSignals.sort(["name"]).reset_index(drop=True)
aveSignals.drop_duplicates(inplace=True)
aveSignals["type"] = "raw"
aveSignals.reset_index(drop=True, inplace=True)

# Normalize by controls
df = aveSignals.copy()
df2 = pd.DataFrame()

ctrls = [
    "K562_10K_CM_IGG_nan_nan_0_0_hg19",
    "K562_500K_CM_IGG_nan_nan_1_0_hg19",
    "K562_50K_ATAC_nan_nan_nan_0_0_hg19",
    "K562_500K_ATAC_H3K4ME3_nan_01ULTN5_PE_1_1_hg19",
    "K562_500K_ATAC_INPUT_nan_01ULTN5_PE_1_1_hg19"
]

for sample in df.name.unique():
    for strand in ["average", "positive", "negative"]:
        # scale and smooth raw data
        treat = df[(df["name"] == sample)][strand]

        tmp = pd.DataFrame()
        scaled = smooth((treat - treat.min()) / (treat.max() - treat.min()))

        tmp[strand] = scaled
        tmp['name'] = sample
        tmp['x'] = np.array(range(len(tmp))) - (len(tmp) / 2)
        tmp['type'] = "raw"

        df2 = df2.append(tmp, ignore_index=True)

        for control in ctrls:
            treat = df[
                (df["name"] == sample)
            ][strand]

            ctrl = df[
                (df["name"] == control)
            ][strand]

            # standardize: 0-1
            treat = smooth((treat - treat.min()) / (treat.max() - treat.min()))
            ctrl = smooth((ctrl - ctrl.min()) / (ctrl.max() - ctrl.min()))

            # normalize by input
            norm = np.array(treat) - np.array(ctrl)

            tmp = pd.DataFrame()
            tmp[strand] = norm
            tmp['name'] = sample
            tmp['x'] = np.array(range(len(tmp))) - (len(tmp) / 2)
            tmp['type'] = control

            df2 = df2.append(tmp, ignore_index=True)
aveSignals = df2

# Plot average profiles
# Grid plots
# raw
sub = aveSignals[aveSignals['type'] == "raw"]
grid = sns.FacetGrid(sub, col="name", hue="name", sharey=False, col_wrap=4, size=4, aspect=2)
grid.map(plt.plot, "x", "positive")
grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
grid.add_legend()
plt.savefig(os.path.join(plotsDir, "all.tss2.raw.pdf"), bbox_inches='tight')

# norm - all in one
grid = sns.FacetGrid(aveSignals, col="name", hue="type", sharey=False, col_wrap=4, size=4, aspect=2)
grid.map(plt.plot, "x", "positive")
grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
grid.add_legend()
plt.savefig(os.path.join(plotsDir, "all.tss2.norm.pdf"), bbox_inches='tight')

# raw
for sample in aveSignals['name'].unique():
    for signal in aveSignals[aveSignals['name'] == sample]["type"].unique():
        sub = aveSignals[
            (aveSignals['name'] == sample) &
            (aveSignals['type'] == signal)
        ]
        plt.plot(sub['x'], sub['positive'], label=signal)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(os.path.join(plotsDir, sample + ".tss2.pdf"), bbox_inches='tight')
    plt.close()

# overlaid plots
# dnase and H3K4ME3
sub = aveSignals[
    (aveSignals['type'] == "raw") &
    (aveSignals['x'] >= -80) &
    (aveSignals['x'] <= 80) &
    (aveSignals['name'].str.contains(
        "K562_500K_CM_H3K4ME3_nan_nan_0_0_hg19|" +
        "DNase"))
]
sub.drop(["positive", "negative", "type"], axis=1, inplace=True)
sub.dropna(inplace=True)

for name in sub.name.unique():
    x = sub["x"].unique()
    y = normalize(sub[sub['name'] == name]["average"])
    plt.plot(x, y)
plt.savefig(os.path.join(plotsDir, "all.tss2.raw.overlaid.open.pdf"), bbox_inches='tight')
plt.close("all")

# mnase and H3K27ME3
sub = aveSignals[
    (aveSignals['type'] == "raw") &
    (aveSignals['x'] >= -80) &
    (aveSignals['x'] <= 80) &
    (aveSignals['name'].str.contains(
        "K562_500K_CM_H3K27ME3_nan_nan_0_0_hg19|" +
        "K562_10M_CM_H3K4ME1_nan_nan_0_0_hg19|" +
        "MNase"))
]
sub.drop(["positive", "negative", "type"], axis=1, inplace=True)
sub.dropna(inplace=True)

for name in sub.name.unique():
    x = sub["x"].unique()
    y = normalize(sub[sub['name'] == name]["average"])
    plt.plot(x, y, label=name)
plt.legend()
plt.savefig(os.path.join(plotsDir, "all.tss2.raw.overlaid.closed+1K4.pdf"), bbox_inches='tight')
plt.close("all")


# HEATMAPS
rawSignals = pd.read_pickle(os.path.join(resultsDir, "TSS_rawSignals.2.pickle"))

# Loop through raw signals, normalize, k-means cluster, save pickle, plot heatmap, export cdt
signals = ["K562_500K_CHIP_H3K4ME3_nan_nan_0_0_hg19",
           "DNase",
           "MNase",
           "K562_50K_ATAC_nan_nan_nan_0_0_hg19",
           "K562_500K_CM_H3K4ME3_nan_nan_0_0_hg19"]

for name in signals:
    print("Clustering based on %s" % name)
    exportName = "{0}.tssSignal_{1}bp.kmeans_{2}k".format(name, str(windowWidth), n_clusters)

    df = rawSignals[name]

    # join strand name (plus + minus)
    df = df.sum(axis=0, level="tss")

    # scale row signal to 0:1 (standardization)
    dfNorm = df.apply(normalize, axis=1)
    dfNorm.replace(["inf", "NaN"], 0, inplace=True)

    # cluster
    clust = k_means(
        dfNorm,
        n_clusters,
        n_init=25,
        max_iter=10000,
        n_jobs=-1
    )  # returns centroid, label, inertia

    # keep dataframe with assigned cluster and row index
    clustOrder = pd.DataFrame(
        clust[1],
        columns=["cluster"],
        index=dfNorm.index.values
    )

    # save object
    pickle.dump(
        clust,
        open(os.path.join(plotsDir, "..", exportName + ".pickle"), "wb"),
        protocol=pickle.HIGHEST_PROTOCOL
    )

    print("Exporting heatmaps for %s" % name)

    # now sort by clust order
    dfNorm = pd.merge(dfNorm, clustOrder, on=None, left_index=True, right_index=True)  # get common rows (should be all)
    dfNorm.replace(["inf", "NaN"], 0, inplace=True)
    dfNorm.sort_index(by="cluster", axis=0, inplace=True)
    dfNorm.drop("cluster", axis=1, inplace=True)

    # Plot heatmap
    # data = Data([Heatmap(z=np.array(dfNorm), colorscale='Portland')])
    # plotly.image.save_as(data, os.path.join(plotsDir, exportName + ".plotly.pdf"))
    # plotHeatmap(dfNorm, os.path.join(plotsDir, exportName + "." + name2 + "_signal.matplotlib.pdf"))

    # Export as cdt
    # exportToJavaTreeView(dfNorm.copy(), os.path.join(plotsDir, exportName + "_signal.cdt"))

    # Sort all signals by dataframe by cluster order
    for name2 in signals:
        print("Exporting heatmaps for %s" % name2)
        df = rawSignals[name2]
        df = df.xs('+', level="strand") + df.xs('-', level="strand")

        # scale row signal to 0:1 (standardization)
        dfNorm = df.apply(lambda x: (x - min(x)) / (max(x) - min(x)), axis=1)
        dfNorm.replace("inf", 0, inplace=True)

        # now sort by clust order
        dfNorm = pd.merge(dfNorm, clustOrder, on=None, left_index=True, right_index=True)  # get common rows (should be all)
        dfNorm.replace(["inf", "NaN"], 0, inplace=True)
        dfNorm.sort_index(by="cluster", axis=0, inplace=True)
        dfNorm.drop("cluster", axis=1, inplace=True)

        # Plot heatmap
        # data = Data([Heatmap(z=np.array(dfNorm), colorscale='Portland')])
        # plotly.image.save_as(data, os.path.join(plotsDir, exportName + ".plotly.pdf"))
        # plotHeatmap(dfNorm, os.path.join(plotsDir, exportName + "." + name2 + "_signal.matplotlib.pdf"))

        # Export as cdt
        exportToJavaTreeView(dfNorm.copy(), os.path.join(plotsDir, exportName + "." + name2 + "_signal.cdt"))
