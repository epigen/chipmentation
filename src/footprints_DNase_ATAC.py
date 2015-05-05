#!/usr/env python
#############################################################################################
#
# This code was used to produce the plots in the ChIPmentation paper (Schmidl, et al. 2015).
#
# It produces plots of average signal profiles around TF motifs in peaks
# and heatmaps of the same
#
#############################################################################################

import os
from collections import OrderedDict
import re
import HTSeq
import pybedtools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  # changes plt style (+ full plotting library)
import cPickle as pickle

sns.set_style("whitegrid")


def smooth(x, window_len=11, window='hanning'):
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


def savePandas(fname, data):
    '''Save DataFrame or Series

    Parameters
    ----------
    fname : str
        filename to use
    data: Pandas DataFrame or Series
    '''
    np.save(open(fname, 'w'), data)
    if len(data.shape) == 2:
        meta = data.index, data.columns
    elif len(data.shape) == 1:
        meta = (data.index,)
    else:
        raise ValueError('save_pandas: Cannot save this type')
    s = pickle.dumps(meta)
    s = s.encode('string_escape')
    with open(fname, 'a') as f:
        f.seek(0, 2)
        f.write(s)


def loadPandas(fname, mmap_mode='r'):
    '''Load DataFrame or Series

    Parameters
    ----------
    fname : str
        filename
    mmap_mode : str, optional
        Same as np.load option
    '''
    values = np.load(fname, mmap_mode=mmap_mode)
    with open(fname) as f:
        np.lib.format.read_magic(f)
        np.lib.format.read_array_header_1_0(f)
        f.seek(values.dtype.alignment * values.size, 1)
        meta = pickle.loads(f.readline().decode('string_escape'))
    if len(meta) == 2:
        return pd.DataFrame(values, index=meta[0], columns=meta[1])
    elif len(meta) == 1:
        return pd.Series(values, index=meta[0]).copy()


def coverage(bam, intervals, fragmentsize, orientation=True, duplicates=True, strand_specific=False):
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


def bedToolsInterval2GenomicInterval(bedtool):
    """
    Given a pybedtools.BedTool object, returns dictionary of HTSeq.GenomicInterval objects.
    """
    intervals = OrderedDict()
    for iv in bedtool:
        if iv.strand == "+" or iv.strand == 0 or iv.strand == str(0):
            intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "+")
        elif iv.strand == "-" or iv.strand == 0 or iv.strand == str(1):
            intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "-")
        else:
            intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end)
    return intervals


def colourPerFactor(name):
    name = str(name.upper())
    print(name)
    if "H3K4ME3" in name:
        return "#009e73"
    elif "H3K4ME1" in name:
        return "#e69f00"
    elif "H3K27AC" in name:
        return "#D55E29"
    elif "H3K27ME3" in name:
        return "#0072b2"
    elif "H3K36ME3" in name:
        return "#9e2400"
    elif "CTCF" in name:
        return "#534202"
    elif "PU1" in name:
        return "#6E022C"
    elif "GATA1" in name:
        return "#9b0505"
    elif "GATA2" in name:
        return "#510303"
    elif "REST" in name:
        return "#25026d"
    elif "CJUN" in name:
        return "#2a0351"
    elif "FLI1" in name:
        return "#515103"
    elif "IGG" in name:
        return "#d3d3d3"
    elif "INPUT" in name:
        return "#d3d3d3"
    elif "DNASE" in name:
        return "#005212"
    elif "MNASE" in name:
        return "#00523b"
    elif "ATAC" in name:
        return "#001f07"
    else:
        raise ValueError


def getMNaseProfile(chrom, start, end):
    command = "bigWigToWig -chrom={0} -start={1} -end={2} ".format(chrom, start, end)
    command += "/home/arendeiro/encode_mnase_k562/wgEncodeSydhNsomeK562Sig.bigWig tmp.wig; "
    command += "tail -n +2 tmp.wig | cut -f 4 > tmp2.wig"
    os.system(command)

    try:
        df = pd.read_csv("tmp2.wig", sep="\t", header=-1)
        return df[0].tolist()
    except ValueError:
        return [0] * (end - start)


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


def callFootprints(cuts, annot):
    """
    Call footprints.
    Requires dataframe with cuts and dataframe with annotation (2> cols).
    """
    import rpy2.robjects as robj  # for ggplot in R
    import rpy2.robjects.pandas2ri  # for R dataframe conversion

    # Plot with R
    footprint = robj.r("""
    library(CENTIPEDE)

    function(cuts, annot) {
        centFit <- fitCentipede(
            Xlist = list(as.matrix(cuts)),
            Y = as.matrix(annot)
        )
        # imageCutSites(cuts[order(centFit$PostPr),][c(1:100, (dim(cuts)[1]-100):(dim(cuts)[1])),])
        # plotProfile(centFit$LambdaParList[[1]],Mlen=2)
        return(centFit$PostPr)
    }

    """)

    # convert the pandas dataframe to an R dataframe
    robj.pandas2ri.activate()
    cuts_R = robj.conversion.py2ri(cuts)
    annot_R = robj.conversion.py2ri(annot)

    # run the plot function on the dataframe
    return np.ndarray.flatten(robj.conversion.ri2py(footprint(cuts_R, annot_R)))


# Define variables
# projectRoot = "/projects/chipmentation/"
projectRoot = "/media/afr/cemm-backup/chipmentation/"
# projectRoot = "/home/arendeiro/chipmentation/"
bamsDir = os.path.join(projectRoot, "data", "mapped/")
peaksDir = os.path.join(projectRoot, "data", "peaks/")
resultsDir = os.path.join(projectRoot, "results")
plotsDir = os.path.join(resultsDir, "plots")
CM = os.path.join(bamsDir, "CM_H3K4ME1-H3K27AC.bam")
DNase = os.path.join(bamsDir, "wgEncodeUwDnaseK562Aln.merged.bam")
ATAC = os.path.join(bamsDir, "K562_50K_ATAC_nan_nan_nan_0_0_hg19.trimmed.bowtie2.shifted.dups.bam")

# Get samples
samples = {
    "CM_CTCF": (
        [bamsDir + "K562_10M_CM_CTCF_nan_nan_0_0_hg19.trimmed.bowtie2.shifted.dups.bam", CM],
        peaksDir + "K562_10M_CM_CTCF_nan_nan_0_0_hg19/K562_10M_CM_CTCF_nan_nan_0_0_hg19_peaks.motifCentered.bed"
    ),
    "CM_GATA1": (
        [bamsDir + "K562_10M_CM_GATA1_nan_nan_0_0_hg19.trimmed.bowtie2.shifted.dups.bam", CM],
        peaksDir + "K562_10M_CM_GATA1_nan_nan_0_0_hg19/K562_10M_CM_GATA1_nan_nan_0_0_hg19_peaks.motifCentered.bed"
    ),
    "CM_PU1": (
        [bamsDir + "K562_10M_CM_PU1_nan_nan_0_0_hg19.trimmed.bowtie2.shifted.dups.bam", CM],
        peaksDir + "K562_10M_CM_PU1_nan_nan_0_0_hg19/K562_10M_CM_PU1_nan_nan_0_0_hg19_peaks.motifCentered.bed"
    ),
    "CM_REST": (
        [bamsDir + "K562_10M_CM_REST_nan_nan_0_0_hg19.trimmed.bowtie2.shifted.dups.bam", CM],
        peaksDir + "K562_10M_CM_REST_nan_nan_0_0_hg19/K562_10M_CM_REST_nan_nan_0_0_hg19_peaks.motifCentered.bed"
    ),
    "ATAC_PU1": (
        [ATAC], peaksDir + "K562_50K_ATAC_nan_nan_nan_0_0_hg19/K562_50K_ATAC_nan_nan_nan_0_0_hg19_peaks.filtered.PU1-motifCentered.bed",
    ),
    "ATAC_GATA1": (
        [ATAC], peaksDir + "K562_50K_ATAC_nan_nan_nan_0_0_hg19/K562_50K_ATAC_nan_nan_nan_0_0_hg19_peaks.filtered.GATA1-motifCentered.bed"
    ),
    "ATAC_REST": (
        [ATAC], peaksDir + "K562_50K_ATAC_nan_nan_nan_0_0_hg19/K562_50K_ATAC_nan_nan_nan_0_0_hg19_peaks.filtered.REST-motifCentered.bed"
    ),
    "ATAC_CTCF": (
        [ATAC], peaksDir + "K562_50K_ATAC_nan_nan_nan_0_0_hg19/K562_50K_ATAC_nan_nan_nan_0_0_hg19_peaks.filtered.CTCF-motifCentered.bed"
    ),
    "DNase_PU1": (
        [DNase], peaksDir + "wgEncodeUwDnaseK562Aln.merged/wgEncodeUwDnaseK562Aln.merged_peaks.filtered.PU1-motifCentered.bed",
    ),
    "DNase_GATA1": (
        [DNase], peaksDir + "wgEncodeUwDnaseK562Aln.merged/wgEncodeUwDnaseK562Aln.merged_peaks.filtered.GATA1-motifCentered.bed",
    ),
    "DNase_REST": (
        [DNase], peaksDir + "wgEncodeUwDnaseK562Aln.merged/wgEncodeUwDnaseK562Aln.merged_peaks.filtered.REST-motifCentered.bed",
    ),
    "DNase_CTCF": (
        [DNase], peaksDir + "wgEncodeUwDnaseK562Aln.merged/wgEncodeUwDnaseK562Aln.merged_peaks.filtered.CTCF-motifCentered.bed"
    )
}

fragmentsize = 1
duplicates = False
n_clusters = 5

gapsRepeats = pybedtools.BedTool(os.path.join("/media/afr/cemm-backup", "reference/hg19/hg19_gapsRepeats.bed"))

foots = dict()

# Loop through all samples, compute coverage in peak regions centered on motifs
for name, (bams, peaks) in samples.items():
    # Load peak file from bed files centered on motif, make window around
    try:
        peaksInt = pybedtools.BedTool(peaks)  # .slop(genome=genome, b=windowWidth / 2)
    except:
        print("Sample's peaks were not found: %s" % name)
        continue

    # Exclude peaks in gaps or repeats
    peaksInt.intersect(b=gapsRepeats, v=True, wa=True)
    peaksInt = bedToolsInterval2GenomicInterval(peaksInt)

    # Get coverage over all signals (self and histone for CM)
    for i in range(len(bams)):
        if i == 0:
            exportName = name
        else:
            exportName = name + "_histones"

        # get cuts
        if not os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
            print("Getting cuts for %s" % exportName)
            bam = HTSeq.BAM_Reader(bams[i])

            cov = coverage(bam, peaksInt, fragmentsize, strand_specific=False)

            # Make dataframe
            df = pd.DataFrame(np.vstack(cov.values()), index=cov.keys())

            # Save raw data
            savePandas(os.path.join(plotsDir, "pickles", exportName + ".pdy"), df)
        else:
            df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()

        # Centipede footprinting
        print("calling footprints for %s" % exportName)
        annot = pd.read_csv(re.sub("motifCentered", "motifAnnotated", peaks), sep="\t", header=None)

        annot[3] = 1
        annot.index = annot[0]
        annot = annot.ix[df.index]

        # filter
        annot[annot[2] > -1.000000e+10]
        df = df.ix[annot.index]

        # sort both in the same way
        df.sort(inplace=True)
        annot.sort(inplace=True)

        # get 200bp in middle of window of cuts
        n = len(df.columns)
        r = (len(df.columns) / 2 - 200, len(df.columns) / 2 + 199)

        # call
        f = callFootprints(df.loc[:, r[0]:r[1]], annot[[3, 2]])
        foots[exportName] = f
        pickle.dump(foots, open(os.path.join(plotsDir, "pickles", "footprintProbs.pickle"), 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

        # export confident foots as bed
        # get original bed
        p = pd.read_csv(peaks, sep="\t", header=None)
        p.columns = ["chrom", "start", "end", "name", "score", "strand"]

        df = df.loc[f > 0.9, :]
        df["name"] = df.index

        df = p.merge(df)[["chrom", "start", "end"]]

        df.to_csv(os.path.join(peaksDir, exportName + "_footprints.bed"), sep="\t", index=False)

# Plot
foots = pickle.load(open(os.path.join(plotsDir, "pickles", "footprintProbs.pickle"), 'rb'))

# Compare overlap of footprints predicted independently by each technique
# 1. export beds
# 2. overlap
# 3. venn diagrams
fig, axis = plt.subplots(2, 2, sharey=True, figsize=(10, 8))
i = 0
for TF in ["CTCF", "GATA1", "PU1", "REST"]:
    d = pd.DataFrame([
        np.ndarray.flatten(foots["CM_" + TF]),
        np.ndarray.flatten(foots["DNase_" + TF]),
        np.ndarray.flatten(foots["ATAC_" + TF]),
    ]).T
    d.columns = ["CM_" + TF, "DNase_" + TF, "ATAC_" + TF]

    # Plot
    if i is 0:
        j, k = (0, 0)
    elif i is 1:
        j, k = (0, 1)
    elif i is 2:
        j, k = (1, 0)
    elif i is 3:
        j, k = (1, 1)

    axis[j][k].boxplot([d["CM_" + TF], d["DNase_" + TF], d["ATAC_" + TF]])
    axis[j][k].set_title(TF)
    axis[j][k].set_ylabel("Footprint probability")
    i += 1
plt.savefig(os.path.join(plotsDir, "footprint.probs.independent.boxplot.pdf"), bbox_inches='tight')
plt.close()


fig, axis = plt.subplots(2, 2, sharey=True, figsize=(10, 8))
i = 0
for TF in ["CTCF", "GATA1", "PU1", "REST"]:
    d = pd.DataFrame([
        np.ndarray.flatten(foots["CM_" + TF]),
        np.ndarray.flatten(foots["DNase_" + TF]),
        np.ndarray.flatten(foots["ATAC_" + TF]),
    ]).T
    d.columns = ["CM_" + TF, "DNase_" + TF, "ATAC_" + TF]

    # Plot
    if i is 0:
        j, k = (0, 0)
    elif i is 1:
        j, k = (0, 1)
    elif i is 2:
        j, k = (1, 0)
    elif i is 3:
        j, k = (1, 1)

    sns.violinplot(d.dropna(), inner="points", inner_kws={"ms": 1}, ax=axis[j][k])
    sns.despine(left=True, bottom=True, ax=axis[j][k])
    i += 1
plt.savefig(os.path.join(plotsDir, "footprint.probs.independent.violinplot.pdf"), bbox_inches='tight')
plt.close()

# Compare number of footprints called by each technique in each technique
# 1. Save probs as dict for each technique
# 2. put together in dataframe (cols - techniques)
# 3. plot distributions of probabilities as factor
fig, axis = plt.subplots(2, 2, sharey=True, figsize=(10, 8))
i = 0
for TF in ["CM_CTCF", "CM_GATA1", "CM_PU1", "CM_REST"]:
    if i is 0:
        j, k = (0, 0)
    elif i is 1:
        j, k = (0, 1)
    elif i is 2:
        j, k = (1, 0)
    elif i is 3:
        j, k = (1, 1)

    d = pd.DataFrame([
        np.ndarray.flatten(foots[TF]),
        np.ndarray.flatten(foots[TF + "_histones"])
    ])
    d = d.T
    d.columns = [TF, TF + "_histones"]

    axis[j][k].boxplot([d[TF], d[TF + "_histones"]])
    axis[j][k].set_title(TF)
    axis[j][k].set_ylabel("Footprint probability")
    i += 1
plt.savefig(os.path.join(plotsDir, "footprint.probs.CM_vs_CMhistones.pdf"), bbox_inches='tight')
plt.close()


# Venn diagrams
"""
for TF in CTCF GATA1 PU1 REST
do
    for TECH1 in ATAC CM DNase
    do
        for TECH2 in ATAC CM DNase
        do
            T=`wc -l ${TECH1}_${TF}_footprints.bed`
            O=`bedtools intersect -u \
            -a ${TECH1}_${TF}_footprints.bed \
            -b ${TECH2}_${TF}_footprints.bed \
            | wc -l`
            echo $TECH1 $TECH2 $TF $O
        done
    done
done

ATAC    ATAC    CTCF    16290
ATAC    CM  CTCF    8481
ATAC    DNase   CTCF    10372
CM  ATAC    CTCF    9017
CM  CM  CTCF    16717
CM  DNase   CTCF    9052
DNase   ATAC    CTCF    10462
DNase   CM  CTCF    8596
DNase   DNase   CTCF    15419
ATAC    ATAC    GATA1   13187
ATAC    CM  GATA1   2817
ATAC    DNase   GATA1   6148
CM  ATAC    GATA1   2847
CM  CM  GATA1   13727
CM  DNase   GATA1   3072
DNase   ATAC    GATA1   6317
DNase   CM  GATA1   3119
DNase   DNase   GATA1   16411
ATAC    ATAC    PU1 3099
ATAC    CM  PU1 945
ATAC    DNase   PU1 781
CM  ATAC    PU1 922
CM  CM  PU1 6429
CM  DNase   PU1 744
DNase   ATAC    PU1 823
DNase   CM  PU1 792
DNase   DNase   PU1 3318
ATAC    ATAC    REST    429
ATAC    CM  REST    260
ATAC    DNase   REST    162
CM  ATAC    REST    291
CM  CM  REST    1201
CM  DNase   REST    243
DNase   ATAC    REST    160
DNase   CM  REST    232
DNase   DNase   REST    430
"""
o = pd.read_clipboard(header=None)
o.columns = ['tec1', 'tec2', 'ip', 'count']

for g in o.groupby(['ip']).groups.items():
    import matplotlib_venn as venn
    plt.figure()
    venn.venn3(subsets=o.ix[g[1]]['count'].tolist(), set_labels=(o.ix[g[1]]['tec1'].unique()))
    plt.savefig(os.path.join(plotsDir, "footprint.overlap.%s.pdf" % g[0]), bbox_inches='tight')


# Pie charts overlap with ChIP
"""
for TF in CTCF GATA1 PU1 REST
do
    for TECH1 in ATAC CM DNase
    do
        T=`wc -l ${TECH1}_${TF}_footprints.bed`
        O=`bedtools intersect -u \
        -a ${TECH1}_${TF}_footprints.bed \
        -b K562_10M_CHIP_${TF}_nan_nan_0_0_hg19/K562_10M_CHIP_${TF}_nan_nan_0_0_hg19_peaks.narrowPeak \
        | wc -l`
        echo $TECH1 $TF $T $O
    done
done


ATAC    CTCF    16290   14900
CM  CTCF    16717   16420
DNase   CTCF    15419   14614
ATAC    GATA1   13187   7357
CM  GATA1   13727   11899
DNase   GATA1   16411   8409
ATAC    PU1 3099    2308
CM  PU1 6429    6383
DNase   PU1 3318    2280
ATAC    REST    429 345
CM  REST    1201    1024
DNase   REST    430 318


"""
o = pd.read_clipboard(header=None)
o.columns = ['tec', 'ip', 'total', 'count']

o['over'] = o['count'] / o['total']
o['leftover'] = 1 - (o['count'] / o['total'])


for ip in o['ip'].unique():
    i = 0
    fig, axis = plt.subplots(1, 3)
    for tec in o['tec'].unique():
        s = o[
            (o['ip'] == ip) &
            (o['tec'] == tec)
        ]
        axis[i].pie(s['over'], s['leftover'], colors=['white', 'black'])
        i += 1
    plt.savefig(os.path.join(plotsDir, "footprint.ChIP-overlap.%s.pie.pdf" % ip), bbox_inches='tight')
    plt.close()

g = sns.FacetGrid(o, col="ip", size=4, aspect=.5)
g.map(sns.barplot, "tec", "over")
plt.savefig(os.path.join(plotsDir, "footprint.ChIP-overlap.bar.pdf"), bbox_inches='tight')
