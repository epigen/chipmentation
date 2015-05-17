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


def bedToolsInterval2GenomicInterval(bedtool):
    """
    Given a pybedtools.BedTool object, returns dictionary of HTSeq.GenomicInterval objects.
    """
    intervals = list()
    for iv in bedtool:
        if iv.name == "":
            iv.name = "".join([str(iv.chrom), str(iv.start), str(iv.end)])

        if iv.strand == "+" or iv.strand == 0 or iv.strand == str(0):
            intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "+")
        elif iv.strand == "-" or iv.strand == 0 or iv.strand == str(1):
            intervals[iv.name] = HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end, "-")
        else:
            intervals.append(HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end))
    return intervals


def coverage(bam, intervals, fragmentsize, orientation=True, duplicates=True):
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
    n = len(intervals)
    # Initialize empty array for this feature
    profile = np.zeros(intervals[0].length, dtype=np.float64)

    try:
        i = 0
        for feature in intervals:
            if i % 1000 == 0:
                print(n - i)
           
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
                profile[start_in_window: end_in_window] += 1

            i += 1
    except KeyboardInterrupt:
        return profile
    return profile


def pairwiseDistances(intervals):
    """
    Brute force
    """
    import itertools
    from collections import Counter
    distances = Counter()
  
    try:
        for i, (d1, d2) in enumerate(itertools.combinations(intervals, 2)):
            if i % 1000000 == 0:
                print(i)
            if d1.chrom != d2.chrom:
                continue
            # distance end-to-end
            if abs(d1.end <= d2.start):
                if abs(d2.start - d1.end) < 500:
                    distances[abs(d2.start - d1.end)] += 1
            else:
                if abs(d1.start - d2.end) < 500:
                    distances[abs(d1.start - d2.end)] += 1

        return distances
    except KeyboardInterrupt:
        return distances


def pairwiseDistances(intervals):
    """
    """
    import itertools
    from collections import Counter
    distances = Counter()
  
    try:
        for i, (d1, d2) in enumerate(itertools.combinations(intervals, 2)):
            if i % 1000000 == 0:
                print(i)
            if d1.chrom != d2.chrom:
                continue
            # distance end-to-end
            if abs(d1.end <= d2.start):
                if abs(d2.start - d1.end) < 500:
                    distances[abs(d2.start - d1.end)] += 1
            else:
                if abs(d1.start - d2.end) < 500:
                    distances[abs(d1.start - d2.end)] += 1

        return distances
    except KeyboardInterrupt:
        return distances


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
        return "#FF4444"
    elif "ATAC" in name:
        return "#001f07"
    else:
        return "#e69f00"


# Define variables
# projectRoot = "/projects/chipmentation/"
projectRoot = "/media/afr/cemm-backup1/chipmentation/"
# projectRoot = "/home/arendeiro/chipmentation/"
bamsDir = os.path.join(projectRoot, "data", "mapped/")
resultsDir = os.path.join(projectRoot, "results")
plotsDir = os.path.join(resultsDir, "plots")

CM = os.path.join(bamsDir, "K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.trimmed.bowtie2.shifted.dups.bam")
ATAC = os.path.join(bamsDir, "K562_50K_ATAC_nan_nan_nan_0_0_hg19.trimmed.bowtie2.shifted.dups.bam")
DNase = os.path.join(bamsDir, "wgEncodeUwDnaseK562Aln.merged.bam")
MNase = os.path.join(bamsDir, "wgEncodeSydhNsomeK562Aln.merged.bam")

samples = {"CM": CM, "ATAC": ATAC, "DNASE": DNase, "MNASE": MNase}

dyadsFile = "/media/afr/cemm-backup1/chipmentation/data/nucleoATAC/10M_CM_H3K4ME1_PE/dyads_3col.bed"

windowRange = (-400, 400)
fragmentsize = 1
duplicates = False
n_clusters = 5
genome = "hg19"

windowWidth = abs(windowRange[0]) + abs(windowRange[1])


dyads = pybedtools.BedTool(dyadsFile)
dyads = dyads.slop(genome=genome, b=windowWidth / 2)
# dyads = bedToolsInterval2GenomicInterval(dyads)

intervals = list()
for iv in dyads:
    intervals.append(HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end))
dyads = intervals

aveSignals = pd.DataFrame(columns=['sample'])

# Loop through all samples, compute coverage in peak regions centered on motifs
for name, bam in samples.items():
    print(name)

    if name not in aveSignals['sample'].unique():
        # Load bam
        bamFile = HTSeq.BAM_Reader(bam)

        cov = coverage(bamFile, dyads, fragmentsize)

        # Average signals
        df = pd.DataFrame({
            "sample": name,
            "average": cov,                            # both strands
            "x": range(windowRange[0], windowRange[1] + 1)
        })
        aveSignals = pd.concat([aveSignals, df])
        pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "dyad_coverage-H3K4me1-all.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)


aveSignals = pickle.load(open(os.path.join(plotsDir, "pickles", "dyad_coverage-H3K4me1-all.pickle"), "r"))

# Plots
aveSignals.groupby(['sample']).plot(['x'], labels=['sample'])

fig, axis = plt.subplots(2, 1, sharex=True, figsize=(10, 8))
for i in range(len(aveSignals['sample'].unique())):
    sample = aveSignals['sample'].unique()[i]
    sub = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['x'] > - 155) &
        (aveSignals['x'] < 145)
    ]

    if i is 0:
        j, k = (0, 0)
    elif i is 1:
        j, k = (0, 1)
    elif i is 2:
        j, k = (1, 0)
    elif i is 3:
        j, k = (1, 1)

    for s in [sample, "DNase", "ATAC"]:
        if s == sample:
            sub2 = sub[sub['sample'] == s]
            colour = colourPerFactor(s)
            axis[i].plot(sub2['x'] + 5, sub2.average, color=colour, linestyle="-", alpha=1)
        else:
            sub2 = sub[sub['sample'].str.contains(s)]
            colour = colourPerFactor(s)
            axis[i].plot(sub2['x'] + 5, sub2.average, color=colour, linestyle="-", alpha=0.35)

    # subplot attributes
    axis[i].set_title(sample)
    axis[i].set_xlabel("Distance to dyad (bp)")
    axis[i].set_ylabel("Insertion frequency")
plt.savefig(os.path.join(plotsDir, "dyad_coverage.signals-alldata.pdf"), bbox_inches='tight')
plt.close()


# Distances between nucleosomes
nucleosomesFile = "/media/afr/cemm-backup1/chipmentation/data/nucleoATAC/10M_CM_H3K4ME1_PE/nucleosomes.bed"
nucleosomes = pybedtools.BedTool(nucleosomesFile)

intervals = list()
for iv in nucleosomes:
    intervals.append(HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end))
nucleosomes = intervals

dists = pairwiseDistances(nucleosomes)

pickle.dump(dists, open(os.path.join(plotsDir, "nucleosome-distances.pickle"), 'wb'))


df = pd.read_csv(nucleosomesFile, sep="\t", header=None)
df = df[[0, 1, 2]]
df.columns = ['chrom', 'start', 'end']

c = list()

for chrom, indices in df.groupby('chrom').groups.items():
    df2 = df.ix[indices].reset_index(drop=True)

    df2.apply(lambda x: [x - x[i] for i in range(len(df2))], axis=0)

    for i in range(len(df2)):
        s = df2['end'][i] - df2['start']
        c += s.tolist()
        

# Positive NucleoATAC signal around TSSs
import BedGraphReader 
import numpy as np
import pybedtools

bg_file = "/media/afr/cemm-backup1/chipmentation/data/nucleoATAC/10M_CM_H3K4ME1_PE/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.nucleoatac_signal.smooth.bedgraph"
bg_file = "/media/afr/cemm-backup1/chipmentation/data/nucleoATAC/10M_CM_H3K4ME1_PE/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.occ.bedgraph"

chromsizes = pybedtools.get_chromsizes_from_ucsc('hg19')

chromsizes = {chrom: end for chrom, (start, end) in chromsizes.items()}
genome = BedGraphReader.load_bedgraph(bg_file, chromsizes)  # actually read bed graph file

ctcf = pybedtools.BedTool("/media/afr/cemm-backup1/chipmentation/data/peaks/K562_10M_CHIP_CTCF_nan_nan_0_0_hg19/K562_10M_CHIP_CTCF_nan_nan_0_0_hg19_peaks.motifCentered.bed")
cage = pybedtools.BedTool("/media/afr/cemm-backup1/chipmentation/annotation/hg19.cage_peak.robust.TATA_Annotated.expr_Annotated.K562_Expressed.1kb.tsv")

profile = np.zeros(len(cage[0]))
for tss in cage:
    if tss.chrom in genome.keys():
        signal = genome[tss.chrom][tss.start: tss.end]
        if len(signal[signal < 0]) > 1:
            print(1)
        signal[signal < 0] = 0
        profile += signal
        if tss.strand == "+":
            profile += signal
        else:
            profile += signal[::-1]


# ChIPmentation coverage in MNase dyads
bam = "/media/afr/cemm-backup1/chipmentation/data/mapped/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.trimmed.bowtie2.shifted.dups.bam"
bamFile = HTSeq.BAM_Reader(bam)

dyads = "/media/afr/cemm-backup1/dyads/dyads.bed"
dyads = pybedtools.BedTool(dyads)
dyads = dyads.slop(b=300, genome="hg19")
intervals = list()
for iv in dyads:
    intervals.append(HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end))
dyads = intervals

cov = coverage(bamFile, dyads, fragmentsize)

pickle.dump(cov, open(os.path.join(plotsDir, "pickles", "dyad_coverage-H3K4me1_in_MNase_dyads-all.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
cov = pickle.load(open(os.path.join(plotsDir, "pickles", "dyad_coverage-H3K4me1_in_MNase_dyads-all.pickle"), "r"))


plt.plot(np.array(range(-150, 150)), cov[145:445], linestyle="-", alpha=1)
plt.title("CM")
plt.xlabel("Distance to MNase-seq dyad (bp)")
plt.ylabel("Insertion frequency")
plt.savefig(os.path.join(plotsDir, "dyad_coverage.H3K4me1_in_MNase_dyads-alldata.pdf"), bbox_inches='tight')
plt.close()



# MNase in ChIPmentation dyads
bam = "/media/afr/cemm-backup1/chipmentation/data/mapped/wgEncodeSydhNsomeK562Aln.merged.bam"
bamFile = HTSeq.BAM_Reader(bam)

dyadsFile = "/media/afr/cemm-backup1/chipmentation/data/nucleoATAC/10M_CM_H3K4ME1_PE/dyads_3col.bed"
dyads = pybedtools.BedTool(dyadsFile)
dyads = dyads.slop(genome=genome, b=windowWidth / 2)
intervals = list()
for iv in dyads:
    intervals.append(HTSeq.GenomicInterval(iv.chrom, iv.start, iv.end))
dyads = intervals

cov = coverage(bamFile, dyads, fragmentsize)
pickle.dump(cov, open(os.path.join(plotsDir, "pickles", "dyad_coverage-MNase_in_H3K4me1_dyads-all.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
mn = pickle.load(open(os.path.join(plotsDir, "pickles", "dyad_coverage-MNase_in_H3K4me1_dyads-all.pickle"), "r"))

plt.plot(range(-200, 200), cov[100:500])
plt.savefig(os.path.join(plotsDir, "dyad_coverage.MNase_in_H3K4me1_dyads-alldata.pdf"), bbox_inches='tight')


cov = coverage(bamFile, dyads, 36)
pickle.dump(cov, open(os.path.join(plotsDir, "pickles", "dyad_coverage-MNase_in_H3K4me1_dyads-36bp-all.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
mn = pickle.load(open(os.path.join(plotsDir, "pickles", "dyad_coverage-MNase_in_H3K4me1_dyads-36bp-all.pickle"), "r"))

