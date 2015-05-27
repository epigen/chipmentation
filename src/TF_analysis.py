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
# import MOODS
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


def coverage(bam, intervals, fragmentsize, orientation=True, duplicates=False, strand_specific=False, switchChromsNames=False):
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
    try:
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
    except KeyboardInterrupt:
        return cov

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
        return "#2d3730ff"
    else:
        raise ValueError


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


def getPWMscores(peaks, PWM, fasta, peaksFile):
    # Get nucleotides
    seq = peaks.sequence(s=True, fi=fasta)
    with open(seq.seqfn) as handle:
        seqs = handle.read().split("\n")[1::2]  # get fasta sequences

    # Match strings with PWM
    scores = list()
    for sequence in seqs:
        result = MOODS.search(sequence, [PWM], 30)
        scores.append(np.array([j for i, j in result[0]]))

    names = open(peaksFile).read().split("\n")
    names.pop(-1)
    names = [line.split("\t")[3] for line in names]

    return pd.DataFrame(scores, index=names, columns=range(-990, 990))


def zscore(x):
    return (x - np.mean(x)) / np.std(x)


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


def plotFootprintModel(cuts, annot):
    """
    Plot footprint model.
    Requires dataframe with cuts and dataframe with annotation (2> cols).
    """
    import rpy2.robjects as robj  # for ggplot in R
    import rpy2.robjects.pandas2ri  # for R dataframe conversion

    # Plot with R
    footprint = robj.r("""
    library(CENTIPEDE)

    function(cuts, annot) {
        imageCutSites(cuts[order(centFit$PostPr),][c(1:100, (dim(cuts)[1]-100):(dim(cuts)[1])),])
        plotProfile(centFit$LambdaParList[[1]],Mlen=2)
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
projectRoot = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/"
# projectRoot = "/media/afr/cemm-backup/chipmentation/"
# projectRoot = "/home/arendeiro/chipmentation/"
bamsDir = os.path.join(projectRoot, "data", "mapped/")
peaksDir = os.path.join(projectRoot, "data", "peaks/")
resultsDir = os.path.join(projectRoot, "results")
plotsDir = os.path.join(resultsDir, "footprints")
DNase = os.path.join(bamsDir, "wgEncodeUwDnaseK562Aln.merged.bam")
MNase = os.path.join(bamsDir, "wgEncodeSydhNsomeK562Aln.merged.bam")
NexteraBackground = "/home/arendeiro/PGA1Nextera/PGA_0001_Nextera-2.bam"

# Get samples
samples = pd.read_csv(os.path.abspath(projectRoot + "chipmentation.replicates.annotation_sheet.csv"))

# Replace input sample
# samples.loc[:, "controlSampleFilePath"] = samples["controlSampleFilePath"].apply(lambda x: re.sub("K562_10M_CM_IGG_nan_nan_0_0_hg19", "K562_10M_CM_IGG_nan_nan_1_0_hg19", str(x)))

# Replace peaks of low cell number TFs
samples.loc[samples['sampleName'] == "K562_500K_CM_CTCF_nan_nan_0_0_hg19", "peakFile"] = peaksDir + "K562_10M_CM_CTCF_nan_nan_0_0_hg19/K562_10M_CM_CTCF_nan_nan_0_0_hg19_peaks.narrowPeak"
samples.loc[samples['sampleName'] == "K562_500K_CM_PU1_nan_nan_0_0_hg19", "peakFile"] = peaksDir + "K562_10M_CM_PU1_nan_nan_0_0_hg19/K562_10M_CM_PU1_nan_nan_0_0_hg19_peaks.narrowPeak"
samples.loc[samples['sampleName'] == "K562_100K_CM_GATA1_nan_nan_0_0_hg19", "peakFile"] = peaksDir + "K562_10M_CM_GATA1_nan_nan_0_0_hg19/K562_10M_CM_GATA1_nan_nan_0_0_hg19_peaks.narrowPeak"
samples.loc[samples['sampleName'] == "K562_100K_CM_REST_nan_nan_0_0_hg19", "peakFile"] = peaksDir + "K562_10M_CM_REST_nan_nan_0_0_hg19/K562_10M_CM_REST_nan_nan_0_0_hg19_peaks.narrowPeak"

# # replace bam path with local
# fhPath = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/"
# samples.loc[:, "filePath"] = samples["filePath"].apply(lambda x: re.sub(fhPath, bamsDir, str(x)))
# samples.loc[:, "controlSampleFilePath"] = samples["controlSampleFilePath"].apply(lambda x: re.sub(fhPath, bamsDir, str(x)))

# # replace peaks path with local
# fhPath = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/"
# samples.loc[:, "peakFile"] = samples["peakFile"].apply(lambda x: re.sub(fhPath, peaksDir, str(x)))

# subset samples
sampleSubset = samples[
    (samples["technique"] == "CM") &
    (samples["ip"].str.contains("CTCF|GATA1|PU1|REST")) &
    (samples["numberCells"] == "10M") &
    (samples["biologicalReplicate"] == 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"])

# subset samples into signals
signals = samples[
    samples["sampleName"].str.contains(
        "K562_50K_ATAC_nan_nan_nan_0_0_hg19"
    )
].reset_index(drop=True)
signals = signals.sort(["ip", "technique"])

signals = signals.append(pd.Series(data=["DNase", DNase], index=["sampleName", "filePath"]), ignore_index=True)
# signals = signals.append(pd.Series(data=["MNase", MNase], index=["sampleName", "filePath"]), ignore_index=True)

windowRange = (-1000, 1000)
fragmentsize = 1
duplicates = False
n_clusters = 5

windowWidth = abs(windowRange[0]) + abs(windowRange[1])

# gapsRepeats = pybedtools.BedTool(os.path.join("/media/afr/cemm-backup", "reference/hg19/hg19_gapsRepeats.bed"))
gapsRepeats = pybedtools.BedTool("/home/arendeiro/reference/Homo_sapiens/hg19_gapsRepeats.bed")

# Get genome fasta file
# fasta = pybedtools.BedTool("/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa")
fasta = pybedtools.BedTool("/home/arendeiro/reference/Homo_sapiens/forBowtie2/hg19.fa")

# Get PWM
PWM = [
    [0.83593967562, 0.9612721605805, 0.924068101918, 1.06685938079, 0.797339497119, 1.18221399883, 0.8166032303705, 0.731357810088, 0.8009085680465, 0.6413220540835, 1.24263470188, 1.27537587614, 1.093537759345, 1.270194775255, 0.607808747652, 1.1045452182715, 0.830270118282, 0.9203172509935, 0.7788696601795, 0.820042458443, 1.071536938715],
    [1.058289907845, 0.9590219891645, 1.307638114905, 0.9474308687955, 1.625234242265, 0.703349925951, 1.126706911395, 1.15491741682, 1.29497649068, 1.446820024825, 0.669382028924, 0.673327304706, 0.85726490824, 0.842842178544, 1.703477693555, 0.883047249808, 0.911871662974, 1.071061400505, 1.12081896184, 1.356389737925, 1.075154356595],
    [1.075154356595, 1.356389737925, 1.12081896184, 1.071061400505, 0.911871662974, 0.883047249808, 1.703477693555, 0.842842178544, 0.85726490824, 0.673327304706, 0.669787971405, 1.446820024825, 1.29497649068, 1.15491741682, 1.126706911395, 0.703349925951, 1.625234242265, 0.9474308687955, 1.307638114905, 0.9590219891645, 1.058289907845],
    [1.071536938715, 0.820042458443, 0.7788696601795, 0.9203172509935, 0.830270118282, 1.1045452182715, 0.607808747652, 1.270194775255, 1.093537759345, 1.27537587614, 1.21553439588, 0.6413220540835, 0.8009085680465, 0.731357810088, 0.8166032303705, 1.18221399883, 0.797339497119, 1.06685938079, 0.924068101918, 0.9612721605805, 0.83593967562]
]


# Loop through all samples, compute coverage in peak regions centered on motifs
aveSignals = pd.DataFrame(columns=["sample", "signal"])
foots = dict()

for i in range(len(sampleSubset)):
    sampleName = sampleSubset['sampleName'][i]
    print(sampleName)

    if "100K" in sampleName:
        peaks = re.sub("100K", "10M", sampleSubset['peakFile'][i]) + "motifCentered.bed"
    elif "500K" in sampleName:
        peaks = re.sub("500K", "10M", sampleSubset['peakFile'][i]) + "motifCentered.bed"
    else:
        peaks = sampleSubset['peakFile'][i]

    motifCentered = re.sub("\..*", ".", peaks) + "motifCentered.bed"

    # replace with CHIP (do it on CHIP peaks)
    motifCentered = re.sub("CM", "CHIP", motifCentered)
    print(motifCentered)

    # Load peak file from bed files centered on motif, make window around
    try:
        peaks = pybedtools.BedTool(motifCentered)  # .slop(genome=genome, b=windowWidth / 2)
    except:
        print("Sample's peaks were not found: %s" % sampleName)
        continue

    peaksInt = bedToolsInterval2GenomicInterval(peaks)

    # Filter peaks near chrm borders
    for peak, interval in peaksInt.iteritems():
        if interval.length < windowWidth:
            peaksInt.pop(peak)

    # Get PWM score
    # signalName = "PWM"
    # print(sampleName, signalName)
    # exportName = "-".join([sampleName, signalName])

    # if signalName not in aveSignals[aveSignals['sample'] == sampleName]['signal'].unique():
    #     if not os.path.isfile(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy")):
    #         ctrl = getPWMscores(peaks, PWM, fasta, motifCentered)
    #         # Save raw data
    #         savePandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy"), ctrl)
    # else:
    #     ctrl = loadPandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy")).copy()

    # df = pd.DataFrame({
    #     "sample": sampleName,
    #     "signal": signalName,
    #     "average": ctrl.apply(np.mean, axis=0),
    #     "x": ctrl.columns
    # })
    # aveSignals = pd.concat([aveSignals, df])
    # pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.subset.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # Get Nextera background
    signalName = "NexteraBackground"
    print(sampleName, signalName)
    exportName = "-".join([sampleName, signalName])

    if signalName not in aveSignals[aveSignals['sample'] == sampleName]['signal'].unique():
        if not os.path.isfile(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy")):

            bamFile = HTSeq.BAM_Reader(NexteraBackground)

            ctrl = coverage(bamFile, peaksInt, fragmentsize, strand_specific=True, switchChromsNames=True)

            # Make multiindex dataframe
            levels = [ctrl.keys(), ["+", "-"]]
            labels = [[y for x in range(len(ctrl)) for y in [x, x]], [y for x in range(len(ctrl.keys())) for y in (0, 1)]]
            index = pd.MultiIndex(labels=labels, levels=levels, names=["peak", "strand"])
            ctrl = pd.DataFrame(np.vstack(ctrl.values()), index=index)
            ctrl.columns = range(windowRange[0], windowRange[1])

            # Save raw data
            savePandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy"), ctrl)
    else:
        ctrl = loadPandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy")).copy()

    df = pd.DataFrame({
        "sample": sampleName,
        "signal": signalName,
        "average": ctrl.apply(np.mean, axis=0),
        "x": ctrl.columns
    })
    aveSignals = pd.concat([aveSignals, df])
    pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.subset.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # Get self coverage
    signalName = sampleName
    exportName = "-".join([sampleName, signalName])
    print(sampleName, signalName)
    if signalName not in aveSignals[aveSignals['sample'] == sampleName]['signal'].unique():
        if not os.path.isfile(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy")):
            # Load bam
            bamFile = HTSeq.BAM_Reader(sampleSubset['filePath'][i])

            cov = coverage(bamFile, peaksInt, fragmentsize, strand_specific=True)

            # Make multiindex dataframe
            levels = [cov.keys(), ["+", "-"]]
            labels = [[y for x in range(len(cov)) for y in [x, x]], [y for x in range(len(cov.keys())) for y in (0, 1)]]
            index = pd.MultiIndex(labels=labels, levels=levels, names=["peak", "strand"])
            df = pd.DataFrame(np.vstack(cov.values()), index=index)
            df.columns = range(windowRange[0], windowRange[1])

            # Save raw data
            savePandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy"), df)
        else:
            df = loadPandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy")).copy()

        # Average signal
        ave = pd.DataFrame({
            "sample": sampleName,
            "signal": signalName,
            "average": df.apply(np.mean, axis=0),                            # both strands
            "positive": df.ix[range(0, len(df), 2)].apply(np.mean, axis=0),  # positive strand
            "negative": df.ix[range(1, len(df), 2)].apply(np.mean, axis=0),  # negative strand
            "x": df.columns
        })
        aveSignals = pd.concat([aveSignals, ave])
        pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.subset.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        # Normalize signal
        idx = ctrl.index.levels[0]
        ctrlT = ctrl.ix[range(0, len(ctrl), 2)].reset_index(drop=True) + ctrl.ix[range(1, len(ctrl), 2)].reset_index(drop=True)
        ctrlT.index = idx

        idx = df.index.levels[0]
        treat = df.ix[range(0, len(df), 2)].reset_index(drop=True) + df.ix[range(1, len(df), 2)].reset_index(drop=True)
        treat.index = idx

        norm = treat / np.exp(ctrlT.apply(zscore, axis=1))

        # Save normalized data
        savePandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks-normalized.pdy"), norm)

        # Average normalized data
        ave = pd.DataFrame({
            "sample": sampleName,
            "signal": signalName + "-normalized",
            "average": norm.apply(sum, axis =0),                            # both strands
            "positive": 0,  # positive strand
            "negative": 0,  # negative strand
            "x": norm.columns
        })
        aveSignals = pd.concat([aveSignals, ave])
        pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.subset.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        # Centipede footprinting
        print("calling footprints for %s" % exportName)
        annot = pd.read_csv(re.sub("motifCentered", "motifAnnotated", motifCentered), sep="\t", header=None)

        annot[3] = 1
        annot.index = annot[0]

        # filter weird stuff out
        annot = annot[annot[2] > -1.000000e+10]
        norm = norm.ix[annot.index]

        # sort both in the same way
        norm.sort(inplace=True)
        annot.sort(inplace=True)

        # get 200bp in middle of window of cuts
        n = len(norm.columns)
        r = (- 200, 199)
        norm = norm.loc[:, r[0]:r[1]]
        norm = norm.dropna()

        # call
        annot = annot.ix[norm.index]
        footprintProbs = callFootprints(norm, annot[[3, 2]])
        foots[exportName] = footprintProbs
        pickle.dump(foots, open(os.path.join(plotsDir, "pickles", "footprintProbs-normalized.pickle"), 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    # Get coverage over other signals
    for j in range(len(signals)):
        signalName = signals['sampleName'][j]
        print(sampleName, signalName)

        exportName = "-".join([sampleName, signalName])
        if signalName not in aveSignals[aveSignals['sample'] == sampleName]['signal'].unique():
            if not os.path.isfile(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy")):
                print(sampleName, signalName)
                bamFile = HTSeq.BAM_Reader(signals['filePath'][j])

                cov = coverage(bamFile, peaksInt, fragmentsize, strand_specific=True)

                # Make multiindex dataframe
                levels = [cov.keys(), ["+", "-"]]
                labels = [[y for x in range(len(cov)) for y in [x, x]], [y for x in range(len(cov.keys())) for y in (0, 1)]]
                index = pd.MultiIndex(labels=labels, levels=levels, names=["peak", "strand"])
                df = pd.DataFrame(np.vstack(cov.values()), index=index)
                df.columns = range(windowRange[0], windowRange[1])

                # Save raw data
                savePandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy"), df)
        else:
            df = loadPandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks.pdy")).copy()

        ave = pd.DataFrame({
            "sample": sampleName,
            "signal": signalName,
            "average": df.apply(np.mean, axis=0),                            # both strands
            "positive": df.ix[range(0, len(df), 2)].apply(np.mean, axis=0),  # positive strand
            "negative": df.ix[range(1, len(df), 2)].apply(np.mean, axis=0),  # negative strand
            "x": df.columns
        })
        aveSignals = pd.concat([aveSignals, ave])
        pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.subset.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        # Normalize signal (not for DNASE!)
        idx = ctrl.index.levels[0]
        ctrlT = ctrl.ix[range(0, len(ctrl), 2)].reset_index(drop=True) + ctrl.ix[range(1, len(ctrl), 2)].reset_index(drop=True)
        ctrlT.index = idx

        idx = df.index.levels[0]
        treat = df.ix[range(0, len(df), 2)].reset_index(drop=True) + df.ix[range(1, len(df), 2)].reset_index(drop=True)
        treat.index = idx

        if signalName != "DNase":
            norm = treat / np.exp(ctrlT.apply(zscore, axis=1))
        else:
            norm = treat

        # Save normalized data
        savePandas(os.path.join(plotsDir, "pickles", exportName + "-ChIPpeaks-normalized.pdy"), norm)

        # Average normalized data
        ave = pd.DataFrame({
            "sample": sampleName,
            "signal": signalName + "-normalized",
            "average": norm.apply(sum, axis =0),                            # both strands
            "positive": 0,  # positive strand
            "negative": 0,  # negative strand
            "x": norm.columns
        })
        aveSignals = pd.concat([aveSignals, ave])
        pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.subset.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

        # Centipede footprinting
        if signalName in ['K562_50K_ATAC_nan_nan_nan_0_0_hg19', 'DNase']:
            print("calling footprints for %s" % exportName)
            annot = pd.read_csv(re.sub("motifCentered", "motifAnnotated", motifCentered), sep="\t", header=None)

            annot[3] = 1
            annot.index = annot[0]

            # filter weird stuff out
            annot = annot[annot[2] > -1.000000e+10]
            norm = norm.ix[annot.index]

            # sort both in the same way
            norm.sort(inplace=True)
            annot.sort(inplace=True)

            # get 200bp in middle of window of cuts
            n = len(norm.columns)
            r = (- 200, 199)
            norm = norm.loc[:, r[0]:r[1]]
            norm = norm.dropna()

            # call
            annot = annot.ix[norm.index]
            footprintProbs = callFootprints(norm, annot[[3, 2]])
            foots[exportName] = footprintProbs
            pickle.dump(foots, open(os.path.join(plotsDir, "pickles", "footprintProbs-normalized.pickle"), 'wb'), protocol=pickle.HIGHEST_PROTOCOL)


# Done with collecting data
aveSignals = pickle.load(open(os.path.join(plotsDir, "pickles", "aveSignals.subset.pickle"), "r"))

# Plot footprint probabilities
foots = pickle.load(open(os.path.join(plotsDir, "pickles", "footprintProbs-normalized.pickle"), 'rb'))

fig, axis = plt.subplots(2, 2, sharey=True, figsize=(10, 8))
for i, TF in enumerate(["CTCF", "GATA1", "PU1", "REST"]):
    TF = "K562_10M_CM_" + TF + "_nan_nan_0_0_hg19"
    d = pd.DataFrame([
        np.ndarray.flatten(foots[TF + "-" + TF]),
        np.ndarray.flatten(foots[TF + "-DNase"]),
        np.ndarray.flatten(foots[TF + "-K562_50K_ATAC_nan_nan_nan_0_0_hg19"]),
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

    sns.violinplot(d.dropna(), ax=axis[j][k])
    sns.despine(left=True, bottom=True, ax=axis[j][k])
plt.savefig(os.path.join(plotsDir, "footprint.probs.independent.violinplot-nodots-normalized.pdf"), bbox_inches='tight')
plt.close()


# Look at aggregate signal
aveSignals = pickle.load(open(os.path.join(plotsDir, "pickles", "aveSignals.subset.pickle"), "r"))
aveSignals.drop_duplicates(inplace=True)
aveSignals["type"] = "raw"
aveSignals.reset_index(drop=True, inplace=True)

# Normalize
df = aveSignals.copy()
df2 = pd.DataFrame()

for sample in df['sample'].unique():
    for signal in df[df['sample'] == sample]['signal'].unique():
        if signal == "NexteraBackground":
            continue
        for strand in ["average"]:
            treat = df[
                (df["sample"] == sample) &
                (df["signal"] == signal) &
                (df['x'] >= -400) &
                (df['x'] <= 400)
            ][strand].tolist()

            ctrl = df[
                (df["sample"] == sample) &
                (df["signal"] == "NexteraBackground") &
                (df['x'] >= -400) &
                (df['x'] <= 400)
            ][strand].tolist()

            # add raw scalled
            tmp = pd.DataFrame()
            tmp[strand] = treat
            tmp['sample'] = sample
            tmp['signal'] = signal
            tmp['x'] = np.array(range(len(tmp))) - (len(tmp) / 2)
            tmp['type'] = "raw"
            df2 = df2.append(tmp, ignore_index=True)

            # smooth
            treatSmoothed = smooth(np.array(treat))

            # add raw scalled and smoothed
            tmp = pd.DataFrame()
            tmp[strand] = treatSmoothed
            tmp['sample'] = sample
            tmp['signal'] = signal
            tmp['x'] = np.array(range(len(tmp))) - (len(tmp) / 2)
            tmp['type'] = "rawSmooth"
            df2 = df2.append(tmp, ignore_index=True)

            # normalize (divide by log(NexteraBackground))
            norm = treat / np.exp(np.array(ctrl))

            # smooth
            normSmoothed = smooth(norm)

            tmp = pd.DataFrame()
            tmp[strand] = norm
            tmp['sample'] = sample
            tmp['signal'] = signal
            tmp['x'] = np.array(range(len(tmp))) - (len(tmp) / 2)
            tmp['type'] = "norm"
            df2 = df2.append(tmp, ignore_index=True)

            tmp = pd.DataFrame()
            tmp[strand] = normSmoothed
            tmp['sample'] = sample
            tmp['signal'] = signal
            tmp['x'] = np.array(range(len(tmp))) - (len(tmp) / 2)
            tmp['type'] = "normSmooth"
            df2 = df2.append(tmp, ignore_index=True)


# append normalized df to df
aveSignals = df2
aveSignals.reset_index(drop=True, inplace=True)


# plot unormalized signals
for sample in aveSignals['sample'].unique():
    fig = plt.figure()
    t = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['signal'] == sample) &
        (aveSignals['x'] >= -200) &
        (aveSignals['x'] <= 200)
    ]['average']

    atac = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['signal'] == "K562_50K_ATAC_nan_nan_nan_0_0_hg19") &
        (aveSignals['x'] >= -200) &
        (aveSignals['x'] <= 200)
    ]['average']

    dnase = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['signal'] == "DNase") &
        (aveSignals['x'] >= -200) &
        (aveSignals['x'] <= 200)
    ]['average']

    plt.plot(range(-200, 200), smooth(zscore(t), 15)[15:], label="ChIPmentation")
    plt.plot(range(-200, 200), smooth(zscore(atac), 15)[15:], label="ATAC-seq")
    plt.plot(range(-200, 200), smooth(zscore(dnase), 15)[15:], label="DNase-seq")
    plt.legend()
    plt.savefig(os.path.join(plotsDir, sample + ".all.footprints.pdf"), bbox_inches='tight')
    plt.close()


# plot normalized signals
for sample in aveSignals['sample'].unique():
    fig = plt.figure()
    t = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['signal'] == sample) &
        (aveSignals['x'] >= -200) &
        (aveSignals['x'] <= 200)
    ]['average']

    c = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['signal'] == "NexteraBackground") &
        (aveSignals['x'] >= -200) &
        (aveSignals['x'] <= 200)
    ]['average']

    atac = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['signal'] == "K562_50K_ATAC_nan_nan_nan_0_0_hg19") &
        (aveSignals['x'] >= -200) &
        (aveSignals['x'] <= 200)
    ]['average']

    norm = smooth(t, 15)[15:] / np.exp(smooth(zscore(c), 15)[15:])
    atacnorm = smooth(atac, 15)[15:] / np.exp(smooth(zscore(c), 15)[15:])

    plt.plot(range(-200, 200), smooth(zscore(c), 15)[15:], label="Nextera genomic DNA library")
    # plt.plot(range(-200, 200), smooth(zscore(t), 15)[15:], label="ChIPmentation")
    plt.plot(range(-200, 200), zscore(norm), label="ChIPmentation normalized")
    # plt.plot(range(-200, 200), smooth(zscore(atac), 15)[15:], label="ATAC-seq")
    plt.plot(range(-200, 200), zscore(atacnorm), label="ATAC-seq normalized")
    plt.legend()
    plt.savefig(os.path.join(plotsDir, sample + ".all.footprints.norm.pdf"), bbox_inches='tight')
    plt.close()



grid = sns.FacetGrid(sub, col="sample", hue="signal", sharey=False, col_wrap=4)
grid.map(plt.plot, "x", "average")
grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
grid.add_legend()
plt.savefig(os.path.join(plotsDir, "footprints.raw.pdf"), bbox_inches='tight')

sub = aveSignals[
    (aveSignals['type'] == "rawSmooth") &
    (aveSignals['x'] >= -400) &
    (aveSignals['x'] <= 400)
]

grid = sns.FacetGrid(sub, col="sample", hue="signal", sharey=False, col_wrap=4)
grid.map(plt.plot, "x", "average")
grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
grid.add_legend()
plt.savefig(os.path.join(plotsDir, "footprints.rawSmooth.pdf"), bbox_inches='tight')


# individual
for sample in aveSignals['sample'].unique():
    sub = aveSignals[
        (aveSignals['sample'] == "K562_10M_CM_CTCF_nan_nan_0_0_hg19") &
        (aveSignals['type'] == "raw") &
        (aveSignals['x'] >= -400) &
        (aveSignals['x'] <= 400)
    ]

    # plot
    grid = sns.FacetGrid(sub, col="signal", hue="signal", sharey=False, col_wrap=4)
    grid.map(plt.plot, "x", "average")
    # grid.set(xlim=(-400, 400))
    grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(os.path.join(plotsDir, sample + ".footprints.raw.pdf"), bbox_inches='tight')


# individual one plot
for sample in aveSignals['sample'].unique():
    fig = plt.figure()
    for signal in ['DNase', sample, 'K562_50K_ATAC_nan_nan_nan_0_0_hg19']:
        sub = aveSignals[
            (aveSignals['sample'] == sample) &
            (aveSignals['signal'] == signal + "-normalized") &
            (aveSignals['type'] == "raw") &
            (aveSignals['x'] >= -400) &
            (aveSignals['x'] <= 400)
        ]
        if len(sub) > 0:
            plt.plot(zscore(smooth(sub['average'])), label=signal)
    plt.legend()


for sample in aveSignals['sample'].unique():
    sub = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['type'] == "rawSmooth") &
        (aveSignals['x'] >= -400) &
        (aveSignals['x'] <= 400)
    ]

    # plot
    grid = sns.FacetGrid(sub, col="signal", hue="signal", sharey=False, col_wrap=4)
    grid.map(plt.plot, "x", "average")
    # grid.set(xlim=(-400, 400))
    grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(os.path.join(plotsDir, sample + ".footprints.rawSmooth.pdf"), bbox_inches='tight')


# several TF footprints
sub = aveSignals[
    (aveSignals['sample'].str.contains("CTCF|GATA1|PU1|REST")) &
    (aveSignals['type'] == "raw") &
    (aveSignals['x'] >= -400) &
    (aveSignals['x'] <= 400)
]

grid = sns.FacetGrid(sub, col="sample", hue="signal", sharey=False, col_wrap=4)
grid.map(plt.plot, "x", "average")
grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
grid.add_legend()
plt.savefig(os.path.join(plotsDir, "footprints.raw.pdf"), bbox_inches='tight')

sub = aveSignals[
    (aveSignals['sample'].str.contains("CTCF|GATA1|PU1|REST")) &
    (aveSignals['signal'].str.contains("CTCF|GATA1|PU1|REST")) &
    (aveSignals['type'] == "rawSmooth") &
    (aveSignals['x'] >= -400) &
    (aveSignals['x'] <= 400)
]

grid = sns.FacetGrid(sub, col="sample", hue="signal", sharey=False, col_wrap=4)
grid.map(plt.plot, "x", "average")
grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
grid.add_legend()
plt.savefig(os.path.join(plotsDir, "footprints.rawSmooth.pdf"), bbox_inches='tight')


# Fig 3:
# compare footprints of various techniques
# plot large area with smoothed signal
aveSignals['ip'] = aveSignals.apply(lambda x: x['sample'].split("_")[3], axis=1)

fig, axis = plt.subplots(2, 2, sharex=True, figsize=(10, 8))
for i in range(len(aveSignals.ip.unique())):
    ip = aveSignals.ip.unique()[i]
    sub = aveSignals[
        (aveSignals['ip'] == ip) &
        (aveSignals['sample'].str.contains("10M")) &
        (aveSignals['type'] == "rawSmooth") &
        (aveSignals['x'] >= -400) &
        (aveSignals['x'] <= 400)
    ]
    sample = sub.sample.unique()[0]

    if i is 0:
        j, k = (0, 0)
    elif i is 1:
        j, k = (0, 1)
    elif i is 2:
        j, k = (1, 0)
    elif i is 3:
        j, k = (1, 1)

    for signal in [sample, "DNase", "ATAC"]:
        if signal == sample:
            sub2 = sub[sub['signal'] == sample]
            colour = colourPerFactor(sample.split("_")[3])
            axis[j][k].plot(sub2.average, color=colour, linestyle="-", label="ChIPmentation", alpha=1)
        else:
            sub2 = sub[sub['signal'].str.contains(signal)]
            colour = colourPerFactor(signal)
            axis[j][k].plot(sub2.average, color=colour, linestyle="-", label=signal, alpha=0.35)

    # subplot attributes
    axis[j][k].set_title(ip)
    axis[j][k].legend(loc="best")
    axis[j][k].set_xlabel("Distance to motif")
plt.savefig(os.path.join(plotsDir, "footprints.rawSmooth.signals.pdf"), bbox_inches='tight')
plt.close()


# compare footprints of various techniques
# plot small area with raw signal
fig, axis = plt.subplots(2, 2, sharex=True, figsize=(10, 8))
for i in range(len(aveSignals.ip.unique())):
    ip = aveSignals.ip.unique()[i]
    sub = aveSignals[
        (aveSignals['ip'] == ip) &
        (aveSignals['sample'].str.contains("10M")) &
        (aveSignals['type'] == "rawSmooth") &
        (aveSignals['x'] >= -100) &
        (aveSignals['x'] <= 100)
    ]
    sample = sub.sample.unique()[0]

    if i is 0:
        j, k = (0, 0)
    elif i is 1:
        j, k = (0, 1)
    elif i is 2:
        j, k = (1, 0)
    elif i is 3:
        j, k = (1, 1)

    for signal in [sample, "DNase", "ATAC"]:
        if signal == sample:
            sub2 = sub[sub['signal'] == sample]
            colour = colourPerFactor(sample.split("_")[3])
            axis[j][k].plot(sub2.average, color=colour, linestyle="-", label="ChIPmentation", alpha=1)
        else:
            sub2 = sub[sub['signal'].str.contains(signal)]
            colour = colourPerFactor(signal)
            axis[j][k].plot(sub2.average, color=colour, linestyle="-", label=signal, alpha=0.35)

    # subplot attributes
    axis[j][k].set_title(ip)
    axis[j][k].legend(loc="best")
    axis[j][k].set_xlabel("Distance to motif")
plt.savefig(os.path.join(plotsDir, "footprints.rawSmooth.signals.zoom.pdf"), bbox_inches='tight')
plt.close()


# Supplement:
# compare footprints of different cell numbers for the same factor
# plot large area with smoothed signal
aveSignals['ip'] = aveSignals.apply(lambda x: x['sample'].split("_")[3], axis=1)
sub = aveSignals[aveSignals['sample'].str.contains("CTCF|GATA1|PU1|REST")]

fig, axis = plt.subplots(2, 2, sharex=True, figsize=(10, 8))

for i in range(len(aveSignals.ip.unique())):
    ip = aveSignals.ip.unique()[i]
    sub = aveSignals[aveSignals['ip'] == ip]

    if i is 0:
        j, k = (0, 0)
    elif i is 1:
        j, k = (0, 1)
    elif i is 2:
        j, k = (1, 0)
    elif i is 3:
        j, k = (1, 1)

    for sample in sub.sample.unique():
        sub = aveSignals[
            (aveSignals['ip'] == ip) &
            (aveSignals['sample'] == sample) &
            (aveSignals['signal'] == sample) &
            (aveSignals['type'] == "rawSmooth") &
            (aveSignals['x'] >= -400) &
            (aveSignals['x'] <= 400)
        ]
        colour = colourPerFactor(sample.split("_")[3])
        if "100K" in sample:
            number = "100K"
        elif "10M" in sample:
            number = "10M"
        elif "500K" in sample:
            number = "500K"

        if number is not "10M":
            axis[j][k].plot(sub.average, color=colour, linestyle="-", label=number, alpha=0.7)
        else:
            axis[j][k].plot(sub.average, color=colour, linestyle="-", label=number, alpha=1)
    # subplot attributes
    axis[j][k].set_title(ip)
    axis[j][k].legend()
    axis[j][k].set_xlabel("Distance to motif")
plt.savefig(os.path.join(plotsDir, "footprints.rawSmooth.sameIP.pdf"), bbox_inches='tight')
plt.close()


# compare footprints of different cell numbers for the same factor
# plot small area with raw signal
fig, axis = plt.subplots(2, 2, sharex=True, figsize=(10, 8))
for i in range(len(aveSignals.ip.unique())):
    ip = aveSignals.ip.unique()[i]
    sub = aveSignals[aveSignals['ip'] == ip]

    if i is 0:
        j, k = (0, 0)
    elif i is 1:
        j, k = (0, 1)
    elif i is 2:
        j, k = (1, 0)
    elif i is 3:
        j, k = (1, 1)

    for sample in sub.sample.unique():
        sub = aveSignals[
            (aveSignals['ip'] == ip) &
            (aveSignals['sample'] == sample) &
            (aveSignals['signal'] == sample) &
            (aveSignals['type'] == "raw") &
            (aveSignals['x'] >= -100) &
            (aveSignals['x'] <= 100)
        ]
        colour = colourPerFactor(sample.split("_")[3])
        if "100K" in sample:
            number = "100K"
        elif "10M" in sample:
            number = "10M"
        elif "500K" in sample:
            number = "500K"

        if number is not "10M":
            axis[j][k].plot(sub.average, color=colour, linestyle="-", label=number, alpha=0.7)
        else:
            axis[j][k].plot(sub.average, color=colour, linestyle="-", label=number, alpha=1)
    # subplot attributes
    axis[j][k].set_title(ip)
    axis[j][k].legend()
    axis[j][k].set_xlabel("Distance to motif")
plt.savefig(os.path.join(plotsDir, "footprints.raw.sameIP.zoom.pdf"), bbox_inches='tight')
plt.close()


# Heatmaps sorted by signal abundance
for sampleName in sampleSubset['sampleName'].unique():
    # Self
    signalName = sampleName
    exportName = "-".join([sampleName, signalName])
    print(exportName)

    if os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
        # if not os.path.isfile(os.path.join(plotsDir, "cdt", exportName + ".cdt")):
        try:
            df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()
        except:
            print("Couldn't load sample %s" % sampleName)
            continue

        pos = df.ix[range(0, len(df), 2)].reset_index(drop=True)  # positive strand
        neg = df.ix[range(1, len(df), 2)].reset_index(drop=True)  # negative strand
        df = (pos + neg) / 2

        df = df[range(-400, 400)]

        s = df.apply(sum, axis=1)
        df["s"] = s.tolist()

        df.sort(['s'], ascending=False, inplace=True)
        df.drop('s', axis=1, inplace=True)
        exportToJavaTreeView(df, os.path.join(plotsDir, "cdt", exportName + ".400bp.cdt"))
    else:
        print("Sample does not have pdy. %s" % exportName)
