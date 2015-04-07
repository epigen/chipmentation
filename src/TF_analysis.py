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
import gc

from divideAndSlurm import DivideAndSlurm, Task
import string
import time
import random
import textwrap

import matplotlib.pyplot as plt
import seaborn as sns  # changes plt style (+ full plotting library)

import cPickle as pickle


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
                python2.7 /home/arendeiro/coverage_parallel.py {0} {1} {2} """.format(inputPickle, outputPickle, self.bam_file)

            if self.strand_wise:
                task += "--strand-wise "
            if self.duplicates:
                task += "--duplicates "
            if self.orientation:
                task += "--orientation "
            if self.permute:
                task += "--permute "
            task += """--fragment-size {0}
            """.format(self.fragment_size)

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


# Define paths
projectRoot = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/"
resultsDir = projectRoot + "results"
plotsDir = resultsDir + "/plots/"
DNase = "/home/arendeiro/data/human/encode/wgEncodeUwDnaseK562Aln.merged.bam"
MNase = "/home/arendeiro/data/human/encode/wgEncodeSydhNsomeK562AlnRep1.bam"

# Get samples
samples = pd.read_csv(os.path.abspath(projectRoot + "chipmentation.replicates.annotation_sheet.csv"))

# Replace input sample
samples.replace(value="K562_10M_CM_IGG_nan_nan_1_0_hg19", to_replace="K562_10M_CM_IGG_nan_nan_0_0_hg19", inplace=True)
samples.replace(
    value="/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/K562_10M_CM_IGG_nan_nan_1_0_hg19.trimmed.bowtie2.shifted.dups.bam",
    to_replace="/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/K562_10M_CM_IGG_nan_nan_0_0_hg19.trimmed.bowtie2.shifted.dups.bam",
    inplace=True
)

# subset samples
sampleSubset = samples[
    (samples["technique"] == "CM") &
    (samples["ip"].str.contains("CTCF|GATA1|PU1|REST")) &
    (samples["biologicalReplicate"] == 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"])

# subset samples into signals
signals = samples[
    samples["sampleName"].str.contains(
        "K562_10M_CHIP_H3K27AC_nan_nan_1_0_hg19|" +
        "K562_10M_CM_H3K27AC_nan_nan_1_0_hg19|" +
        "K562_10M_CHIP_H3K4ME1_nan_nan_1_0_hg19|" +
        "K562_10M_CM_H3K4ME1_nan_nan_0_0_hg19|" +
        "K562_500K_CHIP_H3K27ME3_nan_nan_0_0_hg19|" +
        "K562_500K_CM_H3K27ME3_nan_nan_0_0_hg19|" +
        "K562_500K_CHIP_H3K4ME3_nan_nan_0_0_hg19|" +
        "K562_500K_CM_H3K4ME3_nan_nan_0_0_hg19|" +
        "K562_50K_ATAC_nan_nan_nan_0_0_hg19"  # +
        # "K562_500K_ATAC_INPUT_nan_0.1ULTN5_PE_CM25_1_1"
    )
].reset_index(drop=True)
signals = signals.sort(["ip", "technique"])

signals = signals.append(pd.Series(data=["DNase", DNase], index=["sampleName", "filePath"]), ignore_index=True)
signals = signals.append(pd.Series(data=["MNase", MNase], index=["sampleName", "filePath"]), ignore_index=True)

windowRange = (-1000, 1000)
fragmentsize = 1
duplicates = False
n_clusters = 5

windowWidth = abs(windowRange[0]) + abs(windowRange[1])

slurm = DivideAndSlurm()

gapsRepeats = pybedtools.BedTool("/home/arendeiro/reference/Homo_sapiens/hg19_gapsRepeats.bed")

# Loop through all samples, submit tasks to compute coverage in peak regions centered on motifs
for i in range(len(sampleSubset)):
    sampleName = sampleSubset['sampleName'][i]
    motifCentered = re.sub("\..*", ".", sampleSubset['peakFile'][i]) + "motifCentered.bed"

    # Load peak file from bed files centered on motif, make window around
    try:
        peaks = pybedtools.BedTool(motifCentered)  # .slop(genome=genome, b=windowWidth / 2)
    except ValueError("File not found"), e:
        raise e

    # Exclude peaks in gaps or repeats
    peaks.intersect(b=gapsRepeats, v=True, wa=True)
    peaks = bedToolsInterval2GenomicInterval(peaks)

    # Filter peaks near chrm borders
    for peak, interval in peaks.iteritems():
        if interval.length < windowWidth:
            peaks.pop(peak)

    # Get self coverage
    signalName = sampleName
    exportName = "-".join([sampleName, signalName])
    if not os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
        print(sampleName, signalName)
        bamFile = sampleSubset['filePath'][i]

        task = Coverage(peaks, 4, os.path.join(bamFile),
                        orientation=True, fragment_size=1,
                        strand_wise=True, queue="shortq", cpusPerTask=8,
                        permissive=True
                        )
        task.id = exportName
        slurm.add_task(task)
        slurm.submit(task)

    # PU1 chip-tagmentation
    if sampleName == "K562_10M_CM_PU1_nan_nan_0_0_hg19":
        signalName = "K562_10M_ATAC_PU1_nan_PE_1_1_hg19"
        exportName = "-".join([sampleName, signalName])
        print(exportName)

        if not os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
            print(sampleName, signalName)
            bamFile = sampleSubset['filePath'][i]

            task = Coverage(peaks, 4, os.path.join(bamFile),
                            orientation=True, fragment_size=1,
                            strand_wise=True, queue="shortq", cpusPerTask=8,
                            permissive=True
                            )
            task.id = exportName
            slurm.add_task(task)
            slurm.submit(task)

    # Get control coverage
    signalName = sampleSubset['controlSampleName'][i]
    exportName = "-".join([sampleName, signalName])
    if not os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
        print(sampleName, signalName)
        bamFile = sampleSubset['controlSampleFilePath'][i]

        task = Coverage(peaks, 4, os.path.join(bamFile),
                        orientation=True, fragment_size=1,
                        strand_wise=True, queue="shortq", cpusPerTask=8,
                        permissive=True
                        )
        task.id = exportName
        slurm.add_task(task)
        slurm.submit(task)

    # Get coverage over other signals
    for j in range(len(signals)):
        signalName = signals['sampleName'][j]
        print(sampleName, signalName)

        exportName = "-".join([sampleName, signalName])
        if not os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
            bamFile = signals['filePath'][j]

            # Make task
            task = Coverage(peaks, 4, os.path.join(bamFile),
                            orientation=True, fragment_size=1,
                            strand_wise=True, queue="shortq", cpusPerTask=8,
                            permissive=True
                            )
            task.id = exportName

            slurm.add_task(task)
            slurm.submit(task)

while not all([t.is_ready() for t in slurm.tasks]):
    time.sleep(10)

# Collect task outputs
for task in slurm.tasks:
    print task.id
    if not os.path.isfile(os.path.join(resultsDir, task.id + ".pdy")) and task.is_ready():
        cov = task.collect()

        # Make multiindex dataframe
        levels = [cov.keys(), ["+", "-"]]
        labels = [[y for x in range(len(cov)) for y in [x, x]], [y for x in range(len(cov.keys())) for y in (0, 1)]]
        index = pd.MultiIndex(labels=labels, levels=levels, names=["peak", "strand"])
        df = pd.DataFrame(np.vstack(cov.values()), index=index)
        df.columns = range(windowRange[0], windowRange[1])

        # Save raw data
        savePandas(os.path.join(plotsDir, "pickles", task.id + ".pdy"), df)
        del cov
        del df
        gc.collect()
    else:
        slurm.tasks.pop(slurm.tasks.index(task))
        del task
        gc.collect()

# Average signals
if os.path.exists(os.path.join(plotsDir, "pickles", "aveSignals.pickle")):
    aveSignals = pickle.load(open(os.path.join(plotsDir, "pickles", "aveSignals.pickle"), "r"))
else:
    aveSignals = pd.DataFrame(columns=["sample", "signal"])
names = aveSignals[["sample", "signal"]].apply("-".join, axis=1).unique()

for i in range(len(sampleSubset)):
    sampleName = sampleSubset['sampleName'][i]

    # Self
    signalName = sampleName
    exportName = "-".join([sampleName, signalName])
    print(exportName)

    if exportName not in names:
        if os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
            df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()

            # Get self average profiles and append to dataframe
            df = pd.DataFrame({
                "sample": sampleName,
                "signal": signalName,
                "average": df.apply(np.mean, axis=0),                            # both strands
                "positive": df.ix[range(0, len(df), 2)].apply(np.mean, axis=0),  # positive strand
                "negative": df.ix[range(1, len(df), 2)].apply(np.mean, axis=0),  # negative strand
                "x": df.columns
            })
            aveSignals = pd.concat([aveSignals, df])
            pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # PU1 chip-tagmentation
    if sampleName == "K562_10M_CM_PU1_nan_PE_1_1_hg19":
        signalName = "K562_10M_ATAC_PU1_nan_PE_1_1_hg19"
        exportName = "-".join([sampleName, signalName])
        print(exportName)

        if os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
            df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()

            # Get self average profiles and append to dataframe
            df = pd.DataFrame({
                "sample": sampleName,
                "signal": signalName,
                "average": df.apply(np.mean, axis=0),                            # both strands
                "positive": df.ix[range(0, len(df), 2)].apply(np.mean, axis=0),  # positive strand
                "negative": df.ix[range(1, len(df), 2)].apply(np.mean, axis=0),  # negative strand
                "x": df.columns
            })
            aveSignals = pd.concat([aveSignals, df])
            pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # Control
    signalName = sampleSubset['controlSampleName'][i]
    exportName = "-".join([sampleName, signalName])
    print(exportName)

    if exportName not in names:
        # Get control average profiles and append to dataframe
        if os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
            df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()
            df = pd.DataFrame({
                "sample": sampleName,
                "signal": signalName,
                "average": df.apply(np.mean, axis=0),                            # both strands
                "positive": df.ix[range(0, len(df), 2)].apply(np.mean, axis=0),  # positive strand
                "negative": df.ix[range(1, len(df), 2)].apply(np.mean, axis=0),  # negative strand
                "x": df.columns
            })
            aveSignals = pd.concat([aveSignals, df])
            pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    for j in range(len(signals)):
        signalName = signals['sampleName'][j]
        exportName = "-".join([sampleName, signalName])
        print(exportName)

        if exportName not in names:
            if os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
                df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()
                # Get average profiles and append to dict
                df = pd.DataFrame({
                    "sample": sampleName,
                    "signal": signalName,
                    "average": df.apply(np.mean, axis=0),                            # both strands
                    "positive": df.ix[range(0, len(df), 2)].apply(np.mean, axis=0),  # positive strand
                    "negative": df.ix[range(1, len(df), 2)].apply(np.mean, axis=0),  # negative strand
                    "x": df.columns
                })
                aveSignals = pd.concat([aveSignals, df])
                pickle.dump(aveSignals, open(os.path.join(plotsDir, "pickles", "aveSignals.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

aveSignals = pickle.load(open(os.path.join(plotsDir, "pickles", "aveSignals.pickle"), "r"))
aveSignals.drop_duplicates(inplace=True)
aveSignals["type"] = "raw"
aveSignals.reset_index(drop=True, inplace=True)

# Normalize
df = aveSignals.copy()
df2 = pd.DataFrame()

for sample in df['sample'].unique():
    for signal in df[df['sample'] == sample]['signal'].unique():
        for strand in ["average", "positive", "negative"]:
            treat = df[
                (df["sample"] == sample) &
                (df["signal"] == signal)
            ][strand]

            controlName = sampleSubset.loc[sampleSubset["sampleName"] == sample, "controlSampleName"]
            controlName = controlName.tolist()[0]
            ctrl = df[
                (df["sample"] == sample) &
                (df["signal"] == controlName)
            ][strand]

            # standardize: 0-1
            treat = smooth((treat - treat.min()) / (treat.max() - treat.min()))
            ctrl = smooth((ctrl - ctrl.min()) / (ctrl.max() - ctrl.min()))

            # normalize by input
            norm = np.array(treat) - np.array(ctrl)

            tmp = pd.DataFrame()
            tmp[strand] = norm
            tmp['sample'] = sample
            tmp['signal'] = signal
            tmp['x'] = np.array(range(len(tmp))) - (len(tmp) / 2)
            tmp['type'] = "norm"

            df2 = df2.append(tmp, ignore_index=True)

# append normalized df to df
aveSignals = pd.concat([aveSignals, df2])
aveSignals.reset_index(drop=True, inplace=True)

# Grid plots
# raw
for sample in aveSignals['sample'].unique():
    sub = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['type'] == "raw") &
        (aveSignals['x'] >= -400) &
        (aveSignals['x'] <= 400)
    ]

    # plot
    grid = sns.FacetGrid(sub, col="signal", hue="signal", sharey=False, col_wrap=4)
    grid.map(plt.plot, "x", "positive")
    # grid.set(xlim=(-100, 100))
    grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(os.path.join(plotsDir, sample + ".footprints.raw.pdf"), bbox_inches='tight')

# norm
for sample in aveSignals['sample'].unique():
    sub = aveSignals[
        (aveSignals['sample'] == sample) &
        (aveSignals['type'] == "norm") &
        (aveSignals['x'] >= -400) &
        (aveSignals['x'] <= 400)
    ]

    # plot
    grid = sns.FacetGrid(sub, col="signal", hue="signal", sharey=False, col_wrap=4)
    grid.map(plt.plot, "x", "positive")
    # grid.set(xlim=(-100, 100))
    grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
    plt.savefig(os.path.join(plotsDir, sample + ".footprints.norm.pdf"), bbox_inches='tight')


# Overlaid plots
for sample in aveSignals['sample'].unique():
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')
    for signal in aveSignals[aveSignals['sample'] == sample]['signal'].unique():
        sub = aveSignals[
            (aveSignals['sample'] == sample) &
            (aveSignals['signal'] == signal) &
            (aveSignals['type'] == "raw") &
            (aveSignals['x'] >= -400) &
            (aveSignals['x'] <= 400)
        ]
        # plot raw
        plt.plot(sub["x"], np.log10(sub["average"]), label=signal, color=colourPerFactor(signal))
    plt.title(sample)
    plt.xlabel("distance to motif")
    plt.ylabel("distance to motif")
    plt.legend(loc='best', prop={'size': 4})
    plt.savefig(os.path.join(plotsDir, sample + ".footprints.raw.pdf"), bbox_inches='tight',)
    plt.close()

    plt.figure()
    for signal in aveSignals[aveSignals['sample'] == sample]['signal'].unique():
        # plot norm
        sub = aveSignals[
            (aveSignals['sample'] == sample) &
            (aveSignals['signal'] == signal) &
            (aveSignals['type'] == "norm") &
            (aveSignals['x'] >= -400) &
            (aveSignals['x'] <= 400)
        ]
        # plot norm
        plt.plot(sub["x"], np.log10(sub["average"]), label=signal, color=colourPerFactor(signal))
    plt.title(sample)
    plt.xlim(-400, 400)
    plt.xlabel("distance to motif")
    plt.ylabel("distance to motif")
    plt.legend(loc='best', prop={'size': 4})
    plt.savefig(os.path.join(plotsDir, sample + ".footprints.norm.pdf"), bbox_inches='tight',)

# Heatmaps sorted by signal abundance

# Average signals
aveSignals = pickle.load(open(os.path.join(plotsDir, "pickles", "aveSignals.pickle"), "r"))

for i in range(len(sampleSubset)):
    sampleName = sampleSubset['sampleName'][i]

    # Self
    signalName = sampleName
    exportName = "-".join([sampleName, signalName])
    print(exportName)

    if os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
        df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()

        pos = df.ix[range(0, len(df), 2)].reset_index(drop=True)  # positive strand
        neg = df.ix[range(1, len(df), 2)].reset_index(drop=True)  # negative strand
        df = (pos + neg) / 2

        exportToJavaTreeView(df, os.path.join(plotsDir, "cdt", exportName + ".cdt"))

    # PU1 chip-tagmentation
    if sampleName == "K562_10M_CM_PU1_nan_PE_1_1_hg19":
        signalName = "K562_10M_ATAC_PU1_nan_PE_1_1_hg19"
        exportName = "-".join([sampleName, signalName])
        print(exportName)

        if os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
            df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()

            pos = df.ix[range(0, len(df), 2)].reset_index(drop=True)  # positive strand
            neg = df.ix[range(1, len(df), 2)].reset_index(drop=True)  # negative strand
            df = (pos + neg) / 2

            exportToJavaTreeView(df, os.path.join(plotsDir, "cdt", exportName + ".cdt"))

    # Control
    signalName = sampleSubset['controlSampleName'][i]
    exportName = "-".join([sampleName, signalName])
    print(exportName)

    # Get control average profiles and append to dataframe
    if os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
        df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()

        pos = df.ix[range(0, len(df), 2)].reset_index(drop=True)  # positive strand
        neg = df.ix[range(1, len(df), 2)].reset_index(drop=True)  # negative strand
        df = (pos + neg) / 2

        exportToJavaTreeView(df, os.path.join(plotsDir, "cdt", exportName + ".cdt"))

    for j in range(len(signals)):
        signalName = signals['sampleName'][j]
        exportName = "-".join([sampleName, signalName])
        print(exportName)

        if os.path.isfile(os.path.join(plotsDir, "pickles", exportName + ".pdy")):
            df = loadPandas(os.path.join(plotsDir, "pickles", exportName + ".pdy")).copy()

            pos = df.ix[range(0, len(df), 2)].reset_index(drop=True)  # positive strand
            neg = df.ix[range(1, len(df), 2)].reset_index(drop=True)  # negative strand
            df = (pos + neg) / 2

            exportToJavaTreeView(df, os.path.join(plotsDir, "cdt", exportName + ".cdt"))
