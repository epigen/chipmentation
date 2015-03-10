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
from collections import OrderedDict
import re
import HTSeq
import pybedtools
import numpy as np
import pandas as pd
import gc

from divideAndSlurm import DivideAndSlurm, Task
import string, time, random, textwrap

import matplotlib.pyplot as plt
import seaborn as sns

from ggplot import ggplot, aes, stat_smooth, facet_wrap, xlab, ylab, theme_bw, theme

import rpy2.robjects as robj # for ggplot in R
import rpy2.robjects.pandas2ri # for R dataframe conversion
from rpy2.robjects.packages import importr

from sklearn.cluster import k_means
import cPickle as pickle

import plotly.plotly as plotly
from plotly.graph_objs import Data, Heatmap, Layout, Figure

import numpy.lib
import numpy as np
import pandas as pd
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

    s = numpy.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='valid')
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
        Same as numpy.load option
    '''
    values = np.load(fname, mmap_mode=mmap_mode)
    with open(fname) as f:
        numpy.lib.format.read_magic(f)
        numpy.lib.format.read_array_header_1_0(f)
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
                # Activate virtual environment
                source /home/arendeiro/venv/bin/activate

                python /home/arendeiro/coverage_parallel.py {0} {1} {2} """.format(inputPickle, outputPickle, self.bam_file)

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
projectRoot = "/fhgfs/groups/lab_bock/shared/projects/cm/"
resultsDir = projectRoot + "results"
plotsDir = resultsDir + "/plots"
bedFilePath = "/home/arendeiro/reference/Homo_sapiens/hg19.cage_peak_coord_robust.TATA_Annotated.bed"

# Get samples
samples = pd.read_csv(os.path.abspath(projectRoot + "cm.replicates.annotation_sheet.csv"))

# subset samples
sampleSubset = samples[samples['sampleName'].str.contains(
    "K562_10M_CM_GATA1_nan_nan_1_0_hg19|" +
    "K562_10M_CM_PU1_nan_nan_0_0_hg19|" +
    "K562_10M_CM_REST_nan_nan_1_0_hg19"
)].reset_index(drop=True)

# subset samples into signals
signals = samples[samples['sampleName'].str.contains(
    "K562_10M_CHIP_H3K27AC_nan_nan_1_0_hg19|" +
    "K562_10M_CHIP_H3K4ME1_nan_nan_1_0_hg19|" +
    "K562_10M_CM_H3K27AC_nan_nan_1_0_hg19|" +
    "K562_10M_CM_H3K4ME1_nan_nan_1_0_hg19|" +
    "K562_500K_CHIP_H3K27ME3_nan_nan_0_0_hg19|" +
    "K562_500K_CHIP_H3K4ME3_nan_nan_0_0_hg19|" +
    "K562_500K_CM_H3K27ME3_nan_nan_0_0_hg19|" +
    "K562_500K_CM_H3K4ME3_nan_nan_0_0_hg19"
)].reset_index(drop=True)

genome = "hg19"
windowRange = (-1000, 1000)
fragmentsize = 1
duplicates = False
n_clusters = 5

plotly.sign_in("afrendeiro", "iixmygxac1")
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
    except IOError("File not found"), e:
        raise e

    # exlude peaks in gaps or repeats
    peaks.intersect(b=gapsRepeats, v=True, wa=True)

    peaks = bedToolsInterval2GenomicInterval(peaks)
    # Filter peaks near chrm borders
    for peak, interval in peaks.iteritems():
        if interval.length < windowWidth:
            peaks.pop(peak)

    # Get self coverage
    signalName = sampleName
    exportName = "-".join([sampleName, signalName])
    bamFile = sampleSubset['filePath'][i]

    task = Coverage(peaks, 2, os.path.join(bamFile),
                    orientation=True, fragment_size=1, strand_wise=True, queue="shortq", cpusPerTask=8
                    )
    task.id = exportName
    slurm.add_task(task)
    slurm.submit(task)

    # Get coverage over other signals
    for j in range(len(signals)):
        signalName = signals['sampleName'][j]

        exportName = "-".join([sampleName, signalName])
        if os.path.isfile(os.path.join(resultsDir, exportName + ".pdy")):
            continue
        bamFile = signals['filePath'][j]

        # Get dataframe of signal coverage in bed regions, append to dict
        # Make task
        task = Coverage(peaks, 2, os.path.join(bamFile),
                        orientation=True, fragment_size=1, strand_wise=True, queue="shortq", cpusPerTask=8
                        )
        task.id = exportName

        slurm.add_task(task)
        slurm.submit(task)

while not all([t.is_ready() for t in slurm.tasks]):
    time.sleep(10)

# Collect task outputs
for task in slurm.tasks:
    print task.id
    if os.path.isfile(os.path.join(resultsDir, task.id + ".pdy")):
        continue
    if not task.is_ready():
        continue
    cov = task.collect()

    # Make multiindex dataframe
    levels = [cov.keys(), ["+", "-"]]
    labels = [[y for x in range(len(cov)) for y in [x, x]], [y for x in range(len(cov.keys())) for y in (0, 1)]]
    index = pd.MultiIndex(labels=labels, levels=levels, names=["peak", "strand"])
    df = pd.DataFrame(np.vstack(cov.values()), index=index)
    df.columns = range(windowRange[0], windowRange[1])

    # Save raw data
    savePandas(os.path.join(resultsDir, task.id + ".pdy"), df)
    del cov
    del df
    gc.collect()

# Average signals
aveSignals = dict()
# aveSignals = pickle.load(open(os.path.join(resultsDir, "aveSignals.pickle"), "r"))
for i in range(len(sampleSubset)):
    sampleName = sampleSubset['sampleName'][i]

    signalName = sampleName
    exportName = "-".join([sampleName, signalName])
    if exportName in aveSignals.keys():
        continue

    df = loadPandas(os.path.join(resultsDir, exportName + ".pdy")).copy()
    # Get average profiles and append to dict
    ave = {"sample": sampleName,
           "signal": signalName,
           "average": pd.Series(smooth(df.apply(np.mean, axis=0))[10:], index=df.columns),                           # both strands
           "positive": pd.Series(smooth(df.ix[range(0, len(df), 2)].apply(np.mean, axis=0))[10:], index=df.columns), # positive strand
           "negative": pd.Series(smooth(df.ix[range(1, len(df), 2)].apply(np.mean, axis=0))[10:], index=df.columns)  # negative strand
           }
    aveSignals[exportName] = pd.DataFrame(ave)
    pickle.dump(aveSignals, open(os.path.join(resultsDir, "aveSignals.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    for j in range(len(signals)):
        signalName = signals['sampleName'][j]

        exportName = "-".join([sampleName, signalName])

        if exportName in aveSignals.keys():
            continue

        df = loadPandas(os.path.join(resultsDir, exportName + ".pdy")).copy()
        # Get average profiles and append to dict
        ave = {"sample": sampleName,
               "signal": signalName,
               "average": pd.Series(smooth(df.apply(np.mean, axis=0))[10:], index=df.columns),                           # both strands
               "positive": pd.Series(smooth(df.ix[range(0, len(df), 2)].apply(np.mean, axis=0))[10:], index=df.columns), # positive strand
               "negative": pd.Series(smooth(df.ix[range(1, len(df), 2)].apply(np.mean, axis=0))[10:], index=df.columns)  # negative strand
               }
        aveSignals[exportName] = pd.DataFrame(ave)

        pickle.dump(aveSignals, open(os.path.join(resultsDir, "aveSignals.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

# Plot average profiles for all signals with ggplot
aveSignals = pd.concat(aveSignals, axis=0)
aveSignals["x"] = [i[1] for i in aveSignals.index]
aveSignals = pd.melt(aveSignals, id_vars=["x", "sample", "signal"])

plotFunc = robj.r("""
    library(ggplot2)

    function(df, width, plotsDir){
        p = ggplot(df, aes(x, value, colour = variable)) +
            geom_line() +
            #stat_smooth(se = F, method = "lm", formula = y ~ poly(x, 24)) +
            #stat_smooth(method = "gam", formula = y ~ s(x, k = 80), se = FALSE) +
            #stat_smooth(method = "loess", formula = y ~ x, size = 1, se = FALSE) +
            facet_wrap(sample~signal, scales="free", ncol=9) +
            xlab("Distance to peak") +
            ylab("Average tags per bp") +
            theme_bw() +
            theme(legend.title=element_blank()) +
            scale_colour_manual(values=c("black", "dark blue","dark red")) +
            xlim(c(-width, width))

        ggsave(filename = paste0(plotsDir, "/TFs_signal.", width * 2, ".pdf"), plot = p, height = 20, width = 40)
        ggsave(filename = paste0(plotsDir, "/TFs_signal.", width * 2, ".wide.pdf"), plot = p, height = 20, width = 60, limitsize=FALSE)
    }
""")

gr = importr('grDevices')
# convert the pandas dataframe to an R dataframe
robj.pandas2ri.activate()
aveSignals_R = robj.conversion.py2ri(aveSignals)

# run the plot function on the dataframe
plotFunc(aveSignals_R, 50, plotsDir)
plotFunc(aveSignals_R, 100, plotsDir)
plotFunc(aveSignals_R, 200, plotsDir)
plotFunc(aveSignals_R, 500, plotsDir)










# Clustering

# Loop through raw signals, normalize, k-means cluster, save pickle, plot heatmap, export cdt
for signal in signals:
    exportName = "{0}.tssSignal_{1}bp.kmeans_{2}k".format(signal, str(windowWidth), n_clusters)

    df = rawSignals[signal]

    # join strand signal (plus - minus)
    df = df.xs('+', level="strand") - df.xs('-', level="strand")

    # scale row signal to 0:1 (normalization)
    dfNorm = df.apply(lambda x: (x - min(x)) / (max(x) - min(x)), axis=1)

    clust = k_means(dfNorm,
                    n_clusters,
                    n_init=25,
                    max_iter=10000,
                    n_jobs=2
                    )  # returns centroid, label, inertia

    # save object
    pickle.dump(clust,
                open(os.path.join(bamFilePath, exportName + ".pickl"), "wb"),
                protocol=pickle.HIGHEST_PROTOCOL
                )

    # Sort dataframe by cluster order
    dfNorm["cluster"] = clust[1]  # clust[1] <- label from k-means clustering
    dfNorm.sort_index(by="cluster", axis=0, inplace=True)
    dfNorm.drop("cluster", axis=1, inplace=True)

    # Plot heatmap
    data = Data([Heatmap(z=np.array(dfNorm), colorscale='Portland')])
    plotly.image.save_as(data, os.path.join(plotsDir, exportName + ".pdf"))

    # Export as cdt
    exportToJavaTreeView(df, os.path.join(plotsDir, exportName + ".cdt"))











"""
Investigating the contribution of chromatin loops to TF ChIPmentation signal
"""
slurm = DivideAndSlurm()
signal = "PU1_K562_10mio_CM"
sample = "CTCF_K562_10mio_CM"
sample = "ASP14_50k_ATAC-seq_nan_nan_untreated_ATAC10-7_0_0.trimmed.bowtie2.shifted.dups"
sample = "CTCF_K562_10mio_ChIP"

# load ChIPmentation peaks
peaks = pybedtools.BedTool(os.path.join(peakFilePath, signal + ".motifStrand.bed"))

# intersect with chromatin loops
loops = pybedtools.BedTool("reference/Homo_sapiens/Hi-C/GSE63525_K562_HiCCUPS_looplist.concat.txt")

peaks_in_loops = peaks.intersect(b=loops, wa=True)
peaks_notin_loops = peaks.intersect(b=loops, v=True, wa=True)

peaks_in_loops = bedToolsInterval2GenomicInterval(peaks_in_loops)
peaks_notin_loops = bedToolsInterval2GenomicInterval(peaks_notin_loops)


for name, interval in peaks_in_loops.iteritems():
    if interval.length < windowWidth:
        peaks_in_loops.pop(name)
for name, interval in peaks_notin_loops.iteritems():
    if interval.length < windowWidth:
        peaks_notin_loops.pop(name)

loopsTask = Coverage(peaks_in_loops, 2, os.path.join(bamFilePath, sample + ".bam"),
                     orientation=True, fragment_size=1, strand_wise=True, queue="develop", cpusPerTask=4
                     )
notLoopsTask = Coverage(peaks_notin_loops, 2, os.path.join(bamFilePath, sample + ".bam"),
                        orientation=True, fragment_size=1, strand_wise=True, queue="develop", cpusPerTask=4
                        )
slurm.add_task(loopsTask)
slurm.add_task(notLoopsTask)

slurm.submit(loopsTask)
slurm.submit(notLoopsTask)

loopsCov = loopsTask.collect()
notLoopsCov = notLoopsTask.collect()

levels = [loopsCov.keys(), ["+", "-"]]  # for strand_wise=True
labels = [[y for x in range(len(loopsCov)) for y in [x, x]], [y for x in range(len(loopsCov.keys())) for y in (0, 1)]]
index = pd.MultiIndex(labels=labels, levels=levels, names=["peak", "strand"])
loopsDf = pd.DataFrame(np.vstack(loopsCov.values()), index=index)
loopsDf.columns = range(windowRange[0], windowRange[1])

levels = [notLoopsCov.keys(), ["+", "-"]]
labels = [[y for x in range(len(notLoopsCov)) for y in [x, x]], [y for x in range(len(notLoopsCov.keys())) for y in (0, 1)]]
index = pd.MultiIndex(labels=labels, levels=levels, names=["peak", "strand"])
notLoopsDf = pd.DataFrame(np.vstack(notLoopsCov.values()), index=index)
notLoopsDf.columns = range(windowRange[0], windowRange[1])

# Save as csv
loopsDf.to_csv(os.path.join(coverageFilePath, "{0}_inLoops.tssSignal_{1}bp.csv".format(sample, str(windowWidth))), index=False)
notLoopsDf.to_csv(os.path.join(coverageFilePath, "{0}_notinLoops.tssSignal_{1}bp.csv".format(sample, str(windowWidth))), index=False)
# loopsDf = pd.read_csv(os.path.join(coverageFilePath, "{0}_inLoops.tssSignal_{1}bp.csv".format(sample, str(windowWidth))))
# notLoopsDf = pd.read_csv(os.path.join(coverageFilePath, "{0}_notinLoops.tssSignal_{1}bp.csv".format(sample, str(windowWidth))))

# Pickle
pickle.dump(loopsDf, open(os.path.join(coverageFilePath, "{0}_inLoops.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(notLoopsDf, open(os.path.join(coverageFilePath, "{0}_notinLoops.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
# loopsDf = pickle.load(open(os.path.join(coverageFilePath, "{0}_inLoops.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "r"))
# notLoopsDf = pickle.load(open(os.path.join(coverageFilePath, "{0}_notinLoops.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "r"))


# Get average profiles and append to dict
loopsAve = {"signal": sample,
            "average": loopsDf.apply(np.mean, axis=0),                              # both strands
            "positive": loopsDf.ix[range(0, len(loopsDf), 2)].apply(np.mean, axis=0),    # positive strand
            "negative": loopsDf.ix[range(1, len(loopsDf), 2)].apply(np.mean, axis=0)     # negative strand
            }
notLoopsAve = {"signal": sample,
               "average": notLoopsDf.apply(np.mean, axis=0),                              # both strands
               "positive": notLoopsDf.ix[range(0, len(notLoopsDf), 2)].apply(np.mean, axis=0),    # positive strand
               "negative": notLoopsDf.ix[range(1, len(notLoopsDf), 2)].apply(np.mean, axis=0)     # negative strand
               }
pickle.dump(loopsAve, open(os.path.join(coverageFilePath, "{0}_inLoops_average.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(notLoopsAve, open(os.path.join(coverageFilePath, "{0}_notinLoops_average.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
# loopsAve = pickle.load(open(os.path.join(coverageFilePath, "{0}_inLoops_average.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "r"))
# notLoopsAve = pickle.load(open(os.path.join(coverageFilePath, "{0}_notinLoops_average.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "r"))


# Plot
plt.subplot(121)
plt.plot(loopsAve['average'].index, loopsAve['average'], 'black', linewidth=0.7)
plt.subplot(122)
plt.plot(notLoopsAve['average'].index, notLoopsAve['average'], 'black', linewidth=0.7)
plt.savefig(os.path.join(plotsDir, "%s_peaks.loops.pdf" % sample), bbox_inches='tight')
plt.close()

plt.subplot(121)
plt.plot(loopsAve['positive'].index, loopsAve['positive'], 'r', linewidth=0.7)
plt.plot(loopsAve['negative'].index, loopsAve['negative'], 'b', linewidth=0.7)
plt.subplot(122)
plt.plot(notLoopsAve['positive'].index, notLoopsAve['positive'], 'r', linewidth=0.7)
plt.plot(notLoopsAve['negative'].index, notLoopsAve['negative'], 'b', linewidth=0.7)
plt.savefig(os.path.join(plotsDir, "%s_peaks.loops.strandWise.pdf" % sample), bbox_inches='tight')
plt.close()
