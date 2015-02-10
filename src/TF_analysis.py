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
import HTSeq
import pybedtools
import numpy as np
import pandas as pd

from divideAndSlurm import DivideAndSlurm, Task
import string, time, random, textwrap

from ggplot import ggplot, aes, stat_smooth, facet_wrap, xlab, ylab, theme_bw, theme

import rpy2.robjects as robj # for ggplot in R
import rpy2.robjects.pandas2ri # for R dataframe conversion
from rpy2.robjects.packages import importr

from sklearn.cluster import k_means 
import cPickle as pickle

import plotly.plotly as plotly
from plotly.graph_objs import Data, Heatmap, Layout, Figure

signals = { "CTCF_K562_10mio_CM" : ["CTCF_K562_10mio_CM", "CTCF_K562_10mio_ChIP", "IgG_K562_10mio_CM"],
            "PU1_K562_10mio_CM" : ["PU1_K562_10mio_CM", "PU1_K562_10mio_ChIP", "IgG_K562_10mio_CM"]} # DNase missing
bamFilePath = "/home/arendeiro/data/human/chipmentation/mapped/merged"
peakFilePath = "/home/arendeiro/data/human/chipmentation/bed"
coverageFilePath = "/home/arendeiro/data/human/chipmentation/bed/newPython" # new output dir
plotsDir = "/home/arendeiro/"

genome = "hg19"
windowRange = (-2000, 2000)
fragmentsize = 1
duplicates = True
n_clusters = 5

plotly.sign_in("afrendeiro", "iixmygxac1")
windowWidth = abs(windowRange[0]) + abs(windowRange[1])



class Coverage(Task):
    """
    Task to get read coverage under regions.
    """
    def __init__(self, data, fractions, *args, **kwargs):
        super(Coverage, self).__init__(data, fractions, *args, **kwargs)
        ### Initialize rest
        now = string.join([time.strftime("%Y%m%d%H%M%S", time.localtime()), str(random.randint(1,1000))], sep="_")
        self.name = "coverage_{0}".format(now)

        ### Parse
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
        self.log = os.path.join(self.slurm.logDir, string.join([self.name, "log"], sep=".")) # add abspath

        ### Split data in fractions
        ids, groups, files = self._split_data()

        ### Make jobs with groups of data
        self.jobs = list(); self.jobFiles = list(); self.inputPickles = list(); self.outputPickles = list()

        # for each group of data
        for i in xrange(len(ids)):
            jobFile = files[i] + "_coverage.sh"
            inputPickle = files[i] + ".input.pickle"
            outputPickle = files[i] + ".output.pickle"

            ### assemble job file
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
        if not hasattr(self, "output"): # if output is already stored, just return it
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
                        self.output = output # store output in object
                        self._rm_temps() # delete tmp files
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



slurm = DivideAndSlurm()

# Loop through all samples, compute coverage in peak regions centered on motifs,
# save dicts with coverage and average profiles
rawSignals = OrderedDict()
aveSignals = OrderedDict()
for signal in signals.keys():
    # Load peak file from bed files, centered on motif

    ## TODO: modify pipeline so that slop is done here
    # peaks = pybedtools.BedTool(os.path.join(peakFilePath, signal + ".motifStrand.bed")).slop(genome=genome, b=windowWidth/2)
    peaks = pybedtools.BedTool(os.path.join(peakFilePath, signal + ".motifStrand.bed"))
    peaks = bedToolsInterval2GenomicInterval(peaks)
    # Filter peaks near chrm borders
    for name, interval in peaks.iteritems():
        if interval.length < windowWidth:
            peaks.pop(name)

    for sample in signals[signal]:
        # Get dataframe of signal coverage in bed regions, append to dict
        # Make task
        task = Coverage(peaks, 30, os.path.join(bamFilePath, sample + ".bam"),
            orientation=True, fragment_size=1, strand_specific=False, queue="shortq", cpusPerTask=8
        )
        slurm.add_task(task)
        slurm.submit(task)
        while not task.is_ready():
            time.sleep(10)
        cov = task.collect()

        # Make multiindex dataframe
        levels = [cov.keys(), ["+", "-"]]
        labels = [[y for x in range(len(cov)) for y in [x, x]], [y for x in range(len(cov.keys())) for y in (0,1)]]
        index = pd.MultiIndex(labels=labels, levels=levels, names=["peak", "strand"])
        df = pd.DataFrame(np.vstack(cov.values()), index=index)
        df.columns = range(windowRange[0], windowRange[1])

        # For strand_specific=False
        #cov = coverage(bamfile, peaks, fragmentsize)
        #df = pd.DataFrame(cov).T
        #df.columns = range(windowRange[0], windowRange[1])

        # append to dict   
        rawSignals[sample] = df

        # Save as csv
        df.to_csv(os.path.join(coverageFilePath, "{0}.tssSignal_{1}bp.csv".format(sample, str(windowWidth))), index=False)

        # Get average profiles and append to dict
        ave = {
            "signal" : sample,
            "average" : df.apply(np.mean, axis=0),                              # both strands
            "positive" : df.ix[range(0, len(df), 2)].apply(np.mean, axis=0),    # positive strand
            "negative" : df.ix[range(1, len(df), 2)].apply(np.mean, axis=0)     # negative strand
        }
        aveSignals[sample] = pd.DataFrame(ave)

# serialize with pickle
with open(os.path.join(coverageFilePath, "tssSignals_{0}bp.rawSignals.pickle".format(str(windowWidth))), 'wb') as handle:
    pickle.dump(rawSignals, handle, protocol=pickle.HIGHEST_PROTOCOL)
# with open(os.path.join(coverageFilePath, "tssSignals_{0}bp.rawSignals.pickle".format(str(windowWidth))), 'rb') as handle:
#     rawSignals = pickle.load(handle)

with open(os.path.join(coverageFilePath, "tssSignals_{0}bp.aveSignals.pickle".format(str(windowWidth))), 'wb') as handle:
    pickle.dump(aveSignals, handle, protocol=pickle.HIGHEST_PROTOCOL)
# with open(os.path.join(coverageFilePath, "tssSignals_{0}bp.aveSignals.pickle".format(str(windowWidth))), 'rb') as handle:
#     aveSignals = pickle.load(handle)

# Plot average profiles for all signals with ggplot
aveSignals = pd.concat(aveSignals, axis=0)
aveSignals["x"] = [i[1] for i in aveSignals.index]
aveSignals = pd.melt(aveSignals, id_vars=["x", "signal"])

plotFunc = robj.r("""
    library(ggplot2)

    #cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
   
    function(df, plotsDir){
        p = ggplot(df, aes(x, value, colour = variable)) +
            geom_line() +
            #stat_smooth(method = "gam", formula = y ~ s(x, k = 20), se = FALSE) + 
            facet_wrap(~signal, scales="free", ncol=1) +
            xlab("Distance to peak") +
            ylab("Average tags per bp") +
            theme_bw() +
            theme(legend.title=element_blank()) +
            #scale_colour_manual(values=cbPalette[c("grey", ,3,7)])
            scale_colour_manual(values=c("black", "dark blue","dark red")) +
            xlim(c(-200, 200))

        
        ggsave(filename = paste0(plotsDir, "/tss_signal.pdf"), plot = p, height = 7, width = 7)
        ggsave(filename = paste0(plotsDir, "/tss_signal.wide.pdf"), plot = p, height = 7, width = 10)
    }
""")

gr = importr('grDevices')
# convert the pandas dataframe to an R dataframe
robj.pandas2ri.activate()
aveSignals_R = robj.conversion.py2ri(aveSignals)
 
# run the plot function on the dataframe
plotFunc(aveSignals_R, plotsDir)

# Loop through raw signals, normalize, k-means cluster, save pickle, plot heatmap, export cdt
for signal in signals:
    exportName = "{0}.tssSignal_{1}bp.kmeans_{2}k".format(signal, str(windowWidth), n_clusters)
    
    df = rawSignals[signal]

    # join strand signal (plus - minus)
    df = df.xs('+',level="strand") - df.xs('-',level="strand")

    # scale row signal to 0:1 (normalization)
    dfNorm = df.apply(lambda x : (x - min(x)) / (max(x) - min(x)), axis=1)

    clust = k_means(dfNorm,
        n_clusters,
        n_init=25,
        max_iter=10000,
        n_jobs=2
    ) # returns centroid, label, inertia
    
    # save object
    pickle.dump(clust,
        open(os.path.join(bamFilePath, exportName + ".pickl"), "wb"),
        protocol = pickle.HIGHEST_PROTOCOL
    )

    # Sort dataframe by cluster order 
    dfNorm["cluster"] = clust[1] # clust[1] <- label from k-means clustering
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

loopsTask = Coverage(peaks_in_loops, 30, os.path.join(bamFilePath, sample + ".bam"),
    orientation=True, fragment_size=1, strand_specific=False, queue="shortq", cpusPerTask=2
)
notLoopsTask = Coverage(peaks_notin_loops, 30, os.path.join(bamFilePath, sample + ".bam"),
    orientation=True, fragment_size=1, strand_specific=False, queue="shortq", cpusPerTask=2
)
slurm.add_task(loopsTask)
slurm.add_task(notLoopsTask)

slurm.submit(loopsTask)
slurm.submit(notLoopsTask)

loopsCov = loopsTask.collect()
notLoopsCov = notLoopsTask.collect()

levels = [loopsCov.keys(), ["+", "-"]]
labels = [[y for x in range(len(loopsCov)) for y in [x, x]], [y for x in range(len(loopsCov.keys())) for y in (0,1)]]
index = pd.MultiIndex(labels=labels, levels=levels, names=["peak", "strand"])
loopsDf = pd.DataFrame(np.vstack(loopsCov.values()), index=index)
loopsDf.columns = range(windowRange[0], windowRange[1])

levels = [notLoopsCov.keys(), ["+", "-"]]
labels = [[y for x in range(len(notLoopsCov)) for y in [x, x]], [y for x in range(len(notLoopsCov.keys())) for y in (0,1)]]
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
loopsAve = {
    "signal" : sample,
    "average" : loopsDf.apply(np.mean, axis=0),                              # both strands
    "positive" : loopsDf.ix[range(0, len(loopsDf), 2)].apply(np.mean, axis=0),    # positive strand
    "negative" : loopsDf.ix[range(1, len(loopsDf), 2)].apply(np.mean, axis=0)     # negative strand
}
notLoopsAve = {
    "signal" : sample,
    "average" : notLoopsDf.apply(np.mean, axis=0),                              # both strands
    "positive" : notLoopsDf.ix[range(0, len(notLoopsDf), 2)].apply(np.mean, axis=0),    # positive strand
    "negative" : notLoopsDf.ix[range(1, len(notLoopsDf), 2)].apply(np.mean, axis=0)     # negative strand
}
pickle.dump(loopsAve, open(os.path.join(coverageFilePath, "{0}_inLoops_average.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(notLoopsAve, open(os.path.join(coverageFilePath, "{0}_notinLoops_average.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
# loopsAve = pickle.load(open(os.path.join(coverageFilePath, "{0}_inLoops_average.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "r"))
# notLoopsAve = pickle.load(open(os.path.join(coverageFilePath, "{0}_notinLoops_average.tssSignal_{1}bp.pickle".format(sample, str(windowWidth))), "r"))


# Plot
import matplotlib.pyplot as plt
plt.plot(loopsAve['positive'].index, loopsAve['positive'])
plt.plot(loopsAve['negative'].index, loopsAve['negative'])
plt.plot(loopsAve['average'].index, loopsAve['average'])
plt.plot(notLoopsAve['positive'].index, notLoopsAve['positive'])
plt.plot(notLoopsAve['negative'].index, notLoopsAve['negative'])
plt.plot(notLoopsAve['average'].index, notLoopsAve['average'])
plt.show()
