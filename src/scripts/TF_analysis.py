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

from ggplot import ggplot, aes, stat_smooth, facet_wrap, xlab, ylab, theme_bw, theme

import rpy2.robjects as robj # for ggplot in R
import rpy2.robjects.pandas2ri # for R dataframe conversion
from rpy2.robjects.packages import importr

from sklearn.cluster import k_means 
import cPickle as pickle

import plotly.plotly as plotly
from plotly.graph_objs import Data, Heatmap, Layout, Figure

global v # verbose
# Define variables
v = True
signals = { "CTCF_K562_10mio_CM" : ["CTCF_K562_10mio_CM", "CTCF_K562_10mio_ChIP"],
            "PU1_K562_10mio_CM" : ["PU1_K562_10mio_CM", "PU1_K562_10mio_ChIP"]}
bamFilePath = "/home/arendeiro/data/human/chipmentation/mapped/merged"
peakFilePath = "/home/arendeiro/data/human/chipmentation/bed"
coverageFilePath = "/home/arendeiro/data/human/chipmentation/bed/newPython" # output
plotsDir = "/home/arendeiro/"

genome = "hg19"
windowRange = (-2000, 2000)
fragmentsize = 1
duplicates = True
n_clusters = 5

plotly.sign_in("afrendeiro", "iixmygxac1")
windowWidth = abs(windowRange[0]) + abs(windowRange[1])


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


def coverage(bam, intervals, fragmentsize, orientation=True, duplicates=True, strand_specific=False):
    """
    Gets read coverage in bed regions.
    Returns dict of regionName:numpy.array if strand_specific=False, A dict of "+" and "-" keys with regionName:numpy.array
    bam - HTSeq.BAM_Reader object. Must be sorted and indexed with .bai file!
    intervals - dict with HTSeq.GenomicInterval objects as values
    fragmentsize - integer
    stranded - boolean
    duplicates - boolean.
    """
    # Loop through TSSs, get coverage, append to dict
    cov = OrderedDict()
    i=0
    for name, feature in intervals.iteritems():
        if i%1000 == 0:
            print(i)
        # Initialize empty array for this feature
        if not strand_specific:
            profile = np.zeros(feature.length, dtype=np.int8)
        else:
            profile = np.zeros((2, feature.length), dtype=np.int8)
                    
        # Fetch alignments in feature window
        for aln in bam[feature]:
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
            if start_in_window < 0 or end_in_window > feature.length:
                continue

            # add +1 to all positions overlapped by read within window
            if not strand_specific:
                profile[start_in_window : end_in_window] += 1
            else:
                if aln.iv.strand == "+":
                    profile[0][start_in_window : end_in_window] += 1
                elif aln.iv.strand == "-":
                    profile[1][start_in_window : end_in_window] += 1
        
        # append feature profile to dict
        cov[name] = profile
        i+=1
    return cov


def exportToJavaTreeView(df, filename):
    """
    Export cdt file of cluster to view in JavaTreeView
    df - pandas.DataFrame object
    filename - string.    
    """
    df["X"] = df.index
    df["NAME"] = df.index
    df["GWEIGHT"] = 1
    df.columns = ["X" + str(x) for x in df.columns]
    df = df[["X", "NAME", "GWEIGHT"] + ["X" + str(x) for x in df.columns]]
    df.to_csv(filename, sep="\t", index=False)
    

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

    # Subset
    # i=0
    # d = dict()
    # for name, peak in peaks.iteritems():
    #     if i < 500:
    #         d[name] = peak
    #         i+=1
    # peaks = d

    for sample in signals[signal]:
        # Load bam
        bamfile = HTSeq.BAM_Reader(os.path.join(bamFilePath, sample + ".bam"))

        # Get dataframe of signal coverage in bed regions, append to dict
        cov = coverage(bamfile, peaks, fragmentsize, strand_specific=True)

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
with open(os.path.join(coverageFilePath, "tssSignals_{0}bp.rawSignals.pickl".format(str(windowWidth))), 'wb') as handle:
    pickle.dump(rawSignals, handle, protocol=pickle.HIGHEST_PROTOCOL)
# with open(os.path.join(coverageFilePath, "tssSignals_{0}bp.rawSignals.pickl".format(str(windowWidth))), 'rb') as handle:
#     rawSignals = pickle.load(handle)

with open(os.path.join(coverageFilePath, "tssSignals_{0}bp.aveSignals.pickl".format(str(windowWidth))), 'wb') as handle:
    pickle.dump(aveSignals, handle, protocol=pickle.HIGHEST_PROTOCOL)
# with open(os.path.join(coverageFilePath, "tssSignals_{0}bp.aveSignals.pickl".format(str(windowWidth))), 'rb') as handle:
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
