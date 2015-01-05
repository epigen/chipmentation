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

# Define variables
signals = ["H3K4me3_K562_500k_CM", "H3K4me3_K562_500k_ChIP", "IgG_K562_500k_CM", "DNase..."]
bamFilePath = "/home/arendeiro/data/human/chipmentation/mapped/merged"
bedFilePath = "/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.TATA_Annotated.bed"
plotsDir = "/home/arendeiro/projects/chipmentation/results/plots"

bamFilePath = "/home/afr/"
bedFilePath = "/home/afr/hg19.cage_peak_coord_robust.TATA_Annotated.bed"
signals = ["H3K4me2", "H3K4me3"]
plotsDir = "/home/afr/"

genome = "hg19"
windowRange = (-60, 60)
fragmentsize = 1
duplicates = True
n_clusters = 5

plotly.sign_in("afrendeiro", "iixmygxac1")
windowWidth = abs(windowRange[0]) + abs(windowRange[1])


def loadBed(filename):
    """
    Parses bed file and returns dict of featureName:HTSeq.GenomicInterval objects
    filename - string.
    """
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
    """
    intervals = OrderedDict()
    for iv in bedtool:
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

    for name, feature in intervals.iteritems():
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
                    start_in_window = aln.iv.start - feature.start
                    end_in_window   = aln.iv.end   - feature.start
                else:
                    start_in_window = feature.length - abs(feature.start - aln.iv.end)
                    end_in_window   = feature.length - abs(feature.start - aln.iv.start)
            else:
                start_in_window = aln.iv.start - feature.start 
                end_in_window   = aln.iv.end   - feature.start
            
            # check fragment is within window; this is because of fragmentsize adjustment
            if start_in_window <= 0 or end_in_window > feature.length:
                continue

            # add +1 to all positions overlapped by read within window
            if not strand_specific:
                profile[start_in_window : end_in_window] += 1
            else:
                if aln.iv.strand == "+":
                    profile[0][start_in_window : end_in_window] += 1
                else:
                    profile[1][start_in_window : end_in_window] += 1
        
        # append feature profile to dict
        cov[name] = profile
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
    

# Load TSSs from bed file, transform by window width
tsss = pybedtools.BedTool(bedFilePath).slop(genome=genome, b=windowWidth/2)
tsss = bedToolsInterval2GenomicInterval(tsss)
# Filter tsss near chrm borders
for name, interval in tsss.iteritems():
    if interval.length < windowWidth + 1:
        tsss.pop(name)

# Loop through all signals, compute coverage in bed regions,
# save dicts with coverage and average profiles
rawSignals = OrderedDict()
aveSignals = OrderedDict()
for signal in signals:
    # Load bam
    bamfile = HTSeq.BAM_Reader(os.path.join(bamFilePath, signal + ".bam"))

    # Get dataframe of signal coverage in bed regions, append to dict
    cov = coverage(bamfile, tsss, fragmentsize, strand_specific=True)

    # Make multiindex dataframe
    levels = [cov.keys(), ["+", "-"]] 
    labels = [[y for x in range(len(cov)) for y in [x, x]], [y for x in range(len(cov.keys())) for y in (0,1)]]
    index = pd.MultiIndex(labels=labels, levels=levels, names=["tss", "strand"])
    df = pd.DataFrame(np.vstack(cov.values()), index=index)
    df.columns = range(windowRange[0], windowRange[1] + 1)

    # For strand_specific=False
    #cov = coverage(bamfile, tsss, fragmentsize)
    #df = pd.DataFrame(cov).T
    #df.columns = range(windowRange[0], windowRange[1])
    
    # append to dict   
    rawSignals[signal] = df

    # Save as csv
    #df.to_csv(os.path.join(bamFilePath, "{1}.tssSignal_{2}bp.csv".format(signal, str(windowWidth))))

    # Get average profiles and append to dict
    ave = {
        "signal" : signal,
        "average" : df.apply(np.mean, axis=0),                              # both strands
        "positive" : df.ix[range(0, len(df), 2)].apply(np.mean, axis=0),    # positive strand
        "negative" : df.ix[range(1, len(df), 2)].apply(np.mean, axis=0)     # negative strand
    }
    aveSignals[signal] = pd.DataFrame(ave)


# Plot average profiles for all signals with ggplot
aveSignals = pd.concat(aveSignals, axis=0)
aveSignals["x"] = [i[1] for i in aveSignals.index]
aveSignals = pd.melt(aveSignals, id_vars=["x", "signal"])

plotFunc = robj.r("""
    library(ggplot2)

    #cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
   
    function(df, plotsDir){
        p = ggplot(df, aes(x, value, colour = variable)) +
            stat_smooth(method = "gam", formula = y ~ s(x, k = 20), se = FALSE) + 
            facet_wrap(~signal) +
            xlab("Distance to peak") +
            ylab("Average tags per bp") +
            theme_bw() +
            theme(legend.title=element_blank()) +
            #scale_colour_manual(values=cbPalette[c("grey", ,3,7)])
            scale_colour_manual(values=c("grey", "dark blue","dark red"))
        
        ggsave(filename = paste0(plotsDir, "/tss_signal.pdf"), plot = p, height = 2, width = 7)
        ggsave(filename = paste0(plotsDir, "/tss_signal.wide.pdf"), plot = p, height = 2, width = 4)
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
