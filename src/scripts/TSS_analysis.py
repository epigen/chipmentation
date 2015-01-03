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
import HTSeq
import numpy as np
import pandas as pd

import rpy2.robjects as robj # for ggplot in R
import rpy2.robjects.pandas2ri # for R dataframe conversion

from sklearn.cluster import k_means 
import cPickle as pickle

import plotly.plotly as plotly
from plotly.graph_objs import Data, Heatmap, Layout, Figure

# Define variables
signals = ["H3K4me3_K562_500k_CM", "H3K4me3_K562_500k_ChIP", "IgG_K562_500k_CM", "DNase..."]
bamFilePath = "/home/arendeiro/data/human/chipmentation/mapped/merged"
bedFilePath = "/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.TATA_Annotated.bed"
plotsDir = "/home/arendeiro/projects/chipmentation/results/plots"

windowRange = (-60, 60)
fragmentsize = 1
duplicates = True
n_clusters = 5

plotly.sign_in("afrendeiro", "iixmygxac1")
windowWidth = abs(windowRange[0]) + abs(windowRange[1])


def loadBed(filename, windowWidth=0):
    """
    Parses bed file and returns dict of featureName:HTSeq.GenomicInterval objects
    filename - string.
    """
    import HTSeq.GenomicInterval

    features = dict()
    for line in open(filename):
        fields = line.split("\t")
        feature = HTSeq.GenomicInterval(
            fields[0],                      # chrom
            int(fields[1]) - windowWidth/2, # start
            int(fields[2]) + windowWidth/2, # end
            fields[5]                       # strand
        )
        features[fields[4]] = feature.start_as_pos # append with name
    return features


def coverage(bam, intervals, fragmentsize, stranded=True, duplicates=True):
    """
    Gets read coverage in bed regions, returns dict of regionName:numpy.array
    bam - HTSeq.BAM_Reader object. Must be sorted and indexed with .bai file!
    intervals - dict with HTSeq.GenomicInterval objects as values
    fragmentsize - integer
    stranded - boolean
    duplicates - boolean.
    """
    # Loop through TSSs, get coverage, append to dict
    cov = dict()

    for name, feature in intervals.iteritems():
        # Initialize empty array for this feature
        profile = np.zeros(feature.length, dtype=np.int8)
        
        # Fetch alignments in feature window
        for aln in bam[window]:
            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue

            aln.iv.length = fragmentsize # adjust to size
            
            # get position in relative to window
            if stranded:
                if feature.strand == "+":
                    start_in_window = aln.iv.start - feature.start
                    end_in_window   = aln.iv.end   - feature.start
                else:
                    start_in_window = feature.start - aln.iv.end
                    end_in_window   = feature.start - aln.iv.start
            else:
                start_in_window = aln.iv.start - feature.start 
                end_in_window   = aln.iv.end   - feature.start

            # add +1 to all positions overlapped by read within window
            profile[start_in_window : end_in_window] += 1
        
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
tsss = loadBed(bedFilePath, windowWidth)

# Loop through all signals, compute coverage in bed regions,
# save dicts with coverage and average profiles
rawSignals = dict()
aveSignals = dict()
for signal in signals:
    # Load bam
    bamfile = HTSeq.BAM_Reader(os.path.join(bamFilePath, signal + ".bam"))

    # Get dataframe of signal coverage in bed regions, append to dict
    cov = coverage(bamfile, tsss, windowWidth, fragmentsize)
    df = pd.DataFrame(cov).T
    df.columns = range(windowRange[0], windowRange[1])
    ### TODO: check that df index is tss name!!    
    rawSignals[signal] = df

    # Save as csv
    df.to_csv(os.path.join(bamFilePath, "{1}.tssSignal_{2}bp.csv".format(signal, str(windowWidth))))

    # Get average profile and append to dict
    ave = df.apply(np.mean, axis=0)
    aveSignals[signal] = ave


# Plot average profiles for all signals with ggplot2
aveSignals = pd.DataFrame(aveSignals)
aveSignals["x"] = aveSignals.index
aveSignals = pd.melt(aveSignals, id_vars="x")

plotFunc = robj.r("""
    library(ggplot2)

    function(df, plotsDir){
        p = ggplot(df, aes(x, value, colour = variable)) +
            stat_smooth(method = "gam", formula = y ~ s(x, k = 20), se = FALSE) + 
            facet_wrap( ~ variable) +
            xlab("Distance to peak") +
            ylab("Average tags per bp") +
            theme_bw() +
            theme(legend.title=element_blank())
        
        ggsave(filename = paste0(plotsDir, "/tss_signal.pdf"), plot = p, height = 2, width = 7)
        ggsave(filename = paste0(plotsDir, "/tss_signal.wide.pdf"), plot = p, height = 2, width = 4)
    }
""")

# convert the pandas dataframe to an R dataframe
robj.pandas2ri.activate()
aveSignal_R = robj.conversion.py2ri(aveSignal)
 
# run the plot function on the dataframe
plotFunc(aveSignal_R, plotsDir)


# Loop through raw signals, normalize, k-means cluster, save pickle, plot heatmap, export cdt
for signal in signals:
    exportName = "{0}.tssSignal_{1}bp.kmeans_{2}k".format(signal, str(windowWidth), n_cluster)
    
    df = rawSignals[signal]

    # scale row signal to 0:1 (normalization)
    dfNorm = df.apply(lambda x : (x - min(x)) / (max(x) - min(x)), axis=1)

    clust = k_means(dfNorm,
        n_clusters,
        n_init=25,
        max_iter=10000,
        n_jobs=-1
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
