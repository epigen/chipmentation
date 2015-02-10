#!/usr/env python

#############################################################################################
# 
# This code was used to produce the plots in the ChIPmentation paper (Schmidl, et al. 2015).
# 
# It produces plots of sample correlation over windows genome-wide.
#
#############################################################################################

import os
from pybedtools import BedTool
import HTSeq
import string
import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr

import rpy2.robjects as robj # for ggplot in R
import rpy2.robjects.pandas2ri # for R dataframe conversion

# Define variables
signals = [
    "H3K4me3_K562_500k_CM", "H3K4me3_K562_500k_ChIP", "H3K4me3_K562_10k_CM",
    "H3K27me3_K562_500k_CM", "H3K27me3_K562_500k_ChIP", "H3K27me3_K562_10k_CM",
]
bamFilePath = "/home/arendeiro/data/human/chipmentation/mapped/merged/"
genomeSizes = "/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt"
plotsDir = "home/arendeiro/projects/chipmentation/results/plots"

genome = "hg19"
windowWidth = 1000
fragmentsize = 1
duplicates = False

def coverage(bam, intervals, fragmentsize, duplicates=False):
    """ Gets read coverage in bed regions, returns dict with region:count.
    bam - Bam object from HTSeq.BAM_Reader.
    intervals - dict with HTSeq.GenomicInterval objects as values.
    fragmentsize - integer.
    duplicates - boolean.
    """
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    # Loop through TSSs, get coverage, append to dict
    cov = dict()

    for name, feature in intervals.iteritems():
        if feature.chrom not in chroms:
            continue
        count = 0
        
        # Fetch alignments in feature window
        for aln in bam[feature]:
            # check if duplicate
            if not duplicates and aln.pcr_or_optical_duplicate:
                continue

            # adjust fragment to size
            aln.iv.length = fragmentsize
            
            # add +1 to all positions overlapped by read within window
            count += 1
        
        # append feature profile to dict
        cov[name] = count
    return cov


# Generate 1kb windows genome-wide
w = BedTool.window_maker(BedTool(), genome=genomeSizes, w=windowWidth)
windows = dict()
for interval in w:
    feature = HTSeq.GenomicInterval(
        interval.chrom,
        interval.start,
        interval.end
    )
    name = string.join(interval.fields, sep="_")
    windows[name] = feature

# Loop through all signals, compute coverage in bed regions, append to dict
rawSignals = dict()
for signal in signals:
    # Load bam
    bamfile = HTSeq.BAM_Reader(os.path.join(bamFilePath, signal + ".bam"))

    # Get dataframe of signal coverage in bed regions, append to dict
    rawSignals[signal] = coverage(bamfile, windows, fragmentsize, duplicates)

df = pd.DataFrame(rawSignals)

# Normalize to library size 
dfNorm = df.apply(lambda x: np.log2( 1 + (x / x.sum()) * 1000000))

# Plot with R
plotFunc = robj.r("""
    library(ggplot2)
    library(reshape2)

    function(df, plotsDir){
        # scatterplot 
        pdf(paste(projectDir, "results/plots/correlations_1kb_windows_scatter_mergedReplicates.pdf", sep = ""), height = 11, width = 7)
        m = matrix(c(1,2,3,4,5,5,6,7) ,nrow = 4,ncol = 2, byrow = TRUE)
        layout(mat = m, heights = c(1,1,0.2,1))

        par(oma = c(5,5,2,0) + 0.1,
            mar = c(0,0,0,0) + 0.1,
            pty = "s")

        smoothScatter(log2(rawCounts$H3K4me3_K562_500k_CM), log2(rawCounts$H3K4me3_K562_500k_ChIP),
            col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n')
        text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
            bquote(R^2 == .(round(cor(rawCounts$H3K4me3_K562_500k_CM, rawCounts$H3K4me3_K562_500k_ChIP), 3))),
            cex = 1.6
        )
        mtext(side = 2, "500.000 cells", line = 2)
        mtext(side = 3, "H3K4me3")

        smoothScatter(log2(rawCounts$H3K27me3_K562_500k_CM), log2(rawCounts$H3K27me3_K562_500k_ChIP),
            col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n', yaxt = 'n')
        text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
            bquote(R^2 == .(round(cor(rawCounts$H3K27me3_K562_500k_CM, rawCounts$H3K27me3_K562_500k_ChIP), 3))),
            cex = 1.6
        )
        mtext(side = 3, "H3K27me3")

        smoothScatter(log2(rawCounts$H3K4me3_K562_10k_CM), log2(rawCounts$H3K4me3_K562_500k_ChIP),
            col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
        mtext(side = 2, "10.000 cells", line = 2)
        text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
            bquote(R^2 == .(round(cor(rawCounts$H3K4me3_K562_10k_CM, rawCounts$H3K4me3_K562_500k_ChIP), 3))),
            cex = 1.6
        )

        smoothScatter(log2(rawCounts$H3K27me3_K562_10k_CM), log2(rawCounts$H3K27me3_K562_500k_ChIP),
            col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, yaxt = 'n')
        text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
            bquote(R^2 == .(round(cor(rawCounts$H3K27me3_K562_10k_CM, rawCounts$H3K27me3_K562_500k_ChIP), 3))),
            cex = 1.6
        )
        mtext(side = 1, "ChIP", outer= FALSE, at = -6.2, cex = 1, line = 3)

        plot.new()
        #mtext(side = 1, "ChIP", line = 0, at = 0)

        smoothScatter(log2(rawCounts$H3K4me3_K562_500k_CM), log2(rawCounts$H3K4me3_K562_10k_CM),
            col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
        mtext(side = 2, "500.000 cells", line = 2)
        text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
            bquote(R^2 == .(round(cor(rawCounts$H3K4me3_K562_500k_CM, rawCounts$H3K4me3_K562_10k_CM), 3))),
            cex = 1.6
        )
        mtext(side = 1, "10.000 cells", line = 2, cex.lab = 1.5)

        smoothScatter(log2(rawCounts$H3K27me3_K562_500k_CM), log2(rawCounts$H3K27me3_K562_10k_CM),
            col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, yaxt = 'n')
        text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
            bquote(R^2 == .(round(cor(rawCounts$H3K27me3_K562_500k_CM, rawCounts$H3K27me3_K562_10k_CM), 3))),
            cex = 1.6
        )
        mtext(side = 1, "10.000 cells", line = 2, cex.lab = 1.5)

        title(xlab = "ChIPmentation",
              ylab = "ChIPmentation",
              outer = TRUE, cex.lab = 1.5)
        dev.off()
    }
""")

# convert the pandas dataframe to an R dataframe
robj.pandas2ri.activate()
df_R = robj.conversion.py2ri(dfNorm)
 
# run the plot function on the dataframe
plotFunc(df_R, plotsDir)



# TODO: implement multipleCorrelations function
# def multipleCorrelations(dataframe):
#     """Calculates the Pearson correlation between all columns in a pandas dataframe.
#     """
#     raise NotImplemented
#
# cor = multipleCorrelations(df)



# plotFunc = robj.r("""
#     library(ggplot2)
#     library(reshape2)

#     function(df, plotsDir){
#         # scatterplot 
#         pdf(paste0(plotsDir, "/correlations_1kb_windows.pdf", sep = ""), height = 11, width = 7)
        
#         smoothScatter(log2(df$H3K4me3), log2(df$H3K4me2),
#             col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n')
#         text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
#             bquote(R^2 == .(round(cor(df$H3K4me3, df$H3K4me2), 3))),
#             cex = 1.6
#         )
#         mtext(side = 2, "500.000 cells", line = 2)
#         mtext(side = 3, "H3K4me3")

        
#         title(xlab = "ChIPmentation",
#               ylab = "ChIPmentation",
#               outer = TRUE, cex.lab = 1.5)
#         dev.off()
#     }
# """)
