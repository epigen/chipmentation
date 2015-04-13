#!/usr/env python

#############################################################################################
#
# This code was used to produce the plots in the ChIPmentation paper (Schmidl, et al. 2015).
#
# It produces plots of sample correlation over windows genome-wide.
#
#############################################################################################

import os
import re
import numpy as np
import pandas as pd
import itertools
import rpy2.robjects as robj  # for ggplot in R
import rpy2.robjects.pandas2ri  # for R dataframe conversion

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def getCounts(sampleSubset):
    normCounts = pd.DataFrame()

    for sample in range(len(sampleSubset)):
        name = sampleSubset.ix[sample]["sampleName"]

        # read in bed file with coverage
        try:
            df = pd.read_csv(os.path.join(coverageDir, name + ".cov"), sep="\t", header=None)
        except:
            print("Couldn't open file %s." % os.path.join(coverageDir, name + ".cov"))
            continue

        # get series with counts
        counts = df[3]

        # normalize by total size
        norm = np.log2(1 + (counts / sum(counts)) * 1000000)

        # append to raw counts
        name = re.sub("_hg19", "", name)
        name = re.sub("K562_", "", name)
        normCounts[name] = norm
    return normCounts


def colourPerFactor(name):
    name = str(name.upper())
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
    elif "NAN_NAN_NAN" in name:  # DNAse, ATAC
        return "#00523b"
    else:
        raise ValueError


def plotCorrelations(normCounts, plotName):
    colours = map(colourPerFactor, normCounts.columns)
    cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)

    sns.clustermap(normCounts.corr(), row_colors=colours, method="ward",
                   col_colors=colours, figsize=(15, 15), cmap=cmap)

    plt.savefig(plotName, bbox_inches='tight')


def plotScatterCorrelation(df, path):
    x = df.columns[0]
    y = df.columns[1]

    # Plot with R
    plotFunc = robj.r("""
    library(ggplot2)

    function(df, path, x, y) {

        pdf(path, height = 12, width = 12)

        # scatterplot
        smoothScatter(
            log2(df$A), log2(df$B),
            col = rgb(104,104,104,50, maxColorValue = 255),
            pch = 16,
            xlab = x,
            ylab = y,
            nrpoints = 0,
            ann = FALSE#,
            #xaxt = 'n'
        )
        # pearson's R as text
        text(
            par('usr')[1] + 1.8,
            par('usr')[4] - 0.5,
            bquote(R^2 == .(round(cor(df$A, df$B), 3))),
            cex = 1.6
        )
        dev.off()
    }
    """)
    df.columns = ["A", "B"]

    # convert the pandas dataframe to an R dataframe
    robj.pandas2ri.activate()
    df_R = robj.conversion.py2ri(df)

    # run the plot function on the dataframe
    plotFunc(df_R, path, x, y)


def plotKDE(s1, s2, path):
    sns.jointplot(s1, s2, kind="kde", size=7, space=0)
    plt.savefig(path, bbox_inches='tight')


# Define paths
projectRoot = "/projects/chipmentation/"
coverageDir = projectRoot + "data/coverage"
resultsDir = projectRoot + "results"
plotsDir = resultsDir + "/plots/correlations"

# Get samples
samples = pd.read_csv(os.path.abspath(projectRoot + "chipmentation.replicates.annotation_sheet.csv"))

# set technicalreplicate=0 to samples with only one technical replicate per biological replicate
for n, i in samples.groupby(["cellLine", "numberCells", "technique", "ip",
                             "treatment", "biologicalReplicate", "genome"]).groups.items():
    if len(i) == 1:
        samples.loc[i[0], "technicalReplicate"] = 0

# All
# subset samples
sampleSubset = samples.reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

# normCounts = getCounts(sampleSubset)
# plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.all.pdf"))


# All biological replicates
print("All biologicalReplicates")
# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("H3K4ME3|H3K4ME1|H3K27ME3|H3K36ME3|H3K27AC|CTCF|PU1|GATA1|REST")) &
    (samples["numberCells"].str.contains("10M|500K|100K|10K")) &
    (samples["technicalReplicate"] == 0) &
    (samples["biologicalReplicate"] != 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

normCounts = getCounts(sampleSubset)
normCounts.to_pickle(os.path.join(plotsDir, "correlations.set-reps.pickle"))
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.set-reps.pdf"))

# All final
print("All final")
# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("H3K4ME3|H3K4ME1|H3K27ME3|H3K36ME3|H3K27AC|CTCF|PU1|GATA1|REST")) &
    (samples["numberCells"].str.contains("10M|500K|100K|10K")) &
    (samples["technicalReplicate"] == 0) &
    (samples["biologicalReplicate"] == 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

normCounts = getCounts(sampleSubset)
normCounts.to_pickle(os.path.join(plotsDir, "correlations.set.pickle"))
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.set.pdf"))

# Histones
print("histonesAll")
# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("H3K4ME3|H3K4ME1|H3K27ME3|H3K36ME3|H3K27AC")) &
    (samples["numberCells"].str.contains("10M|500K|10K")) &
    (samples["technicalReplicate"] == 0) &
    (samples["biologicalReplicate"] != 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

normCounts = getCounts(sampleSubset)
normCounts.to_pickle(os.path.join(plotsDir, "correlations.histones-reps.pickle"))
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.histones-reps.pdf"))

# Histones
print("histones biologicalReplicate=0")
# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("H3K4ME3|H3K4ME1|H3K27ME3|H3K36ME3|H3K27AC")) &
    (samples["numberCells"].str.contains("10M|500K|10K")) &
    (samples["technicalReplicate"] == 0) &
    (samples["biologicalReplicate"] == 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

normCounts = getCounts(sampleSubset)
normCounts.to_pickle(os.path.join(plotsDir, "correlations.histones.pickle"))
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.histones.pdf"))

# TFs
print("TFs")
# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("CTCF|PU1|GATA1|REST")) &
    (samples["numberCells"].str.contains("10M|500K|100K|10K")) &
    (samples["technicalReplicate"] == 0) &
    (samples["biologicalReplicate"] != 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

normCounts = getCounts(sampleSubset)
normCounts.to_pickle(os.path.join(plotsDir, "correlations.TFs-reps.pickle"))
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.TFs-reps.pdf"))

# TFs
print("TFs")
# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("CTCF|PU1|GATA1|REST")) &
    (samples["numberCells"].str.contains("10M|500K|100K|10K")) &
    (samples["technicalReplicate"] == 0) &
    (samples["biologicalReplicate"] == 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

normCounts = getCounts(sampleSubset)
normCounts.to_pickle(os.path.join(plotsDir, "correlations.TFs.pickle"))
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.TFs.pdf"))

# Figure 1:
#   correlations between biological == 0 for each factor
#   in combinations of different techniques, number of cells and replicates

print("techniques")
# Different techniques
for ip in samples['ip'].unique():
    for c in samples['numberCells'].unique():
        cells = samples['technique'].unique()
        for (t1, t2) in itertools.combinations(cells, 2):
            s1 = samples[
                (samples["ip"] == ip) &
                (samples["numberCells"] == c) &
                (samples["technique"] == t1) &
                (samples["technicalReplicate"] == 0) &
                (samples["biologicalReplicate"] == 0)
            ]
            s2 = samples[
                (samples["ip"] == ip) &
                (samples["numberCells"] == c) &
                (samples["technique"] == t2) &
                (samples["technicalReplicate"] == 0) &
                (samples["biologicalReplicate"] == 0)
            ]
            if s1.empty or s2.empty:
                continue

            normCounts = getCounts(s1.append(s2).reset_index())
            normCounts.columns = s1.sampleName.tolist() + s2.sampleName.tolist()

            plotScatterCorrelation(
                normCounts,
                os.path.join(
                    plotsDir,
                    "correlation.{0}.{1}.{2}_vs_{3}.pdf".format(ip, c, t1, t2)
                )
            )

print("cells")
# Different number of cells
for ip in samples['ip'].unique():
    for t in samples['technique'].unique():
        cells = samples['numberCells'].unique()
        for (c1, c2) in itertools.combinations(cells, 2):
            s1 = samples[
                (samples["ip"] == ip) &
                (samples["numberCells"] == c1) &
                (samples["technique"] == t) &
                (samples["technicalReplicate"] == 0) &
                (samples["biologicalReplicate"] == 0)
            ]
            s2 = samples[
                (samples["ip"] == ip) &
                (samples["numberCells"] == c2) &
                (samples["technique"] == t) &
                (samples["technicalReplicate"] == 0) &
                (samples["biologicalReplicate"] == 0)
            ]
            if s1.empty or s2.empty:
                continue

            normCounts = getCounts(s1.append(s2).reset_index())
            normCounts.columns = s1.sampleName.tolist() + s2.sampleName.tolist()

            plotScatterCorrelation(
                normCounts,
                os.path.join(
                    plotsDir,
                    "correlation.{0}.{1}.{2}_vs_{3}.pdf".format(ip, t, c1, c2)
                )
            )

print("replicates")
# Replicate 1 vs 2
for ip in samples['ip'].unique():
    for c in samples['numberCells'].unique():
        for t in samples['technique'].unique():
            r1 = 1
            r2 = 2
            s1 = samples[
                (samples["ip"] == ip) &
                (samples["numberCells"] == c) &
                (samples["technique"] == t) &
                (samples["technicalReplicate"] == 0) &
                (samples["biologicalReplicate"] == r1)
            ]
            s2 = samples[
                (samples["ip"] == ip) &
                (samples["numberCells"] == c) &
                (samples["technique"] == t) &
                (samples["technicalReplicate"] == 0) &
                (samples["biologicalReplicate"] == r2)
            ]
            if s1.empty or s2.empty:
                continue

            normCounts = getCounts(s1.append(s2).reset_index())
            normCounts.columns = s1.sampleName.tolist() + s2.sampleName.tolist()

            plotScatterCorrelation(
                normCounts,
                os.path.join(
                    plotsDir,
                    "correlation.{0}.{1}.{2}.{3}_vs_{4}.pdf".format(ip, c, t, r1, r2)
                )
            )
