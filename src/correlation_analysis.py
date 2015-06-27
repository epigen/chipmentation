#!/usr/env python

#############################################################################################
#
# This code was used to produce the plots in the ChIPmentation paper (Schmidl, et al. 2015).
#
# It produces plots of sample correlation over windows genome-wide.
#
#############################################################################################

import os
import yaml
import re
import numpy as np
import pandas as pd
import rpy2.robjects as robj  # for ggplot in R
import rpy2.robjects.pandas2ri  # for R dataframe conversion
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")


def getCounts(sampleSubset):
    normCounts = pd.DataFrame()

    for name in sampleSubset['sampleName']:
        # read in bed file with coverage
        try:
            df = pd.read_csv(os.path.join(coverageDir, name + ".filtered.cov"), sep="\t", header=None)
        except:
            print("Couldn't open file %s." % os.path.join(coverageDir, name + ".filtered.cov"))
            continue

        # get series with counts
        counts = df[3]

        # normalize by total size
        norm = (counts / sum(counts)) * 1000000

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
        return "#00523b"
        #raise ValueError


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    import matplotlib.colors as colors
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def plotCorrelations(normCounts, plotName, method="ward", metric="euclidean", values=False, color="redgreen"):
    corr = normCounts.corr()

    # col/row colours
    colours = map(colourPerFactor, corr.columns)

    # data colour map
    if color == "redgreen":
        cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
    else:
        cmap = plt.get_cmap("YlGn")

    # scale from 0 to 1
    new_cmap = truncate_colormap(cmap, 0, 1)

    if not values:
        sns.clustermap(corr, row_colors=colours, method=method, metric=metric,
                       col_colors=colours, figsize=(15, 15), cmap=new_cmap)
    else:
        sns.clustermap(corr, row_colors=colours, method=method, metric=metric,
                       col_colors=colours, figsize=(15, 15), cmap=new_cmap, annot=True)

    plt.savefig(plotName, bbox_inches='tight')
    plt.close()


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
        # x = y line
        abline(0, 1, lty=3)
        # pearson's R as text
        text(
            par('usr')[1] + 1.8,
            par('usr')[4] - 0.5,
            bquote(r == .(round(cor(df$A, df$B), 3))),
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



# Read configuration file
with open(os.path.join(os.path.expanduser("~"), ".pipelines_config.yaml"), 'r') as handle:
    config = yaml.load(handle)


# Define paths
# coverageDir = os.path.join(config["paths"]["parent"], config["projectname"], "data", "coverage")
# projectRoot = "/media/afr/cemm-backup1/chipmentation/"
projectRoot = "projects/chipmentation/"
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


# Fig 1e: ChIPmentation titration
# subset samples
sampleSubset = samples[
    (samples["technique"] == "CM") &
    (samples["ip"] == "H3K4ME3") &
    (samples["numberCells"] == "500K") &
    (samples["treatment"].notnull())
].reset_index(drop=True)

# get counts per window
normCounts = getCounts(sampleSubset)
# save matrix
c = normCounts.corr()
c['sample1'] = c.index
c = pd.melt(c, id_vars=['sample1'])
c.to_csv(os.path.join(plotsDir, "fig1e-correlations.titration.tsv"))
plotCorrelations(normCounts, os.path.join(plotsDir, "fig1e-correlations.titration.pdf"), method="single", values=True)


# Fig 1g: ChIPmentation histone sample correlation
# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("H3K4ME3|H3K4ME1|H3K27ME3|H3K36ME3|H3K27AC")) &
    (samples["numberCells"].str.contains("10M|500K|10K")) &
    (~samples["treatment"].notnull()) &
    (samples["technicalReplicate"] == 0) &
    (samples["biologicalReplicate"] != 0) &
    (samples["sampleName"] != "K562_500K_CHIP_H3K27AC_nan_nan_1_1_hg19")
].reset_index(drop=True)

# get counts per window
normCounts = getCounts(sampleSubset)
# save matrix
c = normCounts.corr()
c['sample1'] = c.index
c = pd.melt(c, id_vars=['sample1'])
c.to_csv(os.path.join(plotsDir, "fig1e-correlations.histones.tsv"))
plotCorrelations(normCounts, os.path.join(plotsDir, "fig1e-correlations.histones.pdf"), method="single", values=True)


# Fig 1g: ChIPmentation TF correlation
# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("CTCF|GATA1|PU1|REST")) &
    (samples["numberCells"].str.contains("10M|500K|100K")) &
    (~samples["treatment"].notnull()) &
    (samples["technicalReplicate"] == 0)
].reset_index(drop=True)

# get counts per window
normCounts = getCounts(sampleSubset)
# save matrix
c = normCounts.corr()
c['sample1'] = c.index
c = pd.melt(c, id_vars=['sample1'])
c.to_csv(os.path.join(plotsDir, "fig1e-correlations.tfs.tsv"))
plotCorrelations(normCounts, os.path.join(plotsDir, "fig1e-correlations.tfs.pdf"), method="single", values=True)



# Fig S2b: ChIP-tagmentation correlation
# subset samples
sampleSubset = samples[
    samples["sampleName"].str.contains("PBMC_nan_ATAC_H3K4ME3_nan_100PG_1_0_hg19|" +
    "PBMC_nan_ATAC_H3K4ME3_nan_100PG_1_1_hg19|" +
    "PBMC_nan_ATAC_H3K4ME3_nan_100PG_1_2_hg19|" +
    "PBMC_nan_CHIP_H3K4ME3_nan_nan_1_1_hg19")
].reset_index(drop=True)

# get counts per window
normCounts = getCounts(sampleSubset)
# save matrix
c = normCounts.corr()
c['sample1'] = c.index
c = pd.melt(c, id_vars=['sample1'])
c.to_csv(os.path.join(plotsDir, "fig1e-correlations.chip-tagmentation.tsv"))
plotCorrelations(normCounts, os.path.join(plotsDir, "fig1e-correlations.chip-tagmentation.pdf"), method="single", values=True)


for s1, s2 in itertools.combinations(sampleSubset['sampleName']):
    plotScatterCorrelation(
        normCounts[[s1, s2]],
        os.path.join(plotsDir, "figs2b-correlations.chip-tagmentation.scatter.{0}-{1}.pdf".format(s1, s2))
    )
