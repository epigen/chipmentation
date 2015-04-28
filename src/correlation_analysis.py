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

# import matplotlib
# # Force matplotlib to not use any Xwindows backend.
# matplotlib.use('Agg')
# import matplotlib.font_manager as font_manager

# fontpath = '/usr/share/fonts/truetype/Roboto-Regular.ttf'

# prop = font_manager.FontProperties(fname=fontpath)
# matplotlib.rcParams['font.family'] = prop.get_name()

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
        raise ValueError


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    import matplotlib.colors as colors
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def plotCorrelations(normCounts, plotName, method="ward", metric="euclidean", values=False):
    corr = normCounts.corr()

    # col/row colours
    colours = map(colourPerFactor, normCounts.columns)

    # data colour map
    cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
    # cmap = plt.get_cmap("YlGn")

    # scale from 0 to 1
    new_cmap = truncate_colormap(cmap, corr.min().min(), 1)

    if not values:
        sns.clustermap(corr, row_colors=colours, method=method, metric=metric,
                       col_colors=colours, figsize=(15, 15), cmap=new_cmap)
    else:
        sns.clustermap(corr, row_colors=colours, method=method, metric=metric,
                       col_colors=colours, figsize=(15, 15), cmap=new_cmap, annot=True)

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


def plotKDE(s1, s2, path):
    sns.jointplot(s1, s2, kind="kde", size=7, space=0)
    plt.savefig(path, bbox_inches='tight')


# Define paths
# projectRoot = "/projects/chipmentation/"
projectRoot = "/media/afr/cemm-backup/chipmentation/"
# projectRoot = "/home/arendeiro/chipmentation/"
coverageDir = projectRoot + "data/coverage"
resultsDir = projectRoot + "results"
plotsDir = resultsDir + "/plots/correlations"

# Get samples
samples = pd.read_csv(os.path.abspath(projectRoot + "chipmentation.replicates.annotation_sheet.csv"))

# remove missing sample
samples = samples[samples['sampleName'] != "K562_10M_CHIP_H3K4ME1_nan_nan_2_1_hg19"]
samples = samples[samples['sampleName'] != "K562_500K_CHIP_H3K27AC_nan_nan_1_1_hg19"]

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

sampleSubset = sampleSubset[pd.isnull(sampleSubset['treatment'])].reset_index(drop=True)

sampleSubset = sampleSubset.append(
    pd.Series(
        data=["K562_10M_CHIP_H3K36ME3_nan_nan_1_1_hg19", "H3K36ME3", "CHIP"],
        index=["sampleName", "ip", "technique"]
    ),
    ignore_index=True
)
sampleSubset = sampleSubset.append(
    pd.Series(
        data=["K562_10M_CHIP_H3K27AC_nan_nan_1_0_hg19", "H3K27AC", "CHIP"],
        index=["sampleName", "ip", "technique"]
    ),
    ignore_index=True
)

sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)
if os.path.exists(os.path.join(plotsDir, "correlations.set.pickle")):
    normCounts = pd.read_pickle(os.path.join(plotsDir, "correlations.set.pickle"))
else:
    normCounts = getCounts(sampleSubset)
    normCounts.to_pickle(os.path.join(plotsDir, "correlations.set.pickle"))

plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.set.pdf"))
plotCorrelations(normCounts.filter(regex="CM"), os.path.join(plotsDir, "correlations.set.CM.pdf"), method="complete", metric="euclidean")


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

sampleSubset = sampleSubset[pd.isnull(sampleSubset['treatment'])].reset_index(drop=True)

sampleSubset = sampleSubset.append(
    pd.Series(
        data=["K562_10M_CHIP_H3K36ME3_nan_nan_1_1_hg19", "H3K36ME3", "CHIP"],
        index=["sampleName", "ip", "technique"]
    ),
    ignore_index=True
)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

if os.path.exists(os.path.join(plotsDir, "correlations.histones-reps.pickle")):
    normCounts = pd.read_pickle(os.path.join(plotsDir, "correlations.histones-reps.pickle"))
else:
    normCounts = getCounts(sampleSubset)
    normCounts.to_pickle(os.path.join(plotsDir, "correlations.histones-reps.pickle"))

# a = normCounts.index.values
# idx = np.array([a, a, a, a, a]).T.flatten()[:len(a)]
# normCounts = normCounts.groupby(idx).mean()
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.histones-reps.pdf"), method="single")
normCounts = normCounts.filter(regex='10M|500K')
normCounts = normCounts.filter(regex='nan_nan')
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.histones-reps.best.pdf"), method="single")

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

sampleSubset = sampleSubset[pd.isnull(sampleSubset['treatment'])].reset_index(drop=True)

sampleSubset = sampleSubset.append(
    pd.Series(
        data=["K562_10M_CHIP_H3K36ME3_nan_nan_1_1_hg19", "H3K36ME3", "CHIP"],
        index=["sampleName", "ip", "technique"]
    ),
    ignore_index=True
)
sampleSubset = sampleSubset.append(
    pd.Series(
        data=["K562_10M_CHIP_H3K27AC_nan_nan_1_0_hg19", "H3K27AC", "CHIP"],
        index=["sampleName", "ip", "technique"]
    ),
    ignore_index=True
)

sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

if os.path.exists(os.path.join(plotsDir, "correlations.histones.pickle")):
    normCounts = pd.read_pickle(os.path.join(plotsDir, "correlations.histones.pickle"))
else:
    normCounts = getCounts(sampleSubset)
    normCounts.to_pickle(os.path.join(plotsDir, "correlations.histones.pickle"))

# a = normCounts.index.values
# idx = np.array([a, a, a, a, a]).T.flatten()[:len(a)]
# normCounts = normCounts.groupby(idx).mean()
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.histones.pdf"), method="single")
plotCorrelations(normCounts.filter(regex='10M|500K'), os.path.join(plotsDir, "correlations.histones.best.pdf"), method="single")
plotCorrelations(normCounts.filter(regex='CM'), os.path.join(plotsDir, "correlations.histones.CM.pdf"), method="single")


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

if os.path.exists(os.path.join(plotsDir, "correlations.TFs-reps.pickle")):
    normCounts = pd.read_pickle(os.path.join(plotsDir, "correlations.TFs-reps.pickle"))
else:
    normCounts = getCounts(sampleSubset)
    normCounts.to_pickle(os.path.join(plotsDir, "correlations.TFs-reps.pickle"))

# a = normCounts.index.values
# idx = np.array([a, a, a, a, a]).T.flatten()[:len(a)]
# normCounts = normCounts.groupby(idx).mean()
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.TFs-reps.pdf"), method="complete")
plotCorrelations(normCounts.filter(regex='10M'), os.path.join(plotsDir, "correlations.TFs-reps.best.pdf"), method="complete")

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

if os.path.exists(os.path.join(plotsDir, "correlations.TFs.pickle")):
    normCounts = pd.read_pickle(os.path.join(plotsDir, "correlations.TFs.pickle"))
else:
    normCounts = getCounts(sampleSubset)
    normCounts.to_pickle(os.path.join(plotsDir, "correlations.TFs.pickle"))

# a = normCounts.index.values
# idx = np.array([a, a, a, a, a]).T.flatten()[:len(a)]
# normCounts = normCounts.groupby(idx).mean()
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.TFs.pdf"), method="complete")
plotCorrelations(normCounts.filter(regex='10M'), os.path.join(plotsDir, "correlations.TFs.best.pdf"), method="complete")


# H3K4me3 titration
print("Titration")
# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("H3K4ME3")) &
    (samples["numberCells"].str.contains("500K|10K")) &
    (samples["technicalReplicate"] == 0) &
    (samples["biologicalReplicate"] == 0)
].reset_index(drop=True)

s = samples[
    (samples["technique"].str.contains("CM")) &
    (samples["ip"].str.contains("H3K4ME3")) &
    (samples["sampleName"].str.contains("02ULTN5_PE|05ULTN5_PE|1ULTN5_PE|5ULTN5_PE"))
].reset_index(drop=True)

sampleSubset = sampleSubset.append(s).reset_index(drop=True)

sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

if os.path.exists(os.path.join(plotsDir, "correlations.titration.pickle")):
    normCounts = pd.read_pickle(os.path.join(plotsDir, "correlations.titration.pickle"))
else:
    normCounts = getCounts(sampleSubset)
    normCounts.to_pickle(os.path.join(plotsDir, "correlations.titration.pickle"))

# a = normCounts.index.values
# idx = np.array([a, a, a, a, a]).T.flatten()[:len(a)]
# normCounts = normCounts.groupby(idx).mean()
plotCorrelations(normCounts, os.path.join(plotsDir, "correlations.titration.pdf"), method="complete", values=True)


# Scatter plots
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

            d = s1.append(s2).reset_index()
            d = d[pd.isnull(d['treatment'])].reset_index(drop=True)

            normCounts = getCounts(d)
            normCounts.columns = d.sampleName

            pdf = "correlation.{0}.{1}.{2}.{3}_vs_{4}.pdf".format(ip, c, t, r1, r2)

            if not os.path.exists(pdf):
                plotScatterCorrelation(
                    normCounts,
                    os.path.join(
                        plotsDir,
                        pdf
                    )
                )

# Table with TF correlation values (Fig. 1d)
# get all TF samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM")) &
    (samples["ip"].str.contains("CTCF|PU1|GATA1|REST")) &
    (samples["numberCells"].str.contains("10M|500K|100K|10K")) &
    (samples["technicalReplicate"] == 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

# get counts
if os.path.exists(os.path.join(plotsDir, "correlations.table.pickle")):
    normCounts = pd.read_pickle(os.path.join(plotsDir, "correlations.table.pickle"))
else:
    normCounts = getCounts(sampleSubset)
    normCounts.to_pickle(os.path.join(plotsDir, "correlations.table.pickle"))


# group by everything except biological rep
corr = pd.DataFrame()


# get only CM samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM")) &
    (samples["ip"].str.contains("CTCF|PU1|GATA1|REST")) &
    (samples["numberCells"].str.contains("10M|500K|100K|10K")) &
    (samples["biologicalReplicate"] != 0) &
    (samples["technicalReplicate"] == 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)

groups = sampleSubset.groupby(['cellLine', 'numberCells', 'technique', 'ip',  'patient', 'treatment'])

for g in groups.groups.items():
    if g[0][5] == "PE":
        continue
    df = pd.DataFrame()
    for i in g[1]:
        name = sampleSubset.ix[i]['sampleName']
        name = re.sub("_hg19", "", name)
        name = re.sub("K562_", "", name)
        try:
            df = df.append(normCounts[name], ignore_index=True)
        except:
            print("failed sample %s" % name)
    if len(df) == 2:
        corr.loc[name, 'CR'] = df.T.corr()[1][0]

# get only 0_0 samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|CHIP")) &
    (samples["ip"].str.contains("CTCF|PU1|GATA1|REST")) &
    (samples["numberCells"].str.contains("10M|500K|100K|10K")) &
    (samples["biologicalReplicate"] == 0) &
    (samples["technicalReplicate"] == 0)
].reset_index(drop=True)
sampleSubset = sampleSubset.sort(["ip", "technique"]).reset_index(drop=True)


groups = sampleSubset.groupby(['cellLine', 'numberCells', 'ip',  'patient', 'treatment'])

for g in groups.groups.items():
    if len(g[1]) == 2:
        print(g)
        df = pd.DataFrame()
        for i in g[1]:
            name = sampleSubset.ix[i]['sampleName']
            name = re.sub("_hg19", "", name)
            name = re.sub("K562_", "", name)
            print(name)
            try:
                df = df.append(normCounts[name], ignore_index=True)
            except:
                print("failed sample %s" % name)
                break
        if len(df) == 2:
            corr.loc[name, 'CT'] = df.T.corr()[1][0]

# plot
cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
new_cmap = truncate_colormap(cmap, corr.min().min(), 1)

f, ax = plt.subplots()
sns.heatmap(corr.sort(), vmax=.8, linewidths=0, square=True, cmap=new_cmap, annot=True, fmt=".3")
plt.savefig(os.path.join(plotsDir, "correlations.TF.table.pdf"), bbox_inches="tight")
