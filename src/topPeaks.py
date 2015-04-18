#!/usr/env python

import os
import re
import pandas as pd
import itertools  # for R dataframe conversion

import matplotlib.pyplot as plt
import seaborn as sns


def getTopPeakOverlap(a, b, perc=1):
    """
    Gets fraction of top peaks in sample A in top peaks of sample B.
    """
    import subprocess
    import math

    broadFactors = [
        "H3K27ME1", "H3K27ME2", "H3K27ME3",
        "H3K36ME1", "H3K36ME2", "H3K36ME3",
        "H3K9ME1", "H3K9ME2", "H3K9ME3",
        "H3K72ME1", "H3K72ME2", "H3K72ME3"
    ]

    peakA = a['peakFile'].reset_index(drop=True)[0]
    peakB = b['peakFile'].reset_index(drop=True)[0]

    peakA = re.sub("/fhgfs/groups/lab_bock/shared/projects/", "/media/afr/cemm-backup/", peakA)
    peakB = re.sub("/fhgfs/groups/lab_bock/shared/projects/", "/media/afr/cemm-backup/", peakB)

    if a.ip.values[0] not in broadFactors:
        topPeakA = re.sub("narrowPeak", "sorted.narrowPeak", peakA)
    else:
        topPeakA = re.sub("broadPeak", "sorted.broadPeak", peakA)
    if b.ip.values[0] not in broadFactors:
        topPeakB = re.sub("narrowPeak", "sorted.narrowPeak", peakB)
    else:
        topPeakB = re.sub("broadPeak", "sorted.broadPeak", peakB)

    # get total (to get top 1%)
    proc = subprocess.Popen(["wc", "-l", peakA], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    totalA = re.sub("\D.*", "", out)

    proc = subprocess.Popen(["wc", "-l", peakB], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    totalB = re.sub("\D.*", "", out)

    frac = 100 / perc

    topA = str(math.trunc(int(totalA) / frac))
    topB = str(math.trunc(int(totalB) / frac))

    # sort files by score and get top 1%
    ps = subprocess.Popen(('sort', '-k9rn', peakA), stdout=subprocess.PIPE)
    output = subprocess.check_output(('head', '-n', topA), stdin=ps.stdout)

    with open(topPeakA, 'w') as handle:
        handle.write(output)

    ps = subprocess.Popen(('sort', '-k9rn', peakB), stdout=subprocess.PIPE)
    output = subprocess.check_output(('head', '-n', topB), stdin=ps.stdout)

    with open(topPeakB, 'w') as handle:
        handle.write(output)

    # intersect top peaks
    proc = subprocess.Popen(
        ["bedtools", "intersect", "-u", "-a", topPeakA, "-b", topPeakB],
        stdout=subprocess.PIPE
    )
    out, err = proc.communicate()

    # return count
    try:
        print(a.sampleName, b.sampleName)
        print(len(out.split("\n")) / float(topA))
        return len(out.split("\n")) / float(topA)
    except ZeroDivisionError:
        return 0


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

overlaps = pd.DataFrame(index=samples.sampleName.unique(), columns=samples.sampleName.unique())

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
            ctrl = ["IGG", "INPUT"]
            if s1.empty or s2.empty or s1.ip.values[0] in ctrl or s2.ip.values[0] in ctrl:
                continue

            overlaps.loc[s1.sampleName.values[0], s2.sampleName.values[0]] = getTopPeakOverlap(s1, s2)
            overlaps.loc[s2.sampleName.values[0], s1.sampleName.values[0]] = getTopPeakOverlap(s2, s1)

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
            if s1.empty or s2.empty or s1.ip.values[0] in ctrl or s2.ip.values[0] in ctrl:
                continue

            overlaps.loc[s1.sampleName.values[0], s2.sampleName.values[0]] = getTopPeakOverlap(s1, s2)
            overlaps.loc[s2.sampleName.values[0], s1.sampleName.values[0]] = getTopPeakOverlap(s2, s1)

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
            if s1.empty or s2.empty or s1.ip.values[0] in ctrl or s2.ip.values[0] in ctrl:
                continue

            overlaps.loc[s1.sampleName.values[0], s2.sampleName.values[0]] = getTopPeakOverlap(s1, s2)
            overlaps.loc[s2.sampleName.values[0], s1.sampleName.values[0]] = getTopPeakOverlap(s2, s1)

# Plot
# col/row colours
colours = map(colourPerFactor, overlaps.index)

# data colour map
cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)

sns.clustermap(overlaps.fillna(0), row_colors=colours, method="ward", metric="euclidean",
               col_colors=colours, figsize=(15, 15), cmap=cmap)

plt.savefig(os.path.join(plotsDir, "topPeakOverlaps.pdf"), bbox_inches='tight')
