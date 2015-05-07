#!/usr/env python

import os
import re
import pandas as pd
import itertools  # for R dataframe conversion

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")


def getTopPeakOverlap(a, b, percA=100, percB=100):
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

    fracA = 100 / percA
    fracB = 100 / percB

    topA = str(math.trunc(int(totalA) / fracA))
    topB = str(math.trunc(int(totalB) / fracB))

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


samples = samples[samples.ip.notnull()]
samples = samples[samples.ip.str.contains("PU1|CTCF|GATA1|REST")]


overlaps = pd.DataFrame()

for threshold in [1, 3, 5, 10, 12, 18, 25, 50, 100]:
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

                o = max(getTopPeakOverlap(s1, s2, threshold), getTopPeakOverlap(s2, s1, threshold, 100))
                s = pd.Series(
                    [ip, c, c, t1, t2, 0, 0, o, 'techniques', threshold],
                    index=['ip', 'c1', 'c2', 't1', 't2', 'r1', 'r2', 'overlap', 'type', 'threshold'])
                overlaps = overlaps.append(s, ignore_index=True)

    print("cells")
    # Different number of cells
    for ip in samples['ip'].unique():
        for t in ["CM"]:
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

                o = max(getTopPeakOverlap(s1, s2, threshold), getTopPeakOverlap(s2, s1, threshold, 100))
                s = pd.Series(
                    [ip, c1, c2, t, t, 0, 0, o, 'cells', threshold],
                    index=['ip', 'c1', 'c2', 't1', 't2', 'r1', 'r2', 'overlap', 'type', 'threshold'])
                overlaps = overlaps.append(s, ignore_index=True)

    print("replicates")
    # Replicate 1 vs 2
    for ip in samples['ip'].unique():
        for c in samples['numberCells'].unique():
            for t in ["CM"]:
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

                o = max(getTopPeakOverlap(s1, s2, threshold), getTopPeakOverlap(s2, s1, threshold, 100))
                s = pd.Series(
                    [ip, c, c, t, t, r1, r2, o, 'replicates', threshold],
                    index=['ip', 'c1', 'c2', 't1', 't2', 'r1', 'r2', 'overlap', 'type', 'threshold'])
                overlaps = overlaps.append(s, ignore_index=True)

# Remove replicates overlap between low cell number points
oo = overlaps[~(
    (overlaps['type'] == "replicates") &
    (overlaps['c1'] != "10M") & (overlaps['c2'] != "10M")
)]

# set > 1 to 1 (this is due to multiple overlapping)
oo.loc[oo['overlap'] > 1,'overlap'] = 1

# Plot
colours = {ip: colourPerFactor(ip) for ip in overlaps['ip'].unique()}

# All
g = sns.FacetGrid(overlaps, col="type", hue="ip")  # , hue_kws=colours)
# g.map(plt.scatter, "threshold", "overlap")
# g.map(plt.plot, "threshold", "overlap")
g.map(sns.pointplot, "threshold", "overlap")
g.add_legend()
plt.savefig(os.path.join(plotsDir, "topPeakOverlaps.all.pdf"), bbox_inches='tight')

# Major
g = sns.FacetGrid(oo, col="type", hue="ip")
g.map(sns.pointplot, "threshold", "overlap")
g.add_legend()
plt.savefig(os.path.join(plotsDir, "topPeakOverlaps.pdf"), bbox_inches='tight')

# Subset
oo = oo[oo['threshold'] < 50]

g = sns.FacetGrid(oo, col="type", hue="ip")
g.map(sns.pointplot, "threshold", "overlap")
g.add_legend()

plt.savefig(os.path.join(plotsDir, "topPeakOverlaps.subset.pdf"), bbox_inches='tight')



# Pie charts with 5 and 25 thresholds
ooo = oo[oo['threshold'] == 5]
g = sns.FacetGrid(ooo, col="type", hue="ip")
g.map(plt.scatter, "threshold", "overlap")


ooo = oo[oo['threshold'] == 25]


