#!/usr/env python

import os
import re
import numpy as np
import pandas as pd
import itertools
import rpy2.robjects as robj  # for ggplot in R
import rpy2.robjects.pandas2ri  # for R dataframe conversion

import matplotlib.pyplot as plt
import seaborn as sns


def getTopPeaks(a, b):
    import subprocess
    import math

    peakA = a['peakFile'].reset_index(drop=True)[0]
    peakB = b['peakFile'].reset_index(drop=True)[0]

    peakA = re.sub("/fhgfs/groups/lab_bock/shared/projects/", "/media/afr/cemm-backup/chipmentation/", peakA)
    peakB = re.sub("/fhgfs/groups/lab_bock/shared/projects/", "/media/afr/cemm-backup/chipmentation/", peakB)

    topPeakA = re.sub("narrowPeak", "sorted.narrowPeak", peakA)
    topPeakB = re.sub("narrowPeak", "sorted.narrowPeak", peakB)

    # get total (to get top 1%)
    proc = subprocess.Popen(["wc", "-l", peakA], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    totalA = re.sub("\D.*", "", out)

    proc = subprocess.Popen(["wc", "-l", peakB], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    totalB = re.sub("\D.*", "", out)

    topA = math.trunc(totalA / 100)
    topB = math.trunc(totalB / 100)

    # sort files by score and get top 1%
    proc = subprocess.Popen(
        ["sort", "-k9", "-n", peakA, "|", "head", "-n", topA, ">", topPeakA],
        stdout=subprocess.PIPE, shell=True
    )
    out, err = proc.communicate()
    proc = subprocess.Popen(
        ["sort", "-k9", "-n", peakB, "|", "head", "-n", topB, ">", topPeakB],
        stdout=subprocess.PIPE, shell=True
    )
    out, err = proc.communicate()

    # intersect bed
    proc = subprocess.Popen(
        ["bedtools", "intersect", "-c", "-a", topPeakA, "-b", topPeakB],
        stdout=subprocess.PIPE, shell=True
    )
    out, err = proc.communicate()

    intersect = re.sub("\D.*", "", out)

    return intersect



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

            q = (s1, s2)


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

             = getTopPeaks(s1, s2)
