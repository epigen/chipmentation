#!/usr/bin/env python

#############################################################################################
#
# This code was used to produce the plots in the ChIPmentation paper (Schmidl, et al. 2015).
#
# It produces plots of sample correlation over windows genome-wide.
#
#############################################################################################

import pandas as pd
import os
import pysam
import matplotlib.pyplot as plt
import seaborn as sns


sns.set_style("whitegrid")


def getFragmentLenths(bam):
    from collections import Counter

    count = Counter()

    bam = pysam.AlignmentFile(bam)

    try:
        for read in bam:
            count[abs(read.tlen)] += 1
    except KeyboardInterrupt:
        return count

    return count


dataDir = "/media/afr/cemm-backup1/chipmentation/data/mapped/"
plotsDir = "/media/afr/cemm-backup1/chipmentation/plots"
samples = pd.read_csv("workspace/chipmentation/chipmentation.replicates.annotation_sheet.csv")

sampleNames = [
    "K562_500K_CM_H3K4ME3_nan_02ULTN5_PE_1_1_hg19", "K562_500K_CM_H3K4ME3_nan_05ULTN5_PE_1_1_hg19",
    "K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19", "K562_500K_CM_H3K4ME3_nan_5ULTN5_PE_1_1_hg19"
]

dists = pd.DataFrame()

# Plot all samples individually
fig = plt.figure(1)
for i, name in enumerate(sampleNames):
    print(name)
    bam = os.path.join(dataDir, name + ".trimmed.bowtie2.shifted.dups.bam")

    dist = getFragmentLenths(bam)
    # remove "0" (signletons) and values higher than 300bp
    # this is different than saying `plt.xlim(0, 300)` which still saves it into the pdf
    # somehow as transparent lines or something, Christian told me
    dist.pop(0)
    [dist.pop(k) for k in dist.keys() if k > 300]

    dists[name] = pd.Series(dist)

    plt.plot(dist.keys(), (pd.Series(dist.values()) / sum(dist)) * 100, '-', label=name)

dists.to_csv(os.path.join(plotsDir, "fig1b_fragmentDistribution.csv"))

plt.legend(loc='upper right', shadow=False)
plt.ylabel(r'Normalized read density x $10^-3$')

plt.savefig(os.path.join(plotsDir, "fig1b_fragmentDistribution.pdf"))
plt.close()
