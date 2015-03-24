#!/usr/bin/env python

import os
import pickle
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind

import pandas as pd

plotsDir = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/results/pe/plots"
samples = pd.read_csv("/fhgfs/groups/lab_bock/shared/projects/chipmentation/chipmentation.read_stats.csv")


# cmd = """
# getFragments() {
#     NAME=$1
#     samtools view /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.bam | \
#     python /home/arendeiro/fragmentLengths.py $NAME
# }
# export -f getFragments

# SAMPLES=(K562_10M_ATAC_H3K4ME1_nan_PE_1_1_hg19
#     K562_10M_ATAC_PU1_nan_PE_1_1_hg19
#     K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19
#     K562_10M_CM_PU1_nan_PE_1_1_hg19
#     K562_500K_ATAC_H3K4ME3_nan_01ULTN5_PE_1_1_hg19
#     K562_500K_CM_H3K4ME3_nan_02ULTN5_PE_1_1_hg19
#     K562_500K_CM_H3K4ME3_nan_05ULTN5_PE_1_1_hg19
#     K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19
#     K562_500K_CM_H3K4ME3_nan_5ULTN5_PE_1_1_hg19
#     K562_500K_CM_IGG_nan_1ULTN5_PE_1_1_hg19)
# parallel getFragments ::: ${SAMPLES[*]}
# """

# os.system(cmd)

samples = [
    "K562_10M_ATAC_H3K4ME1_nan_PE_1_1_hg19", "K562_10M_ATAC_PU1_nan_PE_1_1_hg19",
    "K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19", "K562_10M_CM_PU1_nan_PE_1_1_hg19",
    "K562_500K_ATAC_H3K4ME3_nan_01ULTN5_PE_1_1_hg19", "K562_500K_CM_H3K4ME3_nan_02ULTN5_PE_1_1_hg19",
    "K562_500K_CM_H3K4ME3_nan_05ULTN5_PE_1_1_hg19", "K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19",
    "K562_500K_CM_H3K4ME3_nan_5ULTN5_PE_1_1_hg19", "K562_500K_CM_IGG_nan_1ULTN5_PE_1_1_hg19"
]

# Plot all samples individually
for i in range(len(samples)):
    name = samples[i]
    print(name)
    dist = pickle.load(open(name + ".pickle", "r"))
    total = sum(dist.values())

    # remove "0" (signletons) and values higher than 300bp
    # this is different than saying `plt.xlim(0, 300)` which still saves it into the pdf
    # somehow as transparent lines or something, Christian told me
    dist.pop(0)
    [dist.pop(k) for k in dist.keys() if k > 300]

    plt.plot(dist.keys(), (pd.Series(dist.values()) / total) * 100, '-', label=name)

    plt.legend(loc='upper right', shadow=False)
    # plt.xlim(0, 300)
    plt.ylabel(r'Normalized read density x $10^-3$')

    plt.savefig(os.path.join(plotsDir, name + "_fragmentDistribution.pdf"))
    plt.close()

    # Get percentiles
    # import numpy as np
    # l = [size for size in dist.keys() for _ in range(dist[size])]
    # print(np.percentile(l, [20, 40, 60, 80]))


# Plot all samples on one plot
for i in range(len(samples)):
    name = samples[i]
    dist = pickle.load(open(name + ".pickle", "r"))
    total = sum(dist.values())

    dist.pop(0)
    [dist.pop(k) for k in dist.keys() if k > 300]

    plt.plot(dist.keys(), (pd.Series(dist.values()) / total) * 100, '-', label=name)

plt.legend(loc='upper right', shadow=False)
plt.ylabel(r'Normalized read density x $10^-3$')
plt.savefig(os.path.join(plotsDir, "all" + "_fragmentDistribution.pdf"))
plt.close()


# Plot Only H3K4me3
subset = [sample for sample in samples if "H3K4ME3" in sample]

for i in range(len(subset)):
    name = subset[i]
    dist = pickle.load(open(name + ".pickle", "r"))
    total = sum(dist.values())
    
    dist.pop(0)
    [dist.pop(k) for k in dist.keys() if k > 300]

    plt.plot(dist.keys(), (pd.Series(dist.values()) / total) * 100, '-', label=name)

plt.legend(loc='upper right', shadow=False)
plt.ylabel(r'Normalized read density x $10^-3$')
plt.savefig(os.path.join(plotsDir, "H3K4ME3" + "_fragmentDistribution.pdf"))
plt.close()


# Plot Only H3K4me1
subset = [sample for sample in samples if "H3K4ME1" in sample]

dists = list()
for i in range(len(subset)):
    name = subset[i]
    dist = pickle.load(open(name + ".pickle", "r"))
    total = sum(dist.values())

    dist.pop(0)
    [dist.pop(k) for k in dist.keys() if k > 300]

    plt.plot(dist.keys(), (pd.Series(dist.values()) / total) * 100, '-', label=name)

    dists.append([size for size in dist.keys() for _ in range(dist[size])])

plt.legend(loc='upper right', shadow=False)
t, p = ttest_ind(dists[0], dists[1])
plt.ylabel(r'Normalized read density x $10^-3$')
# plt.text(200, 300000, "p-value: %f" % p, fontsize=12)
plt.savefig(os.path.join(plotsDir, "H3K4ME1" + "_fragmentDistribution.pdf"))
plt.close()


# Plot Only PU1
subset = [sample for sample in samples if "PU1" in sample]

dists = list()
for i in range(len(subset)):
    name = subset[i]
    dist = pickle.load(open(name + ".pickle", "r"))
    total = sum(dist.values())

    dist.pop(0)
    [dist.pop(k) for k in dist.keys() if k > 300]

    dists.append([size for size in dist.keys() for _ in range(dist[size])])

    plt.plot(dist.keys(), (pd.Series(dist.values()) / total) * 100, '-', label=name)

plt.legend(loc='upper right', shadow=False)
t, p = ttest_ind(dists[0], dists[1])
# plt.text(200, 300000, "p-value: %f" % p, fontsize=12)
plt.savefig(os.path.join(plotsDir, "PU1" + "_fragmentDistribution.pdf"))
plt.close()
