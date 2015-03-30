#!/usr/env python

#############################################################################################
#
# This code was used to produce the plots in the ChIPmentation paper (Schmidl, et al. 2015).
#
#############################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define variables
projectRoot = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/"
resultsDir = projectRoot + "results"
plotsDir = resultsDir + "/plots"

# Get samples
samples = pd.read_csv(os.path.abspath(projectRoot + "chipmentation.read_stats.csv"))

# subset samples
sampleSubset = samples[
    (samples["technique"].str.contains("CM|ATAC")) &
    (samples["ip"].str.contains("H3K4ME3")) &
    (samples["readType"].str.contains("PE"))
].reset_index(drop=True)

sampleSubset['Tn5'] = [0.1, 1, 5, 0.5, 0.2]
sampleSubset['percentAligned'] = [98.07, 98.16, 98.22, 98.06, 97.98]
sampleSubset['percentUnmapped'] = 100 - sampleSubset['percentAligned']
sampleSubset.sort(['Tn5'], inplace=True)

# plot
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(sampleSubset['Tn5'], sampleSubset['percentDuplication'], 'og-', label="% duplicates")
ax2.plot(sampleSubset['Tn5'], sampleSubset['percentUnmapped'], 'ob-', label="% unmapped")

ax1.set_xlabel(r'Tn5 (uL)')
ax1.set_ylabel('% duplicates', color='g')
ax2.set_ylabel('% unmapped', color='b')

ax1.set_xticks(sorted(sampleSubset['Tn5'].tolist()))
ax1.set_xticklabels(sampleSubset['Tn5'].tolist())

plt.xlim(0, 5.5)
ax1.set_ylim(0, 100)
ax2.set_ylim(0, 100)

plt.savefig(os.path.join(plotsDir, "Tn5_titration_percentAlignedUnique.all.pdf"))
plt.close()

# plot only ChIPmentation
sampleSubset[sampleSubset.technique != "ATAC"]

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(sampleSubset['Tn5'], sampleSubset['percentDuplication'], 'og-', label="% duplicates")
ax2.plot(sampleSubset['Tn5'], sampleSubset['percentUnmapped'], 'ob-', label="% unmapped")

ax1.set_xlabel(r'Tn5 (uL)')
ax1.set_ylabel('% duplicates', color='g')
ax2.set_ylabel('% unmapped', color='b')

ax1.set_xticks(sorted(sampleSubset['Tn5'].tolist()))
ax1.set_xticklabels(sampleSubset['Tn5'].tolist())

plt.xlim(0, 5.5)
ax1.set_ylim(0, 100)
ax2.set_ylim(0, 100)

plt.savefig(os.path.join(plotsDir, "Tn5_titration_percentAlignedUnique.pdf"))
plt.close()

# plot with log scale
fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(np.log2(sampleSubset['Tn5']), sampleSubset['percentDuplication'], 'og-', label="% duplicates")
ax2.plot(np.log2(sampleSubset['Tn5']), sampleSubset['percentUnmapped'], 'ob-', label="% unmapped")

ax1.set_xlabel(r'log2(Tn5) (uL)')
ax1.set_ylabel('% duplicates', color='g')
ax2.set_ylabel('% unmapped', color='b')

ax1.set_xticks(sorted(np.log2(sampleSubset['Tn5'].tolist())))
ax1.set_xticklabels(np.round(np.log2(sampleSubset['Tn5'].tolist()), 2))

plt.savefig(os.path.join(plotsDir, "Tn5_titration_percentAlignedUnique.log2.pdf"))
plt.close()
