#!/usr/env python

import os
import pandas as pd
import matplotplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns
sns.set_style("whitegrid")


def millions(x, pos):
    'The two args are the value and tick position'
    return '%1.1f' % (x * 1e-6)


formatter = FuncFormatter(millions)

root = "/media/afr/cemm-backup/dyads"

dists = pd.read_csv(os.path.join(root, "distogram.txt"), sep="\t", header=None)

fig, ax = plt.subplots()
plt.plot(dists[0], dists[1])
plt.axvline(x=147, linestyle='--', color="grey")
plt.xticks(list(plt.xticks()[0]) + [147] + range(0, 1500, 250))
ax.yaxis.set_major_formatter(formatter)
plt.xlabel("Distance (bp)")
plt.ylabel("Distance counts (millions)")

plt.xlim([-20, 1500])
plt.ylim([1450000, 2000000])

plt.savefig(os.path.join(root, "distogram.pdf"))


phase = pd.read_csv(os.path.join(root, "phasogram.txt"), sep="\t", header=None)

fig, ax = plt.subplots()
plt.plot(phase[0], phase[1])
plt.axvline(x=185, linestyle='--', color="grey")
plt.xticks(list(plt.xticks()[0]) + [185] + range(0, 1500, 250))
ax.yaxis.set_major_formatter(formatter)
plt.xlabel("Phase (bp)")
plt.ylabel("Phase counts (millions)")

plt.xlim([-20, 1500])
plt.ylim([2900000, 3900000])

plt.savefig(os.path.join(root, "phasogram.pdf"))
