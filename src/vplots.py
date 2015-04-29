import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("white")

df = pd.read_csv("/media/afr/cemm-backup/chipmentation/data/nucleoATAC/10M_CM_H3K4ME1_PE/dyads_vplot.tsv", sep="\t", header=None)
df = df.sort(ascending=False)
df.columns = np.array(range(len(df.columns))) - (len(df.columns) / 2)

df = np.log2(1 + df)

# cmap = plt.get_cmap("BuGn")
cmap = plt.get_cmap("pink_r")
plt.imshow(df, cmap=cmap)
plt.yticks(plt.yticks()[0][2:-1], [200, 150, 100, 50])
plt.xticks(plt.xticks()[0][1:-1], [-200, -100, 0, 100, 200])

plt.vlines(200 - 60, 0, 250, linestyles="--")
plt.vlines(200 + 60, 0, 250, linestyles="--")

plt.colorbar()
plt.xlabel("distance to dyad (bp)")
plt.ylabel("fragment size (bp)")

plt.savefig("dyads_vplot2.pdf", bbox_inches='tight')
plt.close()

# Zoom on nucleosomal part
plt.imshow(df.loc[250:105, -60:60], cmap=cmap)
plt.colorbar()
plt.xlabel("distance to dyad (bp)")
plt.ylabel("fragment size (bp)")
plt.xticks([0, 30, 60, 90, 120], [-60, -30, 0, 30, 60])
plt.yticks([0, 50, 100, 140], [250, 200, 150, 110])
plt.savefig("dyads_vplot.zoom2.pdf", bbox_inches='tight')
plt.close()


# Insertion frequency
sns.set_style("whitegrid")
df = pd.read_csv("/media/afr/cemm-backup/chipmentation/data/nucleoATAC/10M_CM_H3K4ME1_PE/10M_CM_H3K4ME1_PE.vplot.InsertionProfile.txt", header=None)
df["x"] = range(-325, 326)

plt.plot(df["x"], df[0], color="#e69f00")

plt.xlabel("distance to dyad (bp)")
plt.ylabel("insertion frequency")
plt.ylim([0, 0.011])
plt.savefig("dyads_insertion_frequency.pdf", bbox_inches='tight')
plt.close()
