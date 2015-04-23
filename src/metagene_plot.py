#!/usr/env python

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.theme("whitegrid")


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
    elif "DNASE" in name:
        return "#005212"
    elif "MNASE" in name:
        return "#00523b"
    elif "ATAC" in name:
        return "#001f07"
    else:
        raise ValueError

rootDir = "/home/afr/workspace/chipmentation/"

# preliminary
df = pd.read_csv(os.path.join(rootDir, "confiMat.csv"))
df = df.drop(["Dnase"], axis=1)
df = df.drop(["MNase"], axis=1)
df = df.drop(["Unnamed: 0"], axis=1)

# plot
fig = plt.figure()
for name in df.columns:
    plt.plot(
        df[name],
        label=name,
        color=colourPerFactor(name.split("_")[-1]),
        linestyle="-" if ("CM" in name) or ("ATAC" in name) else "--",
        linewidth=2.5
    )
plt.xticks(list(plt.xticks()[0]), [-2000, "TSS", "33%", "66%", "TES", 2000])
plt.ylabel("Read count per million mapped reads")
plt.legend(loc="best")
fig.set_size_inches(10, 4)
plt.savefig(os.path.join(rootDir, "results", "plots", "metagene2.pdf"), bbox_inches='tight')


# final
df = pd.read_csv(os.path.join(rootDir, "regcovMat.csv"))
df = df.drop(["Dnase"], axis=1)
df = df.drop(["MNase"], axis=1)
df = df.drop(["Unnamed: 0"], axis=1)

# plot
fig = plt.figure()
for name in df.columns:
    plt.plot(
        df[name],
        label=name,
        color=colourPerFactor(name.split("_")[-1]),
        linestyle="-" if ("CM" in name) or ("ATAC" in name) else "--",
        linewidth=2.5
    )
plt.xticks(list(plt.xticks()[0]), [-2000, "TSS", "33%", "66%", "TES", 2000])
plt.ylabel("Read count per million mapped reads")
plt.legend(loc="best")
fig.set_size_inches(10, 4)
plt.savefig(os.path.join(rootDir, "results", "plots", "metagene.pdf"), bbox_inches='tight')
