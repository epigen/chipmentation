import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_style("whitegrid")


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


df = pd.read_csv("chipmentation.stats.csv")
df2 = df[['numberCells', 'technique', 'ip', 'NSC', 'RSC', 'peakNumber', 'FRiP']]

df2 = df2.dropna()
df2 = df2[df2['FRiP'] < 1.05]

df2.loc[df2['FRiP'] > 1, "FRiP"] -= 0.05

f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
for tech in df2.groupby(['technique']):
    for ip in tech[1].groupby(['ip']):
        if tech[0] == "CM":
            ax1.plot(ip[1]['peakNumber'], ip[1]['FRiP'], 'o', label=ip[0], linestyle='--', color=colourPerFactor(ip[0]))
        else:
            ax2.plot(ip[1]['peakNumber'], ip[1]['FRiP'], 'o', label=ip[0], linestyle='--', color=colourPerFactor(ip[0]))
plt.xlabel("Number of called peaks")
plt.ylabel("FRiP")
plt.legend(loc="best")
plt.savefig("numPeaks_FRiP.pdf", bbox_inches='tight')
plt.close()


f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
for tech in df2.groupby(['technique']):
    for ip in tech[1].groupby(['ip']):
        if tech[0] == "CM":
            ax1.plot(np.log10(ip[1]['peakNumber']), ip[1]['FRiP'], 'o', label=ip[0], linestyle='--', color=colourPerFactor(ip[0]))
        else:
            ax2.plot(np.log10(ip[1]['peakNumber']), ip[1]['FRiP'], 'o', label=ip[0], linestyle='--', color=colourPerFactor(ip[0]))
plt.xlabel("Number of called peaks (log10)")
plt.ylabel("FRiP")
plt.legend(loc="best")
plt.savefig("numPeaks_FRiP.log10.pdf", bbox_inches='tight')
plt.close()


# # using grid
# grid = sns.FacetGrid(df2, col="technique", hue="ip", sharex=True, col_wrap=4)
# grid.map(plt.plot, "peakNumber", "FRiP")
# grid.fig.subplots_adjust(wspace=0.5, hspace=0.5)
# grid.add_legend()
# # plt.savefig("frip.pdf", bbox_inches='tight')
