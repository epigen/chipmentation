import sys
import csv
from matplotlib import pyplot
import pandas as pd
from itertools import imap
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.interpolate import interp1d as interp


sample = sys.argv[1]
sample = 'H3K4me3_K562_500k_CM'
projectDir = '/home/afr/Downloads/data/'


def parseBedCoverage(bedFile):
    """
    Reads tab-delimiter file with two rows describing the chromossomes and its lengths.
    Returns dictionary of chr:sizes.
    """
    prev_peak = ''
    peaks = {}
    with open(bedFile, "rb") as f:
        reader = csv.reader(f, delimiter="\t")
        for line, row in enumerate(reader):
            peak = str(row[3])
            bp = int(row[12])
            count = int(row[13])
            strand = str(row[5])
            if peak != prev_peak:
                # new peak
                peaks[peak] = {}
                if strand == "+":
                    peaks[peak][bp - 2001] = count
                elif strand == "-":
                    peaks[peak][2001 - bp] = count
                prev_peak = peak
            else:
                if strand == "+":
                    peaks[peak][bp - 2001] = count
                elif strand == "-":
                    peaks[peak][2001 - bp] = count
    sortedPeaks = {}
    for key in sorted(peaks.iterkeys()):
        sortedPeaks[key] = peaks[key]
    return sortedPeaks


def parseBedNucComposition(bedFile):
    """
    Reads tab-delimiter file with two rows describing the chromossomes and its lengths.
    Returns dictionary of chr:sizes.
    """
    prev_peak = ''
    peaksA = {}
    peaksC = {}
    peaksG = {}
    peaksT = {}
    with open(bedFile, "rb") as f:
        reader = csv.reader(f, delimiter="\t")
        for line, row in enumerate(reader):
            if str(row[0])[0] != '#': # skip header
                peak = str(row[3])
                bp = int(row[12])
                strand = str(row[5])
                total = float(13)
                A = int(row[16]) / total
                C = int(row[17]) / total
                G = int(row[18]) / total
                T = int(row[19]) / total
                if peak != prev_peak:
                    # new peak
                    peaksA[peak] = {}
                    peaksC[peak] = {}
                    peaksG[peak] = {}
                    peaksT[peak] = {}
                    if strand == "+":
                        peaksA[peak][bp - 2001] = A
                        peaksC[peak][bp - 2001] = C
                        peaksG[peak][bp - 2001] = G
                        peaksT[peak][bp - 2001] = T
                    elif strand == "-":
                        peaksA[peak][2001 - bp] = A
                        peaksC[peak][2001 - bp] = C
                        peaksG[peak][2001 - bp] = G
                        peaksT[peak][2001 - bp] = T
                    prev_peak = peak
                else:
                    if strand == "+":
                        peaksA[peak][bp - 2001] = A
                        peaksC[peak][bp - 2001] = C
                        peaksT[peak][bp - 2001] = G
                        peaksG[peak][bp - 2001] = T
                    elif strand == "-":
                        peaksA[peak][2001 - bp] = A
                        peaksC[peak][2001 - bp] = C
                        peaksG[peak][2001 - bp] = G
                        peaksT[peak][2001 - bp] = T
                    prev_peak = peak
    return (peaksA, peaksC, peaksG, peaksT)


def parseBedConservation(bedFile):
    """
    """
    prev_peak = ''
    peaks = {}
    with open(bedFile, "rb") as f:
        reader = csv.reader(f, delimiter="\t")
        for line, row in enumerate(reader):
            peak = str(row[3])
            bp = int(row[8])
            cons = float(row[14])
            strand = str(row[5])
            if peak != prev_peak:
                # new peak
                peaks[peak] = {}
                if strand == "+":
                    peaks[peak][bp - 2001] = cons
                elif strand == "-":
                    peaks[peak][2001 - bp] = cons
                prev_peak = peak
            else:
                if strand == "+":
                    peaks[peak][bp - 2001] = cons
                elif strand == "-":
                    peaks[peak][2001 - bp] = cons
                prev_peak = peak
    sortedPeaks = {}
    for key in sorted(peaks.iterkeys()):
        sortedPeaks[key] = peaks[key]
    return sortedPeaks


# Parse files into dataframes
covCM = parseBedCoverage(projectDir + 'bed/' + sample + '_tss_coverage.bed')
covCMDF = pd.DataFrame(covCM).T
covChIP = parseBedCoverage(projectDir + 'bed/' + sample + '_tss_coverage_ChIP.bed')
covChIPDF = pd.DataFrame(covChIP).T
covIgG = parseBedCoverage(projectDir + 'bed/' + sample + '_tss_coverage_IgG.bed')
covIgGDF = pd.DataFrame(covIgG).T
covDNase = parseBedCoverage(projectDir + 'bed/' + sample + '_tss_coverage_DNase.bed')
covDNaseDF = pd.DataFrame(covDNase).T

A,C,G,T = parseBedNucComposition(projectDir + 'bed/' + sample + '_tss_nucleotide_compos.bed')
nucCompA = pd.DataFrame(A).T
nucCompC = pd.DataFrame(C).T
nucCompG = pd.DataFrame(G).T
nucCompT = pd.DataFrame(T).T

cons = parseBedConservation(projectDir + 'bed/' + sample + '_tss_conservation.bed')
consDF = pd.DataFrame(cons).T

# Pickle them all!
#cons = parseBedConservation(projectDir + 'bed/' + sample + '_tss_conservation.bed')
#with open(projectDir + 'pickles/conservation_' + sample + '.pickle', 'wb') as handle:
#    pickle.dump(cons, handle)

# reading in
#with open(projectDir + 'pickles/coverage_' + sample + 'average.pickle', 'rb') as handle:
#    aveCov = pickle.load(handle)

# plot overal TSS 
plt.plot(covCMDF.columns, covCMDF.mean(), '#F57D05', label = "ChIPmentation")
plt.plot(covChIPDF.columns, covChIPDF.mean(), '#0591F5', label = "ChIP")
plt.plot(covIgGDF.columns, covIgGDF.mean(), 'k', label = "IgG")
plt.plot(covDNaseDF.columns, covDNaseDF.mean(), '#45B06C', label = "DNase")
plt.plot(consDF.columns, consDF.mean(), '#F50559', label = "Conservation")
plt.legend(loc='upper left')
pylab.xlim([-2000,2000])
pylab.ylim([0,1.5])
plt.savefig("TSSs_H3K4me3.pdf", format='pdf')

# focus on TATA
f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
ax1.plot(covCMDF.columns, covCMDF.mean(), '#F57D05', label = "ChIPmentation")
ax1.plot(covChIPDF.columns, covChIPDF.mean(), '#0591F5', label = "ChIP")
ax1.plot(covIgGDF.columns, covIgGDF.mean(), 'k', label = "IgG")
ax1.plot(consDF.columns, consDF.mean(), '#F50559', label = "Conservation")
ax1.set_ylim(0.01, 1.1)
ax1.vlines(0,0,10)
ax1.legend(loc='lower left',prop={'size':8})
ax2.plot(nucCompA.columns, nucCompA.mean(), 'g', label = "A")
ax2.plot(nucCompC.columns, nucCompC.mean(), 'b', label = "C")
ax2.plot(nucCompG.columns, nucCompG.mean(), 'y', label = "G")
ax2.plot(nucCompT.columns, nucCompT.mean(), 'r', label = "T")
ax2.set_ylim(0.15, 0.34)
ax2.set_xlim(-100, 50)
ax2.vlines(0,0,10)
ax2.legend(loc='upper left',prop={'size':8})
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

plt.savefig("TSSs_H3K4me3_spike.pdf", format='pdf')











samplepeak = 'chr11:62495524:62495620:+:0.54:0.999268'


for peak in cons:
    if len(cons[peak]) == 4000:
        print(max([cons[peak][x] for x in range(1,4000)]))


