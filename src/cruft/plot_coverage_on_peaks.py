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
#sample = 'CTCF_K562_10mio_CM'
sample = 'PU1_K562_10mio_CM'
projectDir = '/home/afr/Documents/workspace/chipmentation/data/'

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
            bp = int(row[4])
            count = int(row[5])
            if peak != prev_peak:
                # new peak
                peaks[peak] = {}
                peaks[peak][bp] = count
                prev_peak = peak
            else:
                peaks[peak][bp] = count

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
    peaks = {}
    with open(bedFile, "rb") as f:
        reader = csv.reader(f, delimiter="\t")
        for line, row in enumerate(reader):
            if str(row[0])[0] != '#':
                # skip header
                peak = str(row[3])
                A = int(row[6])
                C = int(row[7])
                G = int(row[8])
                T = int(row[9])
                if peak != prev_peak:
                    # new peak
                    bp = 1
                    peaks[peak] = {}
                    peaks[peak][bp] = [A, C, G, T]
                    prev_peak = peak
                else:
                    bp += 1
                    peaks[peak][bp] = [A, C, G, T]

    sortedPeaks = {}
    for key in sorted(peaks.iterkeys()):
        sortedPeaks[key] = peaks[key]
    return sortedPeaks


def parseBedConservation(bedFile):
    """
    """
    prev_peak = ''
    peaks = {}
    with open(bedFile, "rb") as f:
        reader = csv.reader(f, delimiter="\t")
        for line, row in enumerate(reader):

            # skip header
            peak = str(row[3])
            bp = int(row[4])
            cons = float(row[10])
            if peak != prev_peak:
                # new peak
                peaks[peak] = {}
                peaks[peak][bp] = cons
                prev_peak = peak
            else:
                peaks[peak][bp] = cons
                
    sortedPeaks = {}
    for key in sorted(peaks.iterkeys()):
        sortedPeaks[key] = peaks[key]
    return sortedPeaks

# Parse files into dataframes
covCM = parseBedCoverage(projectDir + 'bed/' + sample + '_peak_coverage.bed')
covCMDF = pd.DataFrame(covCM).T
covChIP = parseBedCoverage(projectDir + 'bed/' + sample + '_peak_coverage_ChIP.bed')
covChIPDF = pd.DataFrame(covChIP).T
covIgG = parseBedCoverage(projectDir + 'bed/' + sample + '_peak_coverage_IgG.bed')
covIgGDF = pd.DataFrame(covIgG).T
covDNase = parseBedCoverage(projectDir + 'bed/' + sample + '_peak_coverage_DNase.bed')
covDNaseDF = pd.DataFrame(covDNase).T

nucComp = parseBedNucComposition(projectDir + 'bed/' + sample + '_peak_nucleotide_compos.bed')
nucCompDF = pd.DataFrame(nucComp).T

cons = parseBedConservation(projectDir + 'bed/' + sample + '_peak_conservation.bed')
consDF = pd.DataFrame(cons).T

# Pickle them all!
#cons = parseBedConservation(projectDir + 'bed/' + sample + '_peak_conservation.bed')
#with open(projectDir + 'pickles/conservation_' + sample + '.pickle', 'wb') as handle:
#    pickle.dump(cons, handle)

# reading in
#with open(projectDir + 'pickles/coverage_' + sample + 'average.pickle', 'rb') as handle:
#    aveCov = pickle.load(handle)


# plot
plt.plot(covCMDF.columns, covCMDF.mean(), label = "ChIPmentation")
plt.plot(covChIPDF.columns, covChIPDF.mean(), label = "ChIP")
plt.plot(covIgGDF.columns, covIgGDF.mean(), label = "IgG")
plt.plot(covDNaseDF.columns, covDNaseDF.mean(), label = "DNase")
plt.plot(nucCompDF.columns, nucCompDF.mean(), label = "Nucleotide composition")
plt.plot(consDF.columns, consDF.mean(), label = "Conservation")
plt.legend(loc='upper left')

plt.savefig("/home/arendeiro/projects/chipmentation/results/plots/peaks_PU1.png", dpi = 600)



# focus on footprint
f, (ax1, ax2) = plt.subplots(2, sharex=False, sharey=False)
ax1.plot(range(-2000,2000), covCMDF.mean(), '#F57D05', label = "ChIPmentation")
ax1.plot(range(-2000,2000), covChIPDF.mean(), '#0591F5', label = "ChIP")
ax1.plot(range(-2000,2000), covDNaseDF.mean(), '#45B06C', label = "DNase")
#ax1.plot(range(-2000,2000), covIgGDF.mean(), 'k', label = "IgG")
ax1.set_ylim(0.01, 0.24)
ax1.set_xlim(-500, 500)
ax1.legend(loc='lower left',prop={'size':8})
ax2.plot(range(-2000,2000), covCMDF.mean(), '#F57D05', label = "ChIPmentation")
ax2.plot(range(-2000,2000), covChIPDF.mean(), '#0591F5', label = "ChIP")
ax2.plot(range(-2000,2000), covDNaseDF.mean(), '#45B06C', label = "DNase")
#ax2.plot(range(-2000,2000), covIgGDF.mean(), 'k', label = "IgG")
#ax2.plot(range(-2000,2000), consDF.mean(), '#F50559', label = "Conservation")
ax2.set_ylim(0.01, 0.24)
ax2.set_xlim(-100, 100)
ax2.legend(loc='upper left',prop={'size':8})
# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

plt.savefig("TSSs_H3K4me3_spike.pdf", format='pdf')

