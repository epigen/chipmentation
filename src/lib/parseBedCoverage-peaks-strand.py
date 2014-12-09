from __future__ import print_function
import sys
import csv
import pandas as pd

sample = str(sys.argv[1])

def parseBedCoverage(bedFile, peakCol, strandCol, bpCol, countCol):
    """
    Reads tab-delimiter file with two rows describing the chromossomes and its lengths.
    Returns dictionary of chr:sizes.
    """
    prev_peak = ''
    peaks = {}
    with open(bedFile, "rb") as f:
        reader = csv.reader(f, delimiter="\t")
        for line, row in enumerate(reader):
            peak = str(row[peakCol])
            bp = int(row[bpCol])
            count = int(row[countCol])
            strand = int(row[strandCol])
            if strand == 0:
                peak_name = peak + "_+"
            elif strand == 1:
                peak_name = peak + "_-"
            if peak != prev_peak:
                # new peak                
                peaks[peak_name] = {}
                peaks[peak_name][bp] = count
                prev_peak = peak
            else:
                peaks[peak_name][bp] = count
    sortedPeaks = {}
    for key in sorted(peaks.iterkeys()):
        sortedPeaks[key] = peaks[key]
    return sortedPeaks


cov = parseBedCoverage(sample, 3, 5, 6, 7)

# add to pandas df
covDF = pd.DataFrame(cov).T

# invert negative strand peaks
pos = pd.DataFrame([covDF.ix[row] for row in range(len(covDF)) if "+" in covDF.index[row]])
neg = pd.DataFrame([covDF.ix[row][::-1] for row in range(len(covDF)) if "-" in covDF.index[row]])
neg.columns = pos.columns
covDFRev = pos.append(neg)
covDFRev.columns = ["X" + str(i) for i in covDFRev.columns]

# write out
sampleOut = '.'.join(sample.split(".")[:-1]) + ".csv"
covDFRev.to_csv(sampleOut)
print("Success")

