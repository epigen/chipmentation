from __future__ import print_function
import sys
import csv
import pandas as pd

sample = str(sys.argv[1])

def parseBedCoverage(bedFile, peakCol, bpCol, countCol):
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


cov = parseBedCoverage(sample, 3, 4, 5)

# add to pandas df
covDF = pd.DataFrame(cov).T

# write out
sampleOut = '.'.join(sample.split(".")[:-1]) + ".csv"
covDF.to_csv(sampleOut)

print("Success")

