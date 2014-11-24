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


if sample == "/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.120bpSlop.CMcoverage.bed":
    cov = parseBedCoverage(sample, 3, 9, 10)
elif "TATA_Annotated.TATA" in sample or "TATA_Annotated.CpG" in sample:
    cov = parseBedCoverage(sample, 4, 11, 12)
elif "K562_expressed" in sample:
    cov = parseBedCoverage(sample, 4, 12, 13)
else:
    sys.exit(1)
    print("FAILED!!", file=sys.stderr)

# add to pandas df
covDF = pd.DataFrame(cov).T

# invert negative strand peaks (only for TSSs!)
pos = pd.DataFrame([covDF.ix[row] for row in range(len(covDF)) if "+" in covDF.index[row]])
neg = pd.DataFrame([covDF.ix[row][::-1] for row in range(len(covDF)) if "-" in covDF.index[row]])
neg.columns = pos.columns
covDFRev = pos.append(neg)
covDFRev.columns = ["X" + str(i) for i in covDFRev.columns]

# write out
sampleOut = '.'.join(sample.split(".")[:-1]) + ".csv"
covDFRev.to_csv(sampleOut)

print("Success")

