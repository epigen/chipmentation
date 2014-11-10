#!/usr/bin/env python
import csv, sys

wr = csv.writer(sys.stdout, delimiter = '\t', lineterminator='\n')

chrSizes = "/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt"

def getChrSizes(chrmFile):
    """
    Reads tab-delimiter file with two rows describing the chromossomes and its lengths.
    Returns dictionary of chr:sizes.
    """
    with open(chrmFile, 'r') as f:
        chrmSizes = {}
        for line in enumerate(f):
            row = line[1].strip().split('\t')
            chrmSizes[str(row[0])] = int(row[1])
    return chrmSizes
    
chrms = getChrSizes(chrSizes)

for row in csv.reader(iter(sys.stdin.readline, ''), delimiter = '\t'):
    chrm = row[0]
    start = int(row[1])
    end = int(row[2])

    if start >= 1 and end <= chrms[chrm] and start < end:
        wr.writerow(row)
    else:
        continue
    