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
    if row[0][0] == "@":
        wr.writerow(row)
    else:
        strand = int(row[1])
        chrm = row[2]

        if strand == 0:
            if int(row[3]) + 4 + 51 < chrms[chrm]:
                row[3] = int(row[3]) + 4
            else:
                continue
        elif strand == 16:
            if int(row[3]) - 5 - 51 >= 1:
                row[3] = int(row[3]) - 5
            else:
                continue
        wr.writerow(row)