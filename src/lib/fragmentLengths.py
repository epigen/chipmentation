#!/usr/bin/env python

"""
Reads from stdin a SAM file, creates Counter with fragment length distribution,
saves pickle with distribution Counter using a provided prefix
"""

import csv
import sys
from collections import Counter
import pickle

lens = Counter()

name = sys.argv[1]

i = 0
for row in csv.reader(iter(sys.stdin.readline, ''), delimiter='\t'):
    if row[0][0] == "@":
        continue
    else:
        lens[abs(int(row[8]))] += 1
    i += 1

pickle.dump(lens, open(name + ".pickle", "wb"), protocol=pickle.HIGHEST_PROTOCOL)
