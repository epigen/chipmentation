#!/usr/bin/env python
import csv, sys

wr = csv.writer(sys.stdout, delimiter = '\t', lineterminator='\n')

for row in csv.reader(iter(sys.stdin.readline, ''), delimiter = '\t'):
    if row[0][0] == "@":
        wr.writerow(row)
    else:
        row[5] = "1M"
        wr.writerow(row)