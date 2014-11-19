#!/usr/bin/env python
import csv, sys

wr = csv.writer(sys.stdout, delimiter = '\t', lineterminator='\n')

for row in csv.reader(iter(sys.stdin.readline, ''), delimiter = '\t'):
    if row[0] != '00Annotation':
        chrm = row[0].split(":")[0]
        tmp = row[0].split(":")[1]
        strand = tmp.split(",")[1]
        tmp2 = tmp.split(",")[0]
        start = tmp2.split("..")[0]
        end = tmp2.split("..")[1]
        peak = row[0]
        wr.writerow([chrm, start, end, peak, strand] + row[1:])
