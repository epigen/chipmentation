#!/usr/bin/env python
import csv, sys

wr = csv.writer(sys.stdout, delimiter = '\t', lineterminator='\n')

for row in csv.reader(iter(sys.stdin.readline, ''), delimiter = '\t'):
    if str(row[0]) != "chrm":
        chrm =  str(row[0])
        start = int(row[1])
        end = int(row[2])
        strand = str(row[3])
        peak = str(row[4])
        cluster = str(row[5])
        
        l = end - start
        for bp in range(l):
            row[1] = start + bp - 3
            row[2] = start + bp + 3
            wr.writerow(row)
