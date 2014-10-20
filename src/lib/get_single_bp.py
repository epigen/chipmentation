#!/usr/bin/env python
import csv, sys

wr = csv.writer(sys.stdout, delimiter = '\t', lineterminator='\n')

for row in csv.reader(iter(sys.stdin.readline, ''), delimiter = '\t'):
    row[1] = int(row[1]) + int(row[4])
    row[2] = int(row[1]) + 1
    wr.writerow(row)

