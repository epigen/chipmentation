#!/usr/bin/env python
import csv, sys


output = []
for row in csv.reader(iter(sys.stdin.readline, ''), delimiter = '\t'):
    if row[0][0] == "@":
        output.append(row)
    else:
        strand = int(row[1])

        if strand == 0:
            row[3] = int(row[3]) + 4
        elif strand == 16:
            row[3] = int(row[3]) - 5
        output.append(row)

wr = csv.writer(sys.stdout, delimiter = '\t', lineterminator='\n')
for line in output:
    wr.writerow(line)