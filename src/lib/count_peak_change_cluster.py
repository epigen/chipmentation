#!/usr/bin/env python
import csv, sys
import pandas as pd

file1 = sys.argv[1]
file2 = sys.argv[2]

#dataDir = '/fhgfs/groups/lab_bock/shared/data/cage_tss/1000it/'
#file1 = "hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage_K_G5.1.kgg"
#file2 = "hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage_K_G5.2.kgg"

p1 = pd.io.parsers.read_csv(file1, sep = "\t")
p1.columns = ["peak", "c1"]
p2 = pd.io.parsers.read_csv(file2, sep = "\t")
p2.columns = ["peak", "c2"]

p = p1.merge(p2, how = "left")

#p.boxplot(by='c2')
score = []
for k in range(5):
    tmp = p.ix[p.c1 == k]
    score.append([])
    for kk in range(5):
        t = tmp.ix[tmp.c2 == kk].shape[0]
        score[k].append(t)
print(score)
