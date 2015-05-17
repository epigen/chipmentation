# Nucleosome coverage
DYADS=/media/afr/cemm-backup1/chipmentation/data/nucleoATAC/10M_CM_H3K4ME1_PE/dyads_3col.bed
GMDYADS=/home/afr/GM.nucpos.bed
CTCF=/media/afr/cemm-backup1/chipmentation/data/peaks/K562_10M_CHIP_CTCF_nan_nan_0_0_hg19/K562_10M_CHIP_CTCF_nan_nan_0_0_hg19_peaks.motifCentered.bed
CAGE=/media/afr/cemm-backup1/chipmentation/annotation/hg19.cage_peak.robust.TATA_Annotated.expr_Annotated.K562_Expressed.1kb.tsv
REFSEQ=/media/afr/cemm-backup1/reference/hg19/GRCh37_hg19_refSeq.tss.bed

bedtools slop -b 500 \
-g /media/afr/cemm-backup1/reference/hg19/hg19.chrom-sizes.tsv \
-i $REFSEQ > refseq.1kb.bed

bedtools coverage -d -a refseq.1kb.bed -b $DYADS > tss.CMdyad.coverage.bed
bedtools coverage -d -a refseq.1kb.bed -b $GMDYADS > tss.GMdyad.coverage.bed
cut -f 6,12,13 tss.CMdyad.coverage.bed > tss.CMdyad.coverage.3col.bed
cut -f 6,12,13 tss.GMdyad.coverage.bed > tss.GMdyad.coverage.3col.bed


ipython
import pandas as pd
import numpy as np
cm = pd.read_csv('tss.CMdyad.coverage.3col.bed', sep="\t", header=None)
cm.columns = ['strand', 'bp', 'count']

cm.loc[cm['strand'] == "-", 'bp'] = abs(cm.loc[cm['strand'] == "-", 'bp'] - 1001)
gm = gm.groupby(['bp']).apply(np.mean)
plt.plot(s['count'])

gm = pd.read_csv('tss.GMdyad.coverage.3col.bed', sep="\t", header=None)
gm.columns = ['strand', 'bp', 'count']

gm.loc[gm['strand'] == "-", 'bp'] = abs(gm.loc[gm['strand'] == "-", 'bp'] - 1001)
gm = gm.groupby(['bp']).apply(np.mean)
plt.plot(s['count'])


# Fraction of overlap with MNase-seq nucleosomes
DYADS=/media/afr/cemm-backup1/chipmentation/data/nucleoATAC/10M_CM_H3K4ME3_PE/dyads_3col.bed
NUCLEO=/media/afr/cemm-backup1/chipmentation/data/nucleoATAC/10M_CM_H3K4ME3_PE/nucleosomes.bed
MNASENUCLEO=/media/afr/cemm-backup1/dyads/K562_nucleosomes.bed

awk '{ if ($5 >= 0.5) print }' $NUCLEO > filtered.bed

bedtools intersect -u -f 0.5 -r \
-a filtered.bed \
-b $MNASENUCLEO \
| wc -l

wc -l filtered.bed
