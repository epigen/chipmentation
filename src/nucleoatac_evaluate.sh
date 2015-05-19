# Post process NucleoATAC output

CHROMSIZES=$HOME/reference/Homo_sapiens/hg19.chrom-sizes.tsv
NUCLEOATACDIR=$HOME/chipmentation/data/nucleoATAC
CTCF=$HOME/chipmentation/data/peaks/K562_10M_CHIP_CTCF_nan_nan_0_0_hg19/K562_10M_CHIP_CTCF_nan_nan_0_0_hg19_peaks.motifCentered.bed
REFSEQ=$HOME/reference/Homo_sapiens/GRCh37_hg19_refSeq.tss.bed

bedtools slop -b 500 \
-g $CHROMSIZES \
-i $REFSEQ > $HOME/reference/Homo_sapiens/GRCh37_hg19_refSeq.tss.1kb.bed

REFSEQ=$HOME/reference/Homo_sapiens/GRCh37_hg19_refSeq.tss.1kb.bed

for SAMPLE in K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19 # K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19 K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19-K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19_peaks
do
    echo $SAMPLE

    ## Filter very stringently
    echo "Filtering " ${SAMPLE}
    awk '{if($5 > 0.8) print}' $NUCLEOATACDIR/${SAMPLE}.nucpos.bed > $NUCLEOATACDIR/${SAMPLE}.nucpos-filtered.bed

    ## Get nucleosomes
    echo "Getting nucleosomes " ${SAMPLE}
    bedtools slop -b 74 -g $CHROMSIZES -i $NUCLEOATACDIR/${SAMPLE}.nucpos.bed > $NUCLEOATACDIR/${SAMPLE}.nucleosomes.bed
    bedtools slop -b 74 -g $CHROMSIZES -i $NUCLEOATACDIR/${SAMPLE}.nucpos-filtered.bed > $NUCLEOATACDIR/${SAMPLE}.nucleosomes-filtered.bed
   
    for TYPE in nucleosomes-filtered
    do
        ## Get coverage around TSSs
        echo "Getting TSS coverage for " ${SAMPLE}.$TYPE
        bedtools coverage -d -a $NUCLEOATACDIR/${SAMPLE}.$TYPE.bed -b $REFSEQ > $NUCLEOATACDIR/${SAMPLE}.$TYPE.TSScoverage.bed
        cut -f 6,7,8 $NUCLEOATACDIR/${SAMPLE}.$TYPE.TSScoverage.bed > t
        mv t $NUCLEOATACDIR/${SAMPLE}.$TYPE.TSScoverage.bed

        ## Get coverage around CTCF
        echo "Getting CTCF coverage for " ${SAMPLE}.$TYPE
        bedtools coverage -d -a $NUCLEOATACDIR/${SAMPLE}.$TYPE.bed -b $CTCF > $NUCLEOATACDIR/${SAMPLE}.$TYPE.CTCFcoverage.bed
        cut -f 6,7,8 $NUCLEOATACDIR/${SAMPLE}.$TYPE.CTCFcoverage.bed > t
        mv t $NUCLEOATACDIR/${SAMPLE}.$TYPE.CTCFcoverage.bed
    done
done


ipython --matplotlib=qt

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

nucleoAtacDir = "/home/arendeiro/chipmentation/data/nucleoATAC"

sample = "GM"
sample = "K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19"

dyadsTSS = pd.read_csv(os.path.join(nucleoAtacDir, sample + '.nucpos.bed'), sep="\t", header=None)
dyadsTSS.columns = ['strand', 'bp', 'count']

dyadsTSS.loc[dyadsTSS['strand'] == "-", 'bp'] = abs(dyadsTSS.loc[dyadsTSS['strand'] == "-", 'bp'] - 1001)

GMs = dyadsTSS.groupby(['bp'])['count'].apply(np.mean)
plt.plot(GMs)

