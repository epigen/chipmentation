# Call peaks
for SAMPLE in wgEncodeUwDnaseK562Aln.merged K562_50K_ATAC_nan_nan_nan_0_0_hg19.trimmed.bowtie2.shifted.dups
do
    NAME=`basename $SAMPLE`
    macs2 callpeak \
    --broad -n $NAME
    -t cemm-backup/chipmentation/data/mapped/${SAMPLE}.bam
    --outdir cemm-backup/chipmentation/data/peaks/$SAMPLE
done

# Prepare mappability files
cd cemm-backup/reference/hg19/

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig

bigWigToBedgraph wgEncodeCrgMapabilityAlign36mer.bigWig wgEncodeCrgMapabilityAlign36mer.bedgraph
bigWigToBedgraph wgEncodeCrgMapabilityAlign50mer.bigWig wgEncodeCrgMapabilityAlign50mer.bedgraph

awk '$4 == 1 {OFS="\t"; print $1, $2, $3}' wgEncodeCrgMapabilityAlign36mer.bedgraph > wgEncodeCrgMapabilityAlign36mer.bed
awk '$4 == 1 {OFS="\t"; print $1, $2, $3}' wgEncodeCrgMapabilityAlign50mer.bedgraph > wgEncodeCrgMapabilityAlign50mer.bed

# Filter for mappability
cd
cd cemm-backup/chipmentation/data/peaks

bedtools intersect -wa -u -f 1 \
-a wgEncodeUwDnaseK562Aln.merged/wgEncodeUwDnaseK562Aln.merged_peaks.broadPeak \
-b ~/cemm-backup/reference/hg19/wgEncodeCrgMapabilityAlign36mer.bed \
> wgEncodeUwDnaseK562Aln.merged/wgEncodeUwDnaseK562Aln.merged_peaks.filtered.broadPeak

bedtools intersect -wa -u -f 1 \
-a K562_50K_ATAC_nan_nan_nan_0_0_hg19/K562_50K_ATAC_nan_nan_nan_0_0_hg19_peaks.broadPeak \
-b ~/cemm-backup/reference/hg19/wgEncodeCrgMapabilityAlign50mer.bed \
> K562_50K_ATAC_nan_nan_nan_0_0_hg19/K562_50K_ATAC_nan_nan_nan_0_0_hg19_peaks.filtered.broadPeak


# Center peaks on motifs
for SAMPLE in K562_50K_ATAC_nan_nan_nan_0_0_hg19 wgEncodeUwDnaseK562Aln.merged
do
    for TF in PU1 GATA1 REST CTCF
    do
        annotatePeaks.pl ${SAMPLE}/${SAMPLE}_peaks.filtered.broadPeak \
        hg19 \
        -size 2000 -center ../motifs/K562_10M_CM_${TF}_nan_nan_0_0_hg19/homerResults/motif1.motif | \
        awk -v OFS='\t' '{print $2, $3, $4, $1, $6, $5}' | \
        awk -v OFS='\t' -F '\t' '{ gsub("0", "+", $6) ; gsub("1", "-", $6) ; print }' | \
        python ~/fix_bedfile_genome_boundaries.py hg19 | \
        sortBed > ${SAMPLE}/${SAMPLE}_peaks.filtered.${TF}-motifCentered.bed

        annotatePeaks.pl ${SAMPLE}/${SAMPLE}_peaks.filtered.${TF}-motifCentered.bed \
        hg19 \
        -mask -mscore -m ../motifs/K562_10M_CM_${TF}_nan_nan_0_0_hg19/homerResults/motif1.motif | \
        tail -n +2 | cut -f 1,5,22 \
        > ${SAMPLE}/${SAMPLE}_peaks.filtered.${TF}-motifAnnotated.bed
        # cut this ^^
    done
done


# 2nd approach: start from ChIP-seq peaks
# Center peaks on motifs
for TF in CTCF PU1 GATA1 REST
do
    annotatePeaks.pl K562_10M_CM_${TF}_nan_nan_0_0_hg19/K562_10M_CM_${TF}_nan_nan_0_0_hg19_peaks.narrowPeak \
    hg19 \
    -size 2000 -center ../motifs/K562_10M_CM_${TF}_nan_nan_0_0_hg19/homerResults/motif1.motif | \
    awk -v OFS='\t' '{print $2, $3, $4, $1, $6, $5}' | \
    awk -v OFS='\t' -F '\t' '{ gsub("0", "+", $6) ; gsub("1", "-", $6) ; print }' | \
    python ~/fix_bedfile_genome_boundaries.py hg19 | \
    sortBed > K562_10M_CM_${TF}_nan_nan_0_0_hg19/K562_10M_CM_${TF}_nan_nan_0_0_hg19_peaks.motifCentered.bed

    annotatePeaks.pl K562_10M_CM_${TF}_nan_nan_0_0_hg19/K562_10M_CM_${TF}_nan_nan_0_0_hg19_peaks.motifCentered.bed \
    hg19 \
    -mask -mscore -m ../motifs/K562_10M_CM_${TF}_nan_nan_0_0_hg19/homerResults/motif1.motif | \
    tail -n +2 | cut -f 1,5,22 \
    > K562_10M_CM_${TF}_nan_nan_0_0_hg19/K562_10M_CM_${TF}_nan_nan_0_0_hg19_peaks.motifAnnotated.bed
    # cut this ^^
done
