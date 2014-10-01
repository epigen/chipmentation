#!/bin/bash

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt

while read SAMPLE_NAME CONTROL_NAME; do
    sbatch /home/arendeiro/projects/chipmentation/src/optimize_peak_calling_macs2.sh \
    $SAMPLE_NAME \
    $CONTROL_NAME \
    --job-name="peakcalling_${SAMPLE_NAME}" \
    --output=/home/arendeiro/logs/mapping_${SAMPLE_NAME}.out
done < $SAMPLES_FILE

# more than fold 10 - dont do
#while read SAMPLE_NAME CONTROL_NAME; do
#    awk '$9 >= 10' $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200/${SAMPLE_NAME}_bw200_peaks.narrowPeak \
#    > $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200/${SAMPLE_NAME}_bw200_peaks_Fold10.narrowPeak
#done < $SAMPLES_FILE

# merge peaks within 200bps between replicates
while read SAMPLE_NAME CONTROL_NAME; do
    bedtools merge -d 200 -i $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200/${SAMPLE_NAME}_bw200_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak
done < $SAMPLES_FILE


BEDTOOLSDIR=/home/arendeiro/.local/software/bedtools2/bin


COUNT=0
while read SAMPLE_NAME CONTROL_NAME; do
    SAMPLE=H3K4me3_K562_500k_ChIP
    SAMPLE1=H3K4me3_K562_500k_ChIP_CM11-9_R1
    SAMPLE2=H3K4me3_K562_500k_ChIP_CM12-9_R2

    bedtools merge -d 100 -i $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_MACS2_bw200/${SAMPLE_NAME}_bw200_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak
    
    if [[ $((COUNTER%2)) -eq 0 ]]
        then
        echo "Replicate 1" $SAMPLE_NAME
    else
        echo "Replicate 2" $SAMPLE_NAME
        NAME2=${SAMPLE_NAME/R1/R2}
        cat $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak $PROJECTDIR/peaks/${SAMPLE2}.narrowPeak | \
        bedtools merge -i stdin | \
        bedtools intersect -wa -a stdin -b $PROJECTDIR/peaks/${SAMPLE_NAME}.narrowPeak > $PROJECTDIR/peaks/${NAME}_peaks.bed
    fi
    COUNTER=$[COUNTER + 1]
done < $SAMPLES_FILE


SAMPLE=H3K4me3_K562_500k_ChIP
SAMPLE1=H3K4me3_K562_500k_ChIP_CM11-9_R1
SAMPLE2=H3K4me3_K562_500k_ChIP_CM12-9_R2
bedtools merge -d 100 -i $PROJECTDIR/peaks/${SAMPLE2}_peaks_MACS2_bw200/${SAMPLE2}_bw200_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE2}.narrowPeak
cat $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak $PROJECTDIR/peaks/${SAMPLE2}.narrowPeak | \
bedtools sort -i stdin | \
bedtools merge -i stdin | \
bedtools intersect -wa -a stdin -b $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak > $PROJECTDIR/peaks/${SAMPLE}_peaks.bed

SAMPLE=H3K4me3_K562_500k_CM
SAMPLE1=H3K4me3_K562_500k_CM_CM11-1_R1
SAMPLE2=H3K4me3_K562_500k_CM_CM12-1_R2
bedtools merge -d 100 -i $PROJECTDIR/peaks/${SAMPLE2}_peaks_MACS2_bw200/${SAMPLE2}_bw200_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE2}.narrowPeak
cat $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak $PROJECTDIR/peaks/${SAMPLE2}.narrowPeak | \
bedtools sort -i stdin | \
bedtools merge -i stdin | \
bedtools intersect -wa -a stdin -b $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak > $PROJECTDIR/peaks/${SAMPLE}_peaks.bed

SAMPLE=H3K4me3_K562_10k_CM
SAMPLE1=H3K4me3_K562_10k_CM_CM11-5_R1
SAMPLE2=H3K4me3_K562_10k_CM_CM12-5_R2
bedtools merge -d 100 -i $PROJECTDIR/peaks/${SAMPLE2}_peaks_MACS2_bw200/${SAMPLE2}_bw200_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE2}.narrowPeak
cat $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak $PROJECTDIR/peaks/${SAMPLE2}.narrowPeak | \
bedtools sort -i stdin | \
bedtools merge -i stdin | \
bedtools intersect -wa -a stdin -b $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak > $PROJECTDIR/peaks/${SAMPLE}_peaks.bed


SAMPLE=H3K4me3_K562_500k_ChIP
SAMPLE1=H3K4me3_K562_500k_ChIP_CM11-9_R1
SAMPLE2=H3K4me3_K562_500k_ChIP_CM12-9_R2
bedtools merge -d 100 -i $PROJECTDIR/peaks/${SAMPLE1}_peaks_MACS2_bw200/${SAMPLE1}_bw200_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak
bedtools merge -d 100 -i $PROJECTDIR/peaks/${SAMPLE2}_peaks_MACS2_bw200/${SAMPLE2}_bw200_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE2}.narrowPeak
cat $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak $PROJECTDIR/peaks/${SAMPLE2}.narrowPeak | \
bedtools sort -i stdin | \
bedtools merge -i stdin | \
bedtools intersect -wa -a stdin -b $PROJECTDIR/peaks/${SAMPLE1}.narrowPeak > $PROJECTDIR/peaks/${SAMPLE}_peaks.bed


wc -l $PROJECTDIR/peaks/H3K4me3_K562_500k_ChIP_peaks.bed
wc -l $PROJECTDIR/peaks/H3K4me3_K562_500k_CM_peaks.bed
wc -l $PROJECTDIR/peaks/H3K4me3_K562_10k_CM_peaks.bed
bedtools intersect -u -wa -a $PROJECTDIR/peaks/H3K4me3_K562_500k_ChIP_peaks.bed -b $PROJECTDIR/peaks/H3K4me3_K562_500k_CM_peaks.bed | wc -l
bedtools intersect -u -wa -a $PROJECTDIR/peaks/H3K4me3_K562_500k_CM_peaks.bed -b $PROJECTDIR/peaks/H3K4me3_K562_10k_CM_peaks.bed | wc -l
bedtools intersect -u -wa -a $PROJECTDIR/peaks/H3K4me3_K562_10k_CM_peaks.bed -b $PROJECTDIR/peaks/H3K4me3_K562_500k_CM_peaks.bed | wc -l