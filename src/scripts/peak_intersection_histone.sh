#!/bin/bash

PROJECTDIR=/home/arendeiro/data/human/chipmentation
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersections.txt

for SAMPLE_NAME in H3K4me3_K562_500k_CM
    do
    SAMPLE=${SAMPLE_NAME/_CM/}
    SAMPLE_CM500k=${SAMPLE}_CM
    SAMPLE_CM10k=${SAMPLE_CM500k/500k/10k}
    SAMPLE_CHIP=${SAMPLE}_ChIP
    SAMPLE_ENCODE=${SAMPLE}_ChIP_Encode

    # Merge peaks within samples that are <100bp apart
    echo "Merging peaks in" $SAMPLE_CM500k
    bedtools merge -d 100 \
    -i $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak \
    > $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.merged.bed

    echo "Merging peaks in" $SAMPLE_CM10k
    bedtools merge -d 100 \
    -i $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak \
    > $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.merged.bed

    echo "Merging peaks in" $SAMPLE_CHIP
    bedtools merge -d 100 \
    -i $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    > $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed

    echo "Merging peaks in" $SAMPLE_ENCODE
    bedtools merge -d 100 \
    -i $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak \
    > $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.merged.bed

    # Count total number of peaks
    echo "Counting peaks in" $SAMPLE_CM500k
    SAMPLE_CM500k_COUNT="$(wc -l < $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.merged.bed)"
    echo "Counting peaks in" $SAMPLE_CM10k
    SAMPLE_CM10k_COUNT="$(wc -l < $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.merged.bed)"
    echo "Counting peaks in" $SAMPLE_CHIP
    SAMPLE_CHIP_COUNT="$(wc -l < $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed)"
    echo "Counting peaks in" $SAMPLE_ENCODE
    SAMPLE_ENC_COUNT="$(wc -l < $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.merged.bed)"

    # Intersect and count
    echo "Intersecting" $SAMPLE_CM500k "with" $SAMPLE_CM10k
    CM500k_CM10k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_CM500k "with" $SAMPLE_CHIP
    CM500k_CHIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_CM500k "with" $SAMPLE_ENCODE
    CM500k_ENC="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_CM10k "with" $SAMPLE_CM500k
    CM10k_CM="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_CM10k "with" $SAMPLE_CHIP
    CM10k_CHIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_CM10k "with" $SAMPLE_ENCODE
    CM10k_ENC="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_CM500k
    CHIP_CM="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_CM10k
    CHIP_CM10k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_ENCODE
    CHIP_ENC="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_ENCODE "with" $SAMPLE_ENCODE
    ENC_CM500k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_ENCODE "with" $SAMPLE_CM10k
    ENC_CM10k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.merged.bed \
    | wc -l)"

    echo "Intersecting" $SAMPLE_ENCODE "with" $SAMPLE_CHIP
    ENC_ChIP500k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.merged.bed \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.merged.bed \
    | wc -l)"

    echo "sample" "CM500k_peak_count" "CM10k_peak_count" "ChIP_peak_count" "Encode_peak_count" "CM500k_CM10k_intersection" "CM_ChIP_intersection" "CM_Encode_intersection" "CM10k_CM_intersection" "CM10k_ChIP_intersection" "CM10k_Encode_intersection" "ChIP_CM_intersection" "ChIP_CM10k_intersection" "ChIP_Encode_intersection" "Encode_CM_intersection" "Encode_CM10k_intersection" "Encode_ChIP_intersection" >> $INTERSECTIONS
    echo $SAMPLE_CM500k $SAMPLE_CM500k_COUNT $SAMPLE_CM10k_COUNT $SAMPLE_CHIP_COUNT $SAMPLE_ENC_COUNT $CM500k_CM10k $CM500k_CHIP $CM500k_ENC $CM10k_CM500k $CM10k_CHIP $CM10k_ENC $CHIP_CM500k $CHIP_CM10k $CHIP_ENC $ENC_CM500k $ENC_CM10k $ENC_CHIP >> $INTERSECTIONS

done

column $INTERSECTIONS > tmp
mv tmp $INTERSECTIONS

date
