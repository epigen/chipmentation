#!/bin/bash

PROJECTDIR=/home/arendeiro/data/human/chipmentation
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersections.txt
rm $INTERSECTIONS


for SAMPLE_NAME in H3K4me3_K562_500k_CM
    do
    SAMPLE=${SAMPLE_NAME/_CM/}
    SAMPLE_CM500k=${SAMPLE}_CM
    SAMPLE_CM10k=${SAMPLE_CM500k/500k/10k}
    SAMPLE_CHIP=${SAMPLE}_ChIP
    SAMPLE_ENCODE=${SAMPLE}_ChIP_Encode

    # Intersect and count
    echo "Intersecting" $SAMPLE_CM500k "with" $SAMPLE_CM10k
    CM500k_CM10k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak \
    | wc -l)"
    
    echo "$SAMPLE_CM500k $SAMPLE_CM10k" \
    $CM500k_CM10k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM500k "with" $SAMPLE_CHIP
    CM500k_CHIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM500k $SAMPLE_CHIP" \
    $CM500k_CHIP \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM500k "with" $SAMPLE_ENCODE
    CM500k_ENC="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM500k $SAMPLE_ENCODE" \
    $CM500k_ENC \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM10k "with" $SAMPLE_CM500k
    CM10k_CM500k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM10k $SAMPLE_CM500k" \
    $CM10k_CM500k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM10k "with" $SAMPLE_CHIP
    CM10k_CHIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM10k $SAMPLE_CHIP" \
    $CM10k_CHIP \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM10k "with" $SAMPLE_ENCODE
    CM10k_ENC="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM10k $SAMPLE_ENCODE" \
    $CM10k_ENC \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_CM500k
    CHIP_CM500k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CHIP $SAMPLE_CM500k" \
    $CHIP_CM500k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_CM10k
    CHIP_CM10k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CHIP $SAMPLE_CM10k" \
    $CHIP_CM10k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_ENCODE
    CHIP_ENC="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CHIP $SAMPLE_ENCODE" \
    $CHIP_ENC \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_ENCODE "with" $SAMPLE_CM500k
    ENC_CM500k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_ENCODE $SAMPLE_CM500k" \
    $ENC_CM500k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_ENCODE "with" $SAMPLE_CM10k
    ENC_CM10k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_ENCODE $SAMPLE_CM10k" \
    $ENC_CM10k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_ENCODE "with" $SAMPLE_CHIP
    ENC_ChIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_ENCODE $SAMPLE_CHIP" \
    $ENC_ChIP \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

done


#### Extract top 40% quantile
R
dir = "/home/arendeiro/data/human/chipmentation/peaks/"
sample = "H3K4me3_K562_500k"
techniques = c("CM", "ChIP", "ChIP_Encode")
for (technique in techniques){
    D <-read.csv(paste0(dir, sample, "_", technique, "_peaks/", sample, "_", technique, "_peaks.narrowPeak"), sep='\t', header=FALSE)
    D <- D[D$V7 >= quantile(D$V7,.4) , ]
    write.table(D, paste0(dir, sample, "_", technique, "_peaks/", sample, "_", technique, "_peaks.top40.narrowPeak"), sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
}

sample = "H3K4me3_K562_10k"
techniques = c("CM")
for (technique in techniques){
    D <-read.csv(paste0(dir, sample, "_", technique, "_peaks/", sample, "_", technique, "_peaks.narrowPeak"), sep='\t', header=FALSE)
    D <- D[D$V7 >= quantile(D$V7,.4) , ]
    write.table(D, paste0(dir, sample, "_", technique, "_peaks/", sample, "_", technique, "_peaks.top40.narrowPeak"), sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
}

### Loop again for top40% peaks.
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersections.top40.txt
rm $INTERSECTIONS

for SAMPLE_NAME in H3K4me3_K562_500k_CM
    do
    SAMPLE=${SAMPLE_NAME/_CM/}
    SAMPLE_CM500k=${SAMPLE}_CM
    SAMPLE_CM10k=${SAMPLE_CM500k/500k/10k}
    SAMPLE_CHIP=${SAMPLE}_ChIP
    SAMPLE_ENCODE=${SAMPLE}_ChIP_Encode

    # Intersect and count
    echo "Intersecting" $SAMPLE_CM500k "with" $SAMPLE_CM10k
    CM500k_CM10k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak \
    | wc -l)"
    
    echo "$SAMPLE_CM500k $SAMPLE_CM10k" \
    $CM500k_CM10k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM500k "with" $SAMPLE_CHIP
    CM500k_CHIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM500k $SAMPLE_CHIP" \
    $CM500k_CHIP \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM500k "with" $SAMPLE_ENCODE
    CM500k_ENC="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM500k $SAMPLE_ENCODE" \
    $CM500k_ENC \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM10k "with" $SAMPLE_CM500k
    CM10k_CM500k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM10k $SAMPLE_CM500k" \
    $CM10k_CM500k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM10k "with" $SAMPLE_CHIP
    CM10k_CHIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM10k $SAMPLE_CHIP" \
    $CM10k_CHIP \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CM10k "with" $SAMPLE_ENCODE
    CM10k_ENC="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CM10k $SAMPLE_ENCODE" \
    $CM10k_ENC \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_CM500k
    CHIP_CM500k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CHIP $SAMPLE_CM500k" \
    $CHIP_CM500k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_CM10k
    CHIP_CM10k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CHIP $SAMPLE_CM10k" \
    $CHIP_CM10k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_ENCODE
    CHIP_ENC="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CHIP $SAMPLE_ENCODE" \
    $CHIP_ENC \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_ENCODE "with" $SAMPLE_CM500k
    ENC_CM500k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_ENCODE $SAMPLE_CM500k" \
    $ENC_CM500k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM500k}_peaks/${SAMPLE_CM500k}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_ENCODE "with" $SAMPLE_CM10k
    ENC_CM10k="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_ENCODE $SAMPLE_CM10k" \
    $ENC_CM10k \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM10k}_peaks/${SAMPLE_CM10k}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_ENCODE "with" $SAMPLE_CHIP
    ENC_ChIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_ENCODE $SAMPLE_CHIP" \
    $ENC_ChIP \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_ENCODE}_peaks/${SAMPLE_ENCODE}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

done



# BETWEEN REPLICATES
DIR=/home/arendeiro/data/human/chipmentation/peaks/individual_samples
# Extract top 40%
R1 <-read.csv("H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks.top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks.top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R1 <-read.csv("H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks.top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks.top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R1 <-read.csv("H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks.top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks.top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)


DIR=/home/arendeiro/data/human/chipmentation/peaks/individual_samples
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersectionsReplicates.txt
rm $INTERSECTIONS

# all
CM500k="$(bedtools intersect -u \
-a $DIR/H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks.narrowPeak \
-b $DIR/H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks.narrowPeak \
| wc -l)"
echo "H3K4me3_K562_500k_CM_CM11-2_R1" "H3K4me3_K562_500k_CM_CM12-2_R2" \
$CM500k \
"$(cat $DIR/H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks.narrowPeak | wc -l)" \
"$(cat $DIR/H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

CM10k="$(bedtools intersect -u \
-a $DIR/H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks.narrowPeak \
-b $DIR/H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks.narrowPeak \
| wc -l)"
echo "H3K4me3_K562_10k_CM_CM11-6_R1" "H3K4me3_K562_10k_CM_CM12-6_R2" \
$CM10k \
"$(cat $DIR/H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks.narrowPeak | wc -l)" \
"$(cat $DIR/H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

CHIP500k="$(bedtools intersect -u \
-a $DIR/H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks.narrowPeak \
-b $DIR/H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks.narrowPeak \
| wc -l)"
echo "H3K4me3_K562_500k_ChIP_CM11-10_R1" "H3K4me3_K562_500k_ChIP_CM12-10_R2" \
$CHIP500k \
"$(cat $DIR/H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks.narrowPeak | wc -l)" \
"$(cat $DIR/H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS


INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersectionsReplicates.top40.txt
rm $INTERSECTIONS

# top40%
CM500k="$(bedtools intersect -u \
-a $DIR/H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks.top40.narrowPeak \
-b $DIR/H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks.top40.narrowPeak \
| wc -l)"
echo "H3K4me3_K562_500k_CM_CM11-2_R1" "H3K4me3_K562_500k_CM_CM12-2_R2" \
$CM500k \
"$(cat $DIR/H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks.top40.narrowPeak | wc -l)" \
"$(cat $DIR/H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

CM10k="$(bedtools intersect -u \
-a $DIR/H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks.top40.narrowPeak \
-b $DIR/H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks.top40.narrowPeak \
| wc -l)"
echo "H3K4me3_K562_10k_CM_CM11-6_R1" "H3K4me3_K562_10k_CM_CM12-6_R2" \
$CM10k \
"$(cat $DIR/H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks.top40.narrowPeak | wc -l)" \
"$(cat $DIR/H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

CHIP500k="$(bedtools intersect -u \
-a $DIR/H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks.top40.narrowPeak \
-b $DIR/H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks.top40.narrowPeak \
| wc -l)"
echo "H3K4me3_K562_500k_ChIP_CM11-10_R1" "H3K4me3_K562_500k_ChIP_CM12-10_R2" \
$CHIP500k \
"$(cat $DIR/H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks.top40.narrowPeak | wc -l)" \
"$(cat $DIR/H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS





date
