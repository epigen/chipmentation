#!/bin/bash

PROJECTDIR=/home/arendeiro/data/human/chipmentation
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersectionsTFs.txt

for SAMPLE_NAME in CTCF_K562_10mio_CM PU1_K562_10mio_CM
    do
    SAMPLE=${SAMPLE_NAME/_CM/}
    SAMPLE_CM=${SAMPLE}_CM
    SAMPLE_CHIP=${SAMPLE}_ChIP

    # Intersect and count
    echo "Intersecting" $SAMPLE_CM "with" $SAMPLE_CHIP
    CM_CHIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    | wc -l)"
    
    echo "$SAMPLE_CM $SAMPLE_CHIP" \
    $CM_CHIP \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_CM
    CHIP_CM="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CHIP $SAMPLE_CM" \
    $CHIP_CM \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

done


#### Extract top 40% quantile
R
dir = "/home/arendeiro/data/human/chipmentation/peaks/"
samples = c("PU1_K562_10mio", "CTCF_K562_10mio")
techniques = c("CM", "ChIP")

for (sample in samples){
    for (technique in techniques){
        D <-read.csv(paste0(dir, sample, "_", technique, "_peaks/", sample, "_", technique, "_peaks.narrowPeak"), sep='\t', header=FALSE)
        D <- D[D$V7 >= quantile(D$V7,.4) , ]
        write.table(D, paste0(dir, sample, "_", technique, "_peaks/", sample, "_", technique, "_peaks.top40.narrowPeak"), sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
    }
}

### Loop again for top40% peaks.
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersectionsTFs.top40.txt
rm $INTERSECTIONS

for SAMPLE_NAME in PU1_K562_10mio_CM CTCF_K562_10mio_CM
    do
    SAMPLE=${SAMPLE_NAME/_CM/}
    SAMPLE_CM=${SAMPLE}_CM
    SAMPLE_CHIP=${SAMPLE}_ChIP

    # Intersect and count
    echo "Intersecting" $SAMPLE_CM "with" $SAMPLE_CHIP
    CM_CHIP="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak \
    | wc -l)"
    
    echo "$SAMPLE_CM $SAMPLE_CHIP" \
    $CM_CHIP \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

    echo "Intersecting" $SAMPLE_CHIP "with" $SAMPLE_CM
    CHIP_CM="$(bedtools intersect -u \
    -a $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak \
    -b $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.top40.narrowPeak \
    | wc -l)"

    echo "$SAMPLE_CHIP $SAMPLE_CM" \
    $CHIP_CM \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CHIP}_peaks/${SAMPLE_CHIP}_peaks.top40.narrowPeak | wc -l)" \
    "$(cat $PROJECTDIR/peaks/${SAMPLE_CM}_peaks/${SAMPLE_CM}_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

done





# BETWEEN REPLICATES
DIR=/home/arendeiro/data/human/chipmentation/peaks/individual_samples
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersectionsTFsReplicates.txt
rm $INTERSECTIONS

# all
CM10mio="$(bedtools intersect -u \
-a $DIR/PU1_K562_10mio_CM_CM15-1_R1_peaks/PU1_K562_10mio_CM_CM15-1_R1_peaks.narrowPeak \
-b $DIR/PU1_K562_10mio_CM_CM16-1_R2_peaks/PU1_K562_10mio_CM_CM16-1_R2_peaks.narrowPeak \
| wc -l)"
echo "PU1_K562_10mio_CM_CM15-1_R1" "PU1_K562_10mio_CM_CM16-1_R2" \
$CM10mio \
"$(cat $DIR/PU1_K562_10mio_CM_CM15-1_R1_peaks/PU1_K562_10mio_CM_CM15-1_R1_peaks.narrowPeak | wc -l)" \
"$(cat $DIR/PU1_K562_10mio_CM_CM16-1_R2_peaks/PU1_K562_10mio_CM_CM16-1_R2_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

ChIP10mio="$(bedtools intersect -u \
-a $DIR/PU1_K562_10mio_ChIP_CM15-5_R1_peaks/PU1_K562_10mio_ChIP_CM15-5_R1_peaks.narrowPeak \
-b $DIR/PU1_K562_10mio_ChIP_CM16-4_R2_peaks/PU1_K562_10mio_ChIP_CM16-4_R2_peaks.narrowPeak \
| wc -l)"
echo "PU1_K562_10mio_ChIP_CM15-5_R1" "PU1_K562_10mio_ChIP_CM16-4_R2" \
$ChIP10mio \
"$(cat $DIR/PU1_K562_10mio_ChIP_CM15-5_R1_peaks/PU1_K562_10mio_ChIP_CM15-5_R1_peaks.narrowPeak | wc -l)" \
"$(cat $DIR/PU1_K562_10mio_ChIP_CM16-4_R2_peaks/PU1_K562_10mio_ChIP_CM16-4_R2_peaks.narrowPeak  | wc -l)" >> $INTERSECTIONS

CM10mio="$(bedtools intersect -u \
-a $DIR/CTCF_K562_10mio_CM_CM15-2_R1_peaks/CTCF_K562_10mio_CM_CM15-2_R1_peaks.narrowPeak \
-b $DIR/CTCF_K562_10mio_CM_CM16-2_R2_peaks/CTCF_K562_10mio_CM_CM16-2_R2_peaks.narrowPeak \
| wc -l)"
echo "CTCF_K562_10mio_CM_CM15-2_R1" "CTCF_K562_10mio_CM_CM16-2_R2" \
$CM10mio \
"$(cat $DIR/CTCF_K562_10mio_CM_CM15-2_R1_peaks/CTCF_K562_10mio_CM_CM15-2_R1_peaks.narrowPeak | wc -l)" \
"$(cat $DIR/CTCF_K562_10mio_CM_CM16-2_R2_peaks/CTCF_K562_10mio_CM_CM16-2_R2_peaks.narrowPeak | wc -l)" >> $INTERSECTIONS

ChIP10mio="$(bedtools intersect -u \
-a $DIR/CTCF_K562_10mio_ChIP_CM15-6_R1_peaks/CTCF_K562_10mio_ChIP_CM15-6_R1_peaks.narrowPeak \
-b $DIR/CTCF_K562_10mio_ChIP_CM16-5_R2_peaks/CTCF_K562_10mio_ChIP_CM16-5_R2_peaks.narrowPeak \
| wc -l)"
echo "CTCF_K562_10mio_ChIP_CM15-6_R1" "CTCF_K562_10mio_ChIP_CM16-5_R2" \
$ChIP10mio \
"$(cat $DIR/CTCF_K562_10mio_ChIP_CM15-6_R1_peaks/CTCF_K562_10mio_ChIP_CM15-6_R1_peaks.narrowPeak | wc -l)" \
"$(cat $DIR/CTCF_K562_10mio_ChIP_CM16-5_R2_peaks/CTCF_K562_10mio_ChIP_CM16-5_R2_peaks.narrowPeak  | wc -l)" >> $INTERSECTIONS



# top40%
samples = c("PU1_K562_10mio_ChIP_CM15-5_R1",
"PU1_K562_10mio_ChIP_CM16-4_R2",
"PU1_K562_10mio_CM_CM15-1_R1",
"PU1_K562_10mio_CM_CM16-1_R2",
"CTCF_K562_10mio_ChIP_CM15-6_R1",
"CTCF_K562_10mio_ChIP_CM16-5_R2",
"CTCF_K562_10mio_CM_CM15-2_R1",
"CTCF_K562_10mio_CM_CM16-2_R2")

# Extract top 40%
for (sample in samples) {
	D <-read.csv(paste0(sample, "_peaks/", sample, "_peaks.narrowPeak"), sep='\t', header=FALSE)
	D <- D[D$V7 >= quantile(D$V7,.4) , ]
	write.table(D, paste0(sample, "_peaks/", sample, "_peaks.top40.narrowPeak"), sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
}


DIR=/home/arendeiro/data/human/chipmentation/peaks/individual_samples
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersectionsTFsReplicates.top40.txt
rm $INTERSECTIONS

CM10mio="$(bedtools intersect -u \
-a $DIR/PU1_K562_10mio_CM_CM15-1_R1_peaks/PU1_K562_10mio_CM_CM15-1_R1_peaks.top40.narrowPeak \
-b $DIR/PU1_K562_10mio_CM_CM16-1_R2_peaks/PU1_K562_10mio_CM_CM16-1_R2_peaks.top40.narrowPeak \
| wc -l)"
echo "PU1_K562_10mio_CM_CM15-1_R1" "PU1_K562_10mio_CM_CM16-1_R2" \
$CM10mio \
"$(cat $DIR/PU1_K562_10mio_CM_CM15-1_R1_peaks/PU1_K562_10mio_CM_CM15-1_R1_peaks.top40.narrowPeak | wc -l)" \
"$(cat $DIR/PU1_K562_10mio_CM_CM16-1_R2_peaks/PU1_K562_10mio_CM_CM16-1_R2_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

ChIP10mio="$(bedtools intersect -u \
-a $DIR/PU1_K562_10mio_ChIP_CM15-5_R1_peaks/PU1_K562_10mio_ChIP_CM15-5_R1_peaks.top40.narrowPeak \
-b $DIR/PU1_K562_10mio_ChIP_CM16-4_R2_peaks/PU1_K562_10mio_ChIP_CM16-4_R2_peaks.top40.narrowPeak \
| wc -l)"
echo "PU1_K562_10mio_ChIP_CM15-5_R1" "PU1_K562_10mio_ChIP_CM16-4_R2" \
$ChIP10mio \
"$(cat $DIR/PU1_K562_10mio_ChIP_CM15-5_R1_peaks/PU1_K562_10mio_ChIP_CM15-5_R1_peaks.top40.narrowPeak | wc -l)" \
"$(cat $DIR/PU1_K562_10mio_ChIP_CM16-4_R2_peaks/PU1_K562_10mio_ChIP_CM16-4_R2_peaks.top40.narrowPeak  | wc -l)" >> $INTERSECTIONS

CM10mio="$(bedtools intersect -u \
-a $DIR/CTCF_K562_10mio_CM_CM15-2_R1_peaks/CTCF_K562_10mio_CM_CM15-2_R1_peaks.top40.narrowPeak \
-b $DIR/CTCF_K562_10mio_CM_CM16-2_R2_peaks/CTCF_K562_10mio_CM_CM16-2_R2_peaks.top40.narrowPeak \
| wc -l)"
echo "CTCF_K562_10mio_CM_CM15-2_R1" "CTCF_K562_10mio_CM_CM16-2_R2" \
$CM10mio \
"$(cat $DIR/CTCF_K562_10mio_CM_CM15-2_R1_peaks/CTCF_K562_10mio_CM_CM15-2_R1_peaks.top40.narrowPeak | wc -l)" \
"$(cat $DIR/CTCF_K562_10mio_CM_CM16-2_R2_peaks/CTCF_K562_10mio_CM_CM16-2_R2_peaks.top40.narrowPeak | wc -l)" >> $INTERSECTIONS

ChIP10mio="$(bedtools intersect -u \
-a $DIR/CTCF_K562_10mio_ChIP_CM15-6_R1_peaks/CTCF_K562_10mio_ChIP_CM15-6_R1_peaks.top40.narrowPeak \
-b $DIR/CTCF_K562_10mio_ChIP_CM16-5_R2_peaks/CTCF_K562_10mio_ChIP_CM16-5_R2_peaks.top40.narrowPeak \
| wc -l)"
echo "CTCF_K562_10mio_ChIP_CM15-6_R1" "CTCF_K562_10mio_ChIP_CM16-5_R2" \
$ChIP10mio \
"$(cat $DIR/CTCF_K562_10mio_ChIP_CM15-6_R1_peaks/CTCF_K562_10mio_ChIP_CM15-6_R1_peaks.top40.narrowPeak | wc -l)" \
"$(cat $DIR/CTCF_K562_10mio_ChIP_CM16-5_R2_peaks/CTCF_K562_10mio_ChIP_CM16-5_R2_peaks.top40.narrowPeak  | wc -l)" >> $INTERSECTIONS



date


