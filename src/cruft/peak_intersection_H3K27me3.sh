
# BETWEEN REPLICATES (SPP)
$PROJECTDIR/spp_peaks/
# Extract top 40%
R1 <-read.csv("H3K27me3_K562_500k_ChIP_CM11-10_R1.broadPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K27me3_K562_500k_ChIP_CM11-10_R1.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K27me3_K562_500k_ChIP_CM12-10_R2.broadPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K27me3_K562_500k_ChIP_CM12-10_R2.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R1 <-read.csv("H3K27me3_K562_500k_CM_CM11-2_R1.broadPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K27me3_K562_500k_CM_CM11-2_R1.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K27me3_K562_500k_CM_CM12-2_R2.broadPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K27me3_K562_500k_CM_CM12-2_R2.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R1 <-read.csv("H3K27me3_K562_10k_CM_CM11-6_R1.broadPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K27me3_K562_10k_CM_CM11-6_R1.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K27me3_K562_10k_CM_CM12-6_R2.broadPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K27me3_K562_10k_CM_CM12-6_R2.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)


PROJECTDIR=/home/arendeiro/data/human/chipmentation
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersectionsReplicates.txt
rm $INTERSECTIONS

# all
CM500k="$(bedtools intersect -u \
-a $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_CM_CM11-2_R1.broadPeak \
-b $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_CM_CM12-2_R2.broadPeak \
| wc -l)"
echo "H3K27me3_K562_500k_CM_CM11-2_R1" "H3K27me3_K562_500k_CM_CM12-2_R2" \
$CM500k \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_CM_CM11-2_R1.broadPeak | wc -l)" \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_CM_CM12-2_R2.broadPeak | wc -l)" >> $INTERSECTIONS

CM10k="$(bedtools intersect -u \
-a $PROJECTDIR/spp_peaks/H3K27me3_K562_10k_CM_CM11-6_R1.broadPeak \
-b $PROJECTDIR/spp_peaks/H3K27me3_K562_10k_CM_CM12-6_R2.broadPeak \
| wc -l)"
echo "H3K27me3_K562_10k_CM_CM11-6_R1" "H3K27me3_K562_10k_CM_CM12-6_R2" \
$CM10k \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_10k_CM_CM11-6_R1.broadPeak | wc -l)" \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_10k_CM_CM12-6_R2.broadPeak | wc -l)" >> $INTERSECTIONS

CHIP500k="$(bedtools intersect -u \
-a $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_ChIP_CM11-10_R1.broadPeak \
-b $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_ChIP_CM12-10_R2.broadPeak \
| wc -l)"
echo "H3K27me3_K562_500k_ChIP_CM11-10_R1" "H3K27me3_K562_500k_ChIP_CM12-10_R2" \
$CHIP500k \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_ChIP_CM11-10_R1.broadPeak | wc -l)" \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_ChIP_CM12-10_R2.broadPeak | wc -l)" >> $INTERSECTIONS


INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersectionsReplicates.top40.txt
rm $INTERSECTIONS

# top40%
CM500k="$(bedtools intersect -u \
-a $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_CM_CM11-2_R1.top40.broadPeak \
-b $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_CM_CM12-2_R2.top40.broadPeak \
| wc -l)"
echo "H3K27me3_K562_500k_CM_CM11-2_R1" "H3K27me3_K562_500k_CM_CM12-2_R2" \
$CM500k \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_CM_CM11-2_R1.top40.broadPeak | wc -l)" \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_CM_CM12-2_R2.top40.broadPeak | wc -l)" >> $INTERSECTIONS

CM10k="$(bedtools intersect -u \
-a $PROJECTDIR/spp_peaks/H3K27me3_K562_10k_CM_CM11-6_R1.top40.broadPeak \
-b $PROJECTDIR/spp_peaks/H3K27me3_K562_10k_CM_CM12-6_R2.top40.broadPeak \
| wc -l)"
echo "H3K27me3_K562_10k_CM_CM11-6_R1" "H3K27me3_K562_10k_CM_CM12-6_R2" \
$CM10k \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_10k_CM_CM11-6_R1.top40.broadPeak | wc -l)" \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_10k_CM_CM12-6_R2.top40.broadPeak | wc -l)" >> $INTERSECTIONS

CHIP500k="$(bedtools intersect -u \
-a $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_ChIP_CM11-10_R1.top40.broadPeak \
-b $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_ChIP_CM12-10_R2.top40.broadPeak \
| wc -l)"
echo "H3K27me3_K562_500k_ChIP_CM11-10_R1" "H3K27me3_K562_500k_ChIP_CM12-10_R2" \
$CHIP500k \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_ChIP_CM11-10_R1.top40.broadPeak | wc -l)" \
"$(cat $PROJECTDIR/spp_peaks/H3K27me3_K562_500k_ChIP_CM12-10_R2.top40.broadPeak | wc -l)" >> $INTERSECTIONS





# BETWEEN TECHNIQUES
samples = c('H3K27me3_K562_10k_CM', 'H3K27me3_K562_500k_ChIP_Encode', 'H3K27me3_K562_500k_ChIP', 'H3K27me3_K562_500k_CM')
for (sample in samples){
R1 <-read.csv(paste0(sample, '.broadPeak'), sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, paste0(sample, '_top40.broadPeak'), sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
}


# H3K27me3
# all 
for SAMPLE1 in H3K27me3_K562_500k_CM_peaks H3K27me3_K562_10k_CM_peaks H3K27me3_K562_500k_ChIP_peaks H3K27me3_K562_500k_ChIP_Encode_peaks
do
    for SAMPLE2 in H3K27me3_K562_500k_CM_peaks H3K27me3_K562_10k_CM_peaks H3K27me3_K562_500k_ChIP_peaks H3K27me3_K562_500k_ChIP_Encode_peaks
    do
        if [ $SAMPLE1 != $SAMPLE2 ]
        then
            A1="$(bedtools intersect -u \
            -a /home/afr/Documents/workspace/chipmentation/data/$SAMPLE1/$SAMPLE1.broadPeak \
            -b /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/$SAMPLE2.broadPeak \
            | wc -l)"
            A2="$(wc -l /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/$SAMPLE2.broadPeak)"
            echo $A1 $A2
        fi
    done
done


# SPP

for SAMPLE1 in H3K27me3_K562_500k_CM H3K27me3_K562_10k_CM H3K27me3_K562_500k_ChIP H3K27me3_K562_500k_ChIP_Encode
do
    #awk -v OFS='\t' '{ print $1, $2, sprintf("%.f", $3), $4, $5, $6, $7, $8, $9; }' $SAMPLE1.broadPeak > tmp
    #mv tmp $SAMPLE1.broadPeak
    for SAMPLE2 in H3K27me3_K562_500k_CM H3K27me3_K562_10k_CM H3K27me3_K562_500k_ChIP H3K27me3_K562_500k_ChIP_Encode
    do
        if [ $SAMPLE1 != $SAMPLE2 ]
        then
            A1="$(bedtools intersect -u \
            -a $SAMPLE1.broadPeak \
            -b $SAMPLE2.broadPeak \
            | wc -l)"
            A2="$(wc -l $SAMPLE1.broadPeak)"
            A3="$(wc -l $SAMPLE2.broadPeak)"
            echo $A1 $A2 $A3
        fi
    done
done


# BETWEEN TECHNIQUES
samples = c('H3K27me3_K562_10k_CM', 'H3K27me3_K562_500k_ChIP_Encode', 'H3K27me3_K562_500k_ChIP', 'H3K27me3_K562_500k_CM')
for (sample in samples){
    R1 <-read.csv(paste0(sample, '.broadPeak'), sep='\t', header=FALSE)
    R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
    write.table(R140, paste0(sample, '_top40.broadPeak'), sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
}

# top40%
for SAMPLE1 in H3K27me3_K562_500k_CM H3K27me3_K562_10k_CM H3K27me3_K562_500k_ChIP H3K27me3_K562_500k_ChIP_Encode
do
    for SAMPLE2 in H3K27me3_K562_500k_CM H3K27me3_K562_10k_CM H3K27me3_K562_500k_ChIP H3K27me3_K562_500k_ChIP_Encode
    do
        if [ $SAMPLE1 != $SAMPLE2 ]
        then
            A1="$(bedtools intersect -u \
            -a ${SAMPLE1}_top40.broadPeak \
            -b ${SAMPLE2}_top40.broadPeak \
            | wc -l)"
            A2="$(wc -l ${SAMPLE1}_top40.broadPeak)"
            A3="$(wc -l ${SAMPLE2}_top40.broadPeak)"
            echo $A1 $A2 $A3
        fi
    done
done







# BETWEEN REPLICATES (MACS)
DIR=/home/arendeiro/data/human/chipmentation/peaks/individual_samples
# Extract top 40%
R1 <-read.csv("H3K27me3_K562_500k_CM_CM11-2_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_CM_CM11-2_R1_std_peaks.broadPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K27me3_K562_500k_CM_CM11-2_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_CM_CM11-2_R1_std_peaks.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K27me3_K562_500k_CM_CM12-2_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_CM_CM12-2_R2_std_peaks.broadPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K27me3_K562_500k_CM_CM12-2_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_CM_CM12-2_R2_std_peaks.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R1 <-read.csv("H3K27me3_K562_10k_CM_CM11-6_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_10k_CM_CM11-6_R1_std_peaks.broadPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K27me3_K562_10k_CM_CM11-6_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_10k_CM_CM11-6_R1_std_peaks.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K27me3_K562_10k_CM_CM12-6_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_10k_CM_CM12-6_R2_std_peaks.broadPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K27me3_K562_10k_CM_CM12-6_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_10k_CM_CM12-6_R2_std_peaks.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R1 <-read.csv("H3K27me3_K562_500k_ChIP_CM11-10_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_ChIP_CM11-10_R1_std_peaks.broadPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K27me3_K562_500k_ChIP_CM11-10_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_ChIP_CM11-10_R1_std_peaks.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K27me3_K562_500k_ChIP_CM12-10_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_ChIP_CM12-10_R2_std_peaks.broadPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K27me3_K562_500k_ChIP_CM12-10_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_ChIP_CM12-10_R2_std_peaks.top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)



DIR=/home/arendeiro/data/human/chipmentation/peaks/individual_samples
INTERSECTIONS=/home/arendeiro/projects/chipmentation/results/peak_intersectionsReplicates.txt
rm $INTERSECTIONS

# all
CM500k="$(bedtools intersect -u \
-a $DIR/H3K27me3_K562_500k_CM_CM11-2_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_CM_CM11-2_R1_std_peaks.broadPeak \
-b $DIR/H3K27me3_K562_500k_CM_CM12-2_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_CM_CM12-2_R2_std_peaks.broadPeak \
| wc -l)"
echo "H3K4me3_K562_500k_CM_CM11-2_R1" "H3K4me3_K562_500k_CM_CM12-2_R2" \
$CM500k \
"$(cat $DIR/H3K27me3_K562_500k_CM_CM11-2_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_CM_CM11-2_R1_std_peaks.broadPeak | wc -l)" \
"$(cat $DIR/H3K27me3_K562_500k_CM_CM12-2_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_CM_CM12-2_R2_std_peaks.broadPeak | wc -l)" >> $INTERSECTIONS

CM10k="$(bedtools intersect -u \
-a $DIR/H3K27me3_K562_10k_CM_CM11-6_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_10k_CM_CM11-6_R1_std_peaks.broadPeak \
-b $DIR/H3K27me3_K562_10k_CM_CM12-6_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_10k_CM_CM12-6_R2_std_peaks.broadPeak \
| wc -l)"
echo "H3K4me3_K562_10k_CM_CM11-6_R1" "H3K4me3_K562_10k_CM_CM12-6_R2" \
$CM10k \
"$(cat $DIR/H3K27me3_K562_10k_CM_CM11-6_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_10k_CM_CM11-6_R1_std_peaks.broadPeak | wc -l)" \
"$(cat $DIR/H3K27me3_K562_10k_CM_CM12-6_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_10k_CM_CM12-6_R2_std_peaks.broadPeak | wc -l)" >> $INTERSECTIONS

CHIP500k="$(bedtools intersect -u \
-a $DIR/H3K27me3_K562_500k_ChIP_CM11-10_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_ChIP_CM11-10_R1_std_peaks.broadPeak \
-b $DIR/H3K27me3_K562_500k_ChIP_CM12-10_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_ChIP_CM12-10_R2_std_peaks.broadPeak \
| wc -l)"
echo "H3K4me3_K562_500k_CM_CM11-2_R1" "H3K4me3_K562_500k_CM_CM12-2_R2" \
$CHIP500k \
"$(cat $DIR/H3K27me3_K562_500k_ChIP_CM11-10_R1_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_ChIP_CM11-10_R1_std_peaks.broadPeak | wc -l)" \
"$(cat $DIR/H3K27me3_K562_500k_ChIP_CM12-10_R2_peaks_MACS2_broad_nomodel_extsize200_pvalue1e3/H3K27me3_K562_500k_ChIP_CM12-10_R2_std_peaks.broadPeak | wc -l)" >> $INTERSECTIONS


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
