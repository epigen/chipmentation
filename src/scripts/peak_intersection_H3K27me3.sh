
# BETWEEN REPLICATES
# CHIP
R1 <-read.csv("H3K27me3_K562_500k_ChIP_CM11-10_R1.broadPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K27me3_K562_500k_ChIP_CM11-10_R1_top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K27me3_K562_500k_ChIP_CM12-10_R2.broadPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K27me3_K562_500k_ChIP_CM12-10_R2_top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)


# all
bedtools intersect -u \
-a H3K27me3_K562_500k_ChIP_CM11-10_R1.broadPeak \
-b H3K27me3_K562_500k_ChIP_CM12-10_R2.broadPeak \
| wc -l
wc -l H3K27me3_K562_500k_ChIP_CM11-10_R1.broadPeak
wc -l H3K27me3_K562_500k_ChIP_CM12-10_R2.broadPeak

# top40%
bedtools intersect -u \
-a H3K27me3_K562_500k_ChIP_CM11-10_R1_top40.broadPeak \
-b H3K27me3_K562_500k_ChIP_CM12-10_R2_top40.broadPeak \
| wc -l
wc -l H3K27me3_K562_500k_ChIP_CM11-10_R1_top40.broadPeak
wc -l H3K27me3_K562_500k_ChIP_CM12-10_R2_top40.broadPeak

# CM 500
R1 <-read.csv("H3K27me3_K562_500k_CM_CM11-2_R1.broadPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K27me3_K562_500k_CM_CM11-2_R1_top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K27me3_K562_500k_CM_CM12-2_R2.broadPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K27me3_K562_500k_CM_CM12-2_R2_top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

# all
bedtools intersect -u \
-a H3K27me3_K562_500k_CM_CM11-2_R1.broadPeak \
-b H3K27me3_K562_500k_CM_CM12-2_R2.broadPeak \
| wc -l

wc -l H3K27me3_K562_500k_CM_CM11-2_R1.broadPeak
wc -l H3K27me3_K562_500k_CM_CM12-2_R2.broadPeak

# 40%
bedtools intersect -u \
-a H3K27me3_K562_500k_CM_CM11-2_R1_top40.broadPeak \
-b H3K27me3_K562_500k_CM_CM12-2_R2_top40.broadPeak \
| wc -l

wc -l H3K27me3_K562_500k_CM_CM11-2_R1_top40.broadPeak
wc -l H3K27me3_K562_500k_CM_CM12-2_R2_top40.broadPeak

# CM 10k
R1 <-read.csv("H3K27me3_K562_10k_CM_CM11-6_R1.broadPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "H3K27me3_K562_10k_CM_CM11-6_R1_top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K27me3_K562_10k_CM_CM12-6_R2.broadPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "H3K27me3_K562_10k_CM_CM12-6_R2_top40.broadPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

# all
bedtools intersect -u \
-a H3K27me3_K562_10k_CM_CM11-6_R1.broadPeak \
-b H3K27me3_K562_10k_CM_CM12-6_R2.broadPeak \
| wc -l

wc -l H3K27me3_K562_10k_CM_CM11-6_R1.broadPeak
wc -l H3K27me3_K562_10k_CM_CM12-6_R2.broadPeak

# top40%
bedtools intersect -u \
-a H3K27me3_K562_10k_CM_CM11-6_R1_top40.broadPeak \
-b H3K27me3_K562_10k_CM_CM12-6_R2_top40.broadPeak \
| wc -l

wc -l H3K27me3_K562_10k_CM_CM11-6_R1_top40.broadPeak
wc -l H3K27me3_K562_10k_CM_CM12-6_R2_top40.broadPeak




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


