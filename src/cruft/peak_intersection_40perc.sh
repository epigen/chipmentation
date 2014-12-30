
# BETWEEN REPLICATES
# CHIP
R1 <-read.csv("H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "/home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks_top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "/home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks_top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)


# all
bedtools intersect -u \
-a /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks.narrowPeak \
-b /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks.narrowPeak \
| wc -l
wc -l /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks.narrowPeak

# top40%
bedtools intersect -u \
-a /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks_top40.narrowPeak \
-b /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_ChIP_CM12-9_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM12-9_R2_bw200_peaks_top40.narrowPeak \
| wc -l
wc -l /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_ChIP_CM11-9_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_ChIP_CM11-9_R1_bw200_peaks_top40.narrowPeak


# CM 500
R1 <-read.csv("H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "/home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks_top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "/home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks_top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

# all
bedtools intersect -u \
-a /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks.narrowPeak \
-b /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks.narrowPeak \
| wc -l

wc -l /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks.narrowPeak


bedtools intersect -u \
-a /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks_top40.narrowPeak \
-b /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_CM_CM12-1_R2_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM12-1_R2_bw200_peaks_top40.narrowPeak \
| wc -l

wc -l /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_500k_CM_CM11-1_R1_peaks_MACS2_bw200/H3K4me3_K562_500k_CM_CM11-1_R1_bw200_peaks_top40.narrowPeak


# CM 10k
R1 <-read.csv("H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, "/home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks_top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

R2 <-read.csv("H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks.narrowPeak", sep='\t', header=FALSE)
R240 <- R2[R2$V7 >= quantile(R2$V7,.4) , ]
write.table(R240, "/home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks_top40.narrowPeak", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

# all
bedtools intersect -u \
-a /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks.narrowPeak \
-b /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks.narrowPeak \
| wc -l

wc -l /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks.narrowPeak

# top40%
bedtools intersect -u \
-a /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks_top40.narrowPeak \
-b /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_10k_CM_CM12-5_R2_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM12-5_R2_bw200_peaks_top40.narrowPeak \
| wc -l

wc -l /home/afr/Documents/workspace/chipmentation/data/H3K4me3_K562_10k_CM_CM11-5_R1_peaks_MACS2_bw200/H3K4me3_K562_10k_CM_CM11-5_R1_bw200_peaks_top40.narrowPeak





# BETWEEN TECHNIQUES
samples = c('H3K4me3_K562_10k_CM_peaks', 'H3K4me3_K562_500k_ChIP_Encode_peaks', 'H3K4me3_K562_500k_ChIP_peaks', 'H3K4me3_K562_500k_CM_peaks',
    'PU1_K562_10mio_CM_peaks', 'PU1_K562_10mio_ChIP_peaks',
    'CTCF_K562_10mio_ChIP_peaks', 'CTCF_K562_10mio_CM_peaks')
for (sample in samples){
R1 <-read.csv(paste('/home/afr/Documents/workspace/chipmentation/data/', sample, '/', sample, '.narrowPeak', sep=''), sep='\t', header=FALSE)
R140 <- R1[R1$V7 >= quantile(R1$V7,.4) , ]
write.table(R140, paste("/home/afr/Documents/workspace/chipmentation/data/", sample, '/', sample, '_top40.narrowPeak', sep=''), sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
}


# H3K4me3
# all 
for SAMPLE1 in H3K4me3_K562_500k_CM_peaks H3K4me3_K562_10k_CM_peaks H3K4me3_K562_500k_ChIP_peaks H3K4me3_K562_500k_ChIP_Encode_peaks
do
    for SAMPLE2 in H3K4me3_K562_500k_CM_peaks H3K4me3_K562_10k_CM_peaks H3K4me3_K562_500k_ChIP_peaks H3K4me3_K562_500k_ChIP_Encode_peaks
    do
        if [ $SAMPLE1 != $SAMPLE2 ]
        then
            A1="$(bedtools intersect -u \
            -a /home/afr/Documents/workspace/chipmentation/data/$SAMPLE1/$SAMPLE1.narrowPeak \
            -b /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/$SAMPLE2.narrowPeak \
            | wc -l)"
            A2="$(wc -l /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/$SAMPLE2.narrowPeak)"
            echo $A1 $A2
        fi
    done
done

# top40%
for SAMPLE1 in H3K4me3_K562_500k_CM_peaks H3K4me3_K562_10k_CM_peaks H3K4me3_K562_500k_ChIP_peaks H3K4me3_K562_500k_ChIP_Encode_peaks
do
    for SAMPLE2 in H3K4me3_K562_500k_CM_peaks H3K4me3_K562_10k_CM_peaks H3K4me3_K562_500k_ChIP_peaks H3K4me3_K562_500k_ChIP_Encode_peaks
    do
        if [ $SAMPLE1 != $SAMPLE2 ]
        then
            A1="$(bedtools intersect -u \
            -a /home/afr/Documents/workspace/chipmentation/data/$SAMPLE1/${SAMPLE1}_top40.narrowPeak \
            -b /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/${SAMPLE2}_top40.narrowPeak \
            | wc -l)"
            A2="$(wc -l /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/${SAMPLE2}_top40.narrowPeak)"
            echo $A1 $A2
        fi
    done
done


# PU1
# all 
for SAMPLE1 in PU1_K562_10mio_CM_peaks PU1_K562_10mio_ChIP_peaks
do
    for SAMPLE2 in PU1_K562_10mio_CM_peaks PU1_K562_10mio_ChIP_peaks
    do
        if [ $SAMPLE1 != $SAMPLE2 ]
        then
            A1="$(bedtools intersect -u \
            -a /home/afr/Documents/workspace/chipmentation/data/$SAMPLE1/$SAMPLE1.narrowPeak \
            -b /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/$SAMPLE2.narrowPeak \
            | wc -l)"
            A2="$(wc -l /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/$SAMPLE2.narrowPeak)"
            echo $A1 $A2
        fi
    done
done

# top40%
for SAMPLE1 in PU1_K562_10mio_CM_peaks PU1_K562_10mio_ChIP_peaks
do
    for SAMPLE2 in PU1_K562_10mio_CM_peaks PU1_K562_10mio_ChIP_peaks
    do
        if [ $SAMPLE1 != $SAMPLE2 ]
        then
            A1="$(bedtools intersect -u \
            -a /home/afr/Documents/workspace/chipmentation/data/$SAMPLE1/${SAMPLE1}_top40.narrowPeak \
            -b /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/${SAMPLE2}_top40.narrowPeak \
            | wc -l)"
            A2="$(wc -l /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/${SAMPLE2}_top40.narrowPeak)"
            echo $A1 $A2
        fi
    done
done

# CTCF
# all 
for SAMPLE1 in CTCF_K562_10mio_ChIP_peaks CTCF_K562_10mio_CM_peaks
do
    for SAMPLE2 in CTCF_K562_10mio_ChIP_peaks CTCF_K562_10mio_CM_peaks
    do
        if [ $SAMPLE1 != $SAMPLE2 ]
        then
            A1="$(bedtools intersect -u \
            -a /home/afr/Documents/workspace/chipmentation/data/$SAMPLE1/$SAMPLE1.narrowPeak \
            -b /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/$SAMPLE2.narrowPeak \
            | wc -l)"
            A2="$(wc -l /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/$SAMPLE2.narrowPeak)"
            echo $A1 $A2
        fi
    done
done

# top40%
for SAMPLE1 in CTCF_K562_10mio_ChIP_peaks CTCF_K562_10mio_CM_peaks
do
    for SAMPLE2 in CTCF_K562_10mio_ChIP_peaks CTCF_K562_10mio_CM_peaks
    do
        if [ $SAMPLE1 != $SAMPLE2 ]
        then
            A1="$(bedtools intersect -u \
            -a /home/afr/Documents/workspace/chipmentation/data/$SAMPLE1/${SAMPLE1}_top40.narrowPeak \
            -b /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/${SAMPLE2}_top40.narrowPeak \
            | wc -l)"
            A2="$(wc -l /home/afr/Documents/workspace/chipmentation/data/$SAMPLE2/${SAMPLE2}_top40.narrowPeak)"
            echo $A1 $A2
        fi
    done
done
