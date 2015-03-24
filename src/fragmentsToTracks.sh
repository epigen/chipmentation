#!/usr/bin/env bash

splitBam() {
    NAME=$1
    samtools view -h /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.shifted.bam | \
    awk '{
        abs=($9<0?-$9:$9)
        if( $1 ~ /^\@/ ) {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.1.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.1.sam" )
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.2.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.2.sam" )
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.3.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.3.sam" )
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.4.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.4.sam" )
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.5.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.5.sam" )
        }
        else if( abs < 63 ) {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.1.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.1.sam" )
        }
        else if( 64 < abs  && abs < 81 ) {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.2.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.2.sam" )
        }
        else if( 82 < abs  && abs < 103 ) {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.3.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.3.sam" )
        }
        else if( 104 < abs  && abs < 137 ) {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.4.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.4.sam" )
        }
        else {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.5.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.5.sam" )
        }
    }'
}
export -f splitBam

# for NAME in ${SAMPLES[*]}
# do
#     head -n 27 /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.4.sam > header
#     cat header /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.5.sam > tmp
#     mv tmp /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.5.sam
# done

replaceHeader() {
    NAME=$1
    sed -i .bak -e '1s/SO:coordinate/SO:unsorted/' /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.1.sam
    sed -i .bak -e '1s/SO:coordinate/SO:unsorted/' /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.2.sam
    sed -i .bak -e '1s/SO:coordinate/SO:unsorted/' /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.3.sam
    sed -i .bak -e '1s/SO:coordinate/SO:unsorted/' /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.4.sam
    sed -i .bak -e '1s/SO:coordinate/SO:unsorted/' /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.5.sam
}
export -f replaceHeader


samToBigWig() {
    NAME=$1
    samtools view -S -b /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/$NAME.sam | \
    bedtools bamtobed | \
    bedtools slop -i stdin -g /fhgfs/groups/lab_bock/arendeiro/share/hg19.chrom.sizes -s -l 0 -r 130 | \
    python /home/arendeiro/chipseq-pipelines/lib/fix_bedfile_genome_boundaries.py hg19 | \
    genomeCoverageBed -i stdin -bg -g /fhgfs/groups/lab_bock/arendeiro/share/hg19.chrom.sizes > \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.cov

    bedGraphToBigWig /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.cov \
    /fhgfs/groups/lab_bock/arendeiro/share/hg19.chrom.sizes \
    /fhgfs/groups/lab_bock/public_html/arendeiro/chipmentation/pe/${NAME}.bigWig

    chmod 655 /fhgfs/groups/lab_bock/public_html/arendeiro/chipmentation/pe/${NAME}.bigWig 
}
export -f samToBigWig


samToUCSC() {
    NAME=$1

    for N in ${NAME}.1 ${NAME}.2 ${NAME}.3 ${NAME}.4 ${NAME}.5
    do
         samToBigWig $N
    done
}
export -f samToUCSC


addTrackToHub() {
    NAME=$1
    echo "track type=bigWig name='"${NAME}"' description='"${NAME}"' height=32 visibility=full maxHeightPixels=32:32:25 bigDataUrl=http://www.biomedical-sequencing.at/bocklab/arendeiro/chipmentation/pe/"${NAME}".bigWig color=88,61,61" >> \
    /fhgfs/groups/lab_bock/public_html/arendeiro/chipmentation/pe/trackHub.txt
}
export -f addTrackToHub


samToBam2() {
    NAME=$1
    samtools view -S -b \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.sam \
    -o /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.bam

    samtools sort /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.bam \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.sorted
    
    samtools index /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.sorted.bam
}
export -f samToBam2


samToBam() {
    NAME=$1

    for N in ${NAME}.1 ${NAME}.2 ${NAME}.3 ${NAME}.4 ${NAME}.5
    do
        samToBam2 $N
    done
}
export -f samToBam


countReadsInPeaks() {
    NAME=$1

    bedtools multicov -bams /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.1.sorted.bam \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.2.sorted.bam \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.3.sorted.bam \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.4.sorted.bam \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.5.sorted.bam \
    -bed /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/${NAME}/${NAME}_peaks.narrowPeak > \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.readCount
}
export -f countReadsInPeaks


declare -a SAMPLES

SAMPLES=(K562_10M_ATAC_H3K4ME1_nan_PE_1_1_hg19 K562_10M_ATAC_PU1_nan_PE_1_1_hg19 
    K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19 K562_10M_CM_PU1_nan_PE_1_1_hg19 
    K562_500K_ATAC_H3K4ME3_nan_01ULTN5_PE_1_1_hg19 K562_500K_CM_H3K4ME3_nan_02ULTN5_PE_1_1_hg19 
    K562_500K_CM_H3K4ME3_nan_05ULTN5_PE_1_1_hg19 K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19 
    K562_500K_CM_H3K4ME3_nan_5ULTN5_PE_1_1_hg19 K562_500K_CM_IGG_nan_1ULTN5_PE_1_1_hg19)

parallel splitBam ::: ${SAMPLES[*]}

parallel replaceHeader ::: ${SAMPLES[*]}

parallel samToUCSC ::: ${SAMPLES[*]}

for SAMPLE in ${SAMPLES[@]}
do
    for N in ${SAMPLE}.1 ${SAMPLE}.2 ${SAMPLE}.3 ${SAMPLE}.4 ${SAMPLE}.5
    do
         addTrackToHub $N
    done
done

chmod 655 /fhgfs/groups/lab_bock/public_html/arendeiro/chipmentation/pe/trackHub.txt

chmod 655 /fhgfs/groups/lab_bock/public_html/arendeiro/chipmentation/pe/

# Go to:
# http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&hgt.customText=http://www.biomedical-sequencing.at/bocklab/arendeiro/chipmentation/pe/trackHub.txt


#### Get FRiP per fragment

parallel samToBam ::: ${SAMPLES[*]}

SAMPLES=(K562_10M_ATAC_H3K4ME1_nan_PE_1_1_hg19 K562_10M_ATAC_PU1_nan_PE_1_1_hg19 
    K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19 K562_10M_CM_PU1_nan_PE_1_1_hg19 
    K562_500K_ATAC_H3K4ME3_nan_01ULTN5_PE_1_1_hg19 K562_500K_CM_H3K4ME3_nan_02ULTN5_PE_1_1_hg19 
    K562_500K_CM_H3K4ME3_nan_05ULTN5_PE_1_1_hg19 K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19 
    K562_500K_CM_H3K4ME3_nan_5ULTN5_PE_1_1_hg19)

parallel countReadsInPeaks ::: ${SAMPLES[*]}


# Make table with FRiP per fragment length
echo name fraction value > /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/FRiP_per_fragment_length.tsv

for NAME in ${SAMPLES[*]}
do
    SUMS=($(awk -v OFS='\t' '{sum1+=$11; sum2+=$12; sum3+=$13; sum4+=$14; sum5+=$15} END {print sum1,sum2,sum3,sum4,sum5}' \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.readCount))
    #echo $SUMS
    for N in 1 2 3 4 5
    do
        TOTAL=`samtools idxstats /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.${N}.sorted.bam | \
        awk '{sum+=$3} END {print sum}'`
        CUR=$((N - 1))
        SUM=${SUMS[$CUR]}
        echo $NAME $N $(awk "BEGIN {printf \"%.4f\",${SUM}/${TOTAL}}") >> \
        /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/FRiP_per_fragment_length.tsv
    done
done

### R code to plot
# library(reshape)
# df = read.table("/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/FRiP_per_fragment_length.tsv", header=TRUE)

# library(ggplot2)

# p = ggplot(df, aes(fraction, value)) +
#     geom_line() +
#     facet_wrap(~name, scales="free") +
#     xlab("Fragment length bin") +
#     ylab("FRiP") +
#     theme_bw() +
#     theme(legend.title=element_blank())

# ggsave(
#     paste0(
#     "/fhgfs/groups/lab_bock/shared/projects/chipmentation/results/pe/plots",
#     "FRiP_per_fragment_length.pdf"),
#     plot = p, height = 12, width = 15
# )
