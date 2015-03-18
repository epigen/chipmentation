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
        }
        else if( abs < 98 ) {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.1.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.1.sam" )
        }
        else if( 99 < abs  && abs < 112 ) {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.2.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.2.sam" )
        }
        else if( 113 < abs  && abs < 140 ) {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.3.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.3.sam" )
        }
        else {
            print >> "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.4.sam"
            close( "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/'$NAME'.4.sam" )
        }
    }'
}
export -f splitBam


replaceHeader() {
    NAME=$1
    sed -i .bak -e '1s/SO:coordinate/SO:unsorted/' /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.1.sam
    sed -i .bak -e '1s/SO:coordinate/SO:unsorted/' /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.2.sam
    sed -i .bak -e '1s/SO:coordinate/SO:unsorted/' /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.3.sam
    sed -i .bak -e '1s/SO:coordinate/SO:unsorted/' /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/pe/${NAME}.4.sam
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

    for N in ${NAME}.1 ${NAME}.2 ${NAME}.3 ${NAME}.4
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
    for N in ${SAMPLE}.1 ${SAMPLE}.2 ${SAMPLE}.3 ${SAMPLE}.4
    do
         addTrackToHub $N
    done
done

chmod 655 /fhgfs/groups/lab_bock/public_html/arendeiro/chipmentation/pe/trackHub.txt

chmod 655 /fhgfs/groups/lab_bock/public_html/arendeiro/chipmentation/pe/

# Go to:
# http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&hgt.customText=http://www.biomedical-sequencing.at/bocklab/arendeiro/chipmentation/pe/trackHub.txt
