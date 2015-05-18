# Split libraries by fragment length

splitHigh() {
    NAME=$1
    samtools view -h /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.bam | \
    awk '{
        abs=($9<0?-$9:$9)
        if( $1 ~ /^\@/ )
            print
        else if( abs > 110)
            print
    }' >> \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.high-fragments.sam

    sed -e '1 s/SO:coordinate/SO:unsorted/' \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.high-fragments.sam \
    > /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.high-fragments.tmp.sam

    samtools view -S -b /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.high-fragments.tmp.sam \
    > /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.high-fragments.bam 

    sambamba sort -t 16 \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.high-fragments.bam \
    > /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.high-fragments.sorted.bam

    sambamba index -t 16 \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.high-fragments.sorted.bam
}
export -f splitHigh

SAMPLES=( 
    K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19 
    K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19 
)
module load parallel

parallel splitHigh ::: ${SAMPLES[*]}


splitLow() {
    NAME=$1
    samtools view -h /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.bam | \
    awk '{
        abs=($9<0?-$9:$9)
        if( $1 ~ /^\@/ )
            print
        else if( abs < 110)
            print
    }' >> \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.low-fragments.sam

    sed -e '1 s/SO:coordinate/SO:unsorted/' \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.low-fragments.sam \
    > /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.low-fragments.tmp.sam

    samtools view -S -b /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.low-fragments.tmp.sam \
    > /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.low-fragments.bam

    sambamba sort -t 16 \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.low-fragments.bam \
    > /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.low-fragments.sorted.bam

    sambamba index -t 16 \
    /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.nodups.low-fragments.sorted.bam
}
export -f splitLow


SAMPLES=( 
    K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19 
    K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19 
)
module load parallel

parallel splitLow ::: ${SAMPLES[*]}


# call broad peaks for samples
macs2 callpeak \
-t /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.trimmed.bowtie2.shifted.dups.bam \
-c /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/K562_500K_CM_IGG_nan_1ULTN5_PE_1_1_hg19.trimmed.bowtie2.shifted.dups.bam \
--broad --nomodel --extsize 73 --pvalue 1e-3 \
-g hs -n K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19-broad --outdir /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19

macs2 callpeak \
-t /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19.trimmed.bowtie2.shifted.dups.bam \
-c /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/K562_500K_CM_IGG_nan_1ULTN5_PE_1_1_hg19.trimmed.bowtie2.shifted.dups.bam \
--broad --nomodel --extsize 73 --pvalue 1e-3 \
-g hs -n K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19-broad --outdir /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19


# Merge unshifted samples
sambamba merge -t 16 \
/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19-K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.trimmed.bowtie2.nodups.bam \
/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.trimmed.bowtie2.nodups.bam \
/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19.trimmed.bowtie2.nodups.bam

# Merge peaks
cat \
/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19_peaks.narrowPeak \
/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19_peaks.narrowPeak \
| bedtools sort | bedtools merge -i - > /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19-K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19_peaks.narrowPeak 

# Merge broad peaks
cat \
/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19-broad_peaks.broadPeak \
/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19/K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19-broad_peaks.broadPeak \
| bedtools sort | bedtools merge -i - > /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19-K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19_peaks.broadPeak 
