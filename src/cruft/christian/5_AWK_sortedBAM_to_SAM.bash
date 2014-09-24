#!/bin/bash

#parameters
NGS_PICARD="/Users/cschmidl/homer_access/picard-tools-1.115/picard-tools-1.115/"

# Take the first argument supplied to the script, remove trailing
# slashes and assign to the samples shell parameter (i.e. variable).

declare samples="${1%%/}";

# For each file that ends with *_sorted.bam in the samples directory repeat the
# following steps.


for file in ${samples}/*_sorted.bam; do

echo "... adjusting position of reads of $input to detect original position of transposase reaction for $AWKin"

# Remove the absolute directory part of the file name to have a
# relative name.
declare AWKin=$(basename "${file}")


# awk to move +/- strand reads to detect the origional tagmentation site

samtools view ${AWKin} | awk -F"\t" 'BEGIN{OFS="\t"}($2==16){ offset=$4-5; print $1,$2,$3, offset,$5,$6,$7,$8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' > ${AWKin%%.bam}_minus.sam

samtools view ${AWKin} | awk -F"\t" 'BEGIN{OFS="\t"}($2==0){ offset=$4+4; print $1,$2,$3, offset,$5,$6,$7,$8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' > ${AWKin%%.bam}_plus.sam

# merge Sam files
echo "... merging + and - strand reads"

cat ${AWKin%%.bam}_minus.sam ${AWKin%%.bam}_plus.sam > ${AWKin%%.bam}_footprint.sam

echo "... $AWKin finished"

# remove *_plus.sam and *_minus.sam files
echo "... removing *_minus.sam and *_plus.sam files and keeping only merged + and - strand read-file"

rm ${AWKin%%.bam}_minus.sam
rm ${AWKin%%.bam}_plus.sam

done







