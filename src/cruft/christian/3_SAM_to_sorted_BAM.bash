#!/bin/bash

# Set-up the environment and define NGS_PICARD
# module load picard-tools

NGS_PICARD="/Users/cschmidl/homer_access/picard-tools-1.115/picard-tools-1.115/"

# Take the first argument supplied to the script, remove trailing
# slashes and assign to the samples shell parameter (i.e. variable).

declare samples="${1%%/}";

for file in ${samples}/*.sam; do


# For each file that ends with .sam in the samples directory repeat the following steps.
# make a sorted bam file with Picard SortSam
declare Picardin=$(basename "${file}")
echo "... converting mapped .SAM file to sorted .BAM file with PICARD SortSAM for fiel $Picardin"

java -Xmx2G -jar ${NGS_PICARD}SortSam.jar \
INPUT=$Picardin \
OUTPUT="${Picardin%%.sam}_sorted.bam" \
SORT_ORDER=coordinate

# remove SAM file (mapped sequences from bowtie). from now on work with the sorted .BAM file
echo "... remove mapped .SAM file and keep only sorted .BAM file"
rm $Picardin


done












