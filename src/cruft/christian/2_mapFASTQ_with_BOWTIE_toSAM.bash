#!/bin/bash

# Set-up the environment and define NGS_PICARD
# module load picard-tools

NGS_PICARD="/Users/cschmidl/homer_access/picard-tools-1.115/picard-tools-1.115/"

# Take the first argument supplied to the script, remove trailing
# slashes and assign to the samples shell parameter (i.e. variable).

declare samples="${1%%/}";

for file in ${samples}/*.fastq; do

# For each file that ends with .fastq in the samples directory repeat the following steps.
declare read1=$(basename "${file}")
bowtieout="${read1%%_R1.fastq}.sam"
logout="${read1%%_R1.fastq}.log"
echo "... mapping $read1 with bowtie"

# align with bowtie
bowtie2 -p 3 -x /Users/cschmidl/bowtie2-2.2.1/Bowtie2Index/genome $read1 -S ${bowtieout} >> $logout 2>&1
echo "... mapping of file $read1 done"

#remove fastq file
echo "... removing original .FASTQ file and keeping mapped .SAM file"
rm $read1

done












