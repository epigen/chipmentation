#!/bin/bash

# Set-up the environment and define NGS_PICARD
# module load picard-tools

if [ -z ${NGS_PICARD} ]; then
    NGS_PICARD="/Users/cschmidl/homer_access/picard-tools-1.115/picard-tools-1.115"
fi

# Take the first argument supplied to the script, remove trailing
# slashes and assign to the samples shell parameter (i.e. variable).

declare samples="${1%%/}";

# echo "Samples directory: ${samples}";
# ls -la "${samples}";

# For each file that ends with .bam in the samples directory repeat the
# following steps.

for file in ${samples}/*.bam; do


    # Remove the absolute directory part of the file name to have a
    # relative name.
    declare read1=$(basename "${file}")
    echo "... converting demultiplexed .BAM file to .FASTQ file for $read1"

    # Match away the .bam part of the relative file name and append .fastq.
    read1="${read1%%.bam}_R1.fastq";

    # Same for read2, whether it exists or not.
    declare read2=$(basename "${file}")
    read2="${read2%%.bam}_R2.fastq";

    # Run Picard SamToFastq on the original BAM file with absolute
    # file name pointing to the BSF directory and expand into R1 and R2
    # in the local directory.

    java -Xmx2G -jar ${NGS_PICARD}/SamToFastq.jar \
        INPUT="${file}" \
        FASTQ="${read1}" \
        # SECOND_END_FASTQ="${read2}"


    #remove original BAM file if copied in folder

done












