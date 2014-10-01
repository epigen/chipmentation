#!/bin/bash

# This is the main script for next-generation-sequencing analysis for Tagmentation based ChIP-seq or ATAC-seq

#parameters

declare samples="${1%%/}";


# converting demultiplexed bam files to fastq file with "1_demultiplex_to_fastq_forBOWTIE.bash"
#1_demultiplex_to_fastq_forBOWTIE.bash ${1%%/};


# map fastq file with bowtie to get sam file with "2_mapFASTQ_with_BOWTIE_toSAM.bash"
#2_mapFASTQ_with_BOWTIE_toSAM.bash ${1%%/};

# concert SAM to sorted BAM file with "3_SAM_to_sorted_BAM.bash"
#3_SAM_to_sorted_BAM.bash ${1%%/};

# Quality Control with preseq using the sorted bam file with "4_sortedSAM_preseq_QC.bash"
#4_sortedSAM_preseq_QC.bash ${1%%/};

# manipulate sequencing files so that position of transposition reaction can be identified correctly with "5_AWK_sortedBAM_to_SAM.bash"
5_AWK_sortedBAM_to_SAM.bash ${1%%/};

# create TagDirectories with "6_TagDirectory.bash"
6_TagDirectory.bash ${1%%/};

# create UCSCfiles 1.normal 2. footprint resolution with "7_makeUCSCfile.bash"
7_makeUCSCfile.bash ${1%%/};

