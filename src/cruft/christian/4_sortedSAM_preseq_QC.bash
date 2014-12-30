#!/bin/bash

# Set-up the environment and define NGS_PICARD
# module load picard-tools

NGS_PICARD="/Users/cschmidl/homer_access/picard-tools-1.115/picard-tools-1.115/"

# Take the first argument supplied to the script, remove trailing
# slashes and assign to the samples shell parameter (i.e. variable).

declare samples="${1%%/}";
declare preseqin=$(basename "${file}")

# For each file that ends with .bam in the samples directory repeat the
# following steps.

for file in ${samples}/*_sorted.bam; do


# preseq c_curve
declare preseqin=$(basename "${file}")
echo "preparing library complexity reports for $preseqin"

preseq c_curve -o qc_c_curve_"${preseqin%%.bam}.txt" \
-B "${preseqin}"


# preseq lc_extrapol

preseq lc_extrap -e 1e8 -s 2e6 -o qc_lc_extrap_"${preseqin%%.bam}.txt" \
-B "${preseqin}"



done














