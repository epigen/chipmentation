#!/bin/bash

#parameters

# Take the first argument supplied to the script, remove trailing
# slashes and assign to the samples shell parameter (i.e. variable).

declare samples="${1%%/}";

# For each file that ends with *_sorted.bam in the samples directory repeat the
# following steps.


for file in ${samples}/*_sorted.bam; do


# makeTagDir and UCSC files
declare TagDirin=$(basename "${file}")
declare BSFless=${TagDirin#*#}
declare BSFless2=${BSFless%%_sorted.bam}

#echo "... creating tagDirectory for $TagDirin"
#echo $TagDirin
#echo $BSFless
#echo $BSFless2
makeTagDirectory $BSFless2 $TagDirin >> ${TagDirin}.log 2>&1
mv ${TagDirin}.log $BSFless2


echo "... $TagDirin finished"

done







