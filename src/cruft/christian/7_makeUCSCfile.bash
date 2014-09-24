#!/bin/bash

#parameters
colorRGB="0,0,170"
colorRGB2="0,0,100"

# Take the first argument supplied to the script, remove trailing
# slashes and assign to the samples shell parameter (i.e. variable).

declare samples="${1%%/}";

# For each file that ends with *_sorted.bam in the samples directory repeat the
# following steps.


for directory in ${samples}/*_footprint/; do


# make UCSC files
declare TagDir=$(basename "${directory}")
declare TagDir_less=${TagDir%%_footprint}

makeUCSCfile ${TagDir} -tbp 1 -fragLength 1 -color $colorRGB -name $TagDir -o $TagDir/${TagDir}.txt

makeUCSCfile ${TagDir} -tbp 1 -color $colorRGB2 -name $TagDir_less -o $TagDir/${TagDir_less}.txt


done







