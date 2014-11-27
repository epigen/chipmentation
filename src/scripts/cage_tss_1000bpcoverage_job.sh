#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --time=48:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32000
#SBATCH --nodes=1

#SBATCH --job-name=cagePeaksCoverage
#SBATCH --output=/home/arendeiro/logs/cagePeaksCoverage_%j.out

# *** setup environment ***
# load the required environmental modules that your script depends upon

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
CAGE=$1
TECH=$2
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss/

### Specify paths

### Start work
# echo "Calculating ChIPmentation coverage in TSS set" $CAGE
# bedtools coverage -d \
# -a /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_CM.bed \
# -b $CAGEDIR/$CAGE.1000bpSlop.bed \
# > $CAGEDIR/$CAGE.1000bpSlop.CMcoverage.bed
if [[ $TECH == CM ]]
	then
	echo "Calculating" $TECH " coverage in TSS set" $CAGE
	bedtools coverage -d \
	-a /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_CM.bed \
	-b $CAGEDIR/$CAGE.1000bpSlop.bed \
	> $CAGEDIR/$CAGE.1000bpSlop.CMcoverage.bed
elif [[ $TECH == IgG ]]
	then
	echo "Calculating" $TECH "coverage in TSS set" $CAGE
	bedtools coverage -d \
	-a /home/arendeiro/data/human/chipmentation/bed/IgG_K562_500k_CM.bed \
	-b $CAGEDIR/$CAGE.1000bpSlop.bed \
	> $CAGEDIR/$CAGE.1000bpSlop.IgGcoverage.bed
elif [[ $TECH == ChIP ]]
	then
	echo "Calculating" $TECH "coverage in TSS set" $CAGE
	bedtools coverage -d \
	-a /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_ChIP.bed \
	-b $CAGEDIR/$CAGE.1000bpSlop.bed \
	> $CAGEDIR/$CAGE.1000bpSlop.ChIPcoverage.bed
elif [[ $TECH == DNase ]]
	then
	echo "Calculating" $TECH "coverage in TSS set" $CAGE
	bedtools coverage -d \
	-a /home/arendeiro/data/human/chipmentation/bed/DNase_UWashington_K562_mergedReplicates.bed \
	-b $CAGEDIR/$CAGE.1000bpSlop.bed \
	> $CAGEDIR/$CAGE.1000bpSlop.DNasecoverage.bed
fi

date
