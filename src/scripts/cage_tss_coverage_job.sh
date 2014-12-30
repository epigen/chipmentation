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
SAMPLE=$1
CAGE=$2
TECH=$3

### Specify paths
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss/

### Start work

# Get 5' position of signal
if [[ $TECH == CM ]]
	then
	echo "Getting 5' read positions for sample: " ${SAMPLE}_${TECH}
	sbatch /home/arendeiro/jobScripts/bamTo5primeBed.job.sh \
	$PROJECTDIR/mapped/merged/${SAMPLE}_${TECH}.bam \
	$PROJECTDIR/mapped/merged/${SAMPLE}_${TECH}.5prime.bed
elif [[ $TECH == IgG ]]
	then
	echo "Getting 5' read positions for sample: " ${CONTROL}_${TECH}
	sbatch /home/arendeiro/jobScripts/bamTo5primeBed.job.sh \
	$PROJECTDIR/mapped/merged/IgG_K562_500k_CM.bam \
	$PROJECTDIR/mapped/merged/IgG_K562_500k_CM.5prime.bed
elif [[ $TECH == ChIP ]]
	then
	echo "Getting 5' read positions for sample: " {SAMPLE}_${TECH}
	sbatch /home/arendeiro/jobScripts/bamTo5primeBed.job.sh \
	$PROJECTDIR/mapped/merged/${SAMPLE}_${TECH}.bam \
	$PROJECTDIR/mapped/merged/${SAMPLE}_${TECH}.5prime.bed
elif [[ $TECH == DNase ]]
	then
	echo "Getting 5' read positions for DNase sample."
	sbatch /home/arendeiro/jobScripts/bamTo5primeBed.job.sh \
	$PROJECTDIR/mapped/DNase_UWashington_K562_mergedReplicates.bam \
	$PROJECTDIR/mapped/merged/DNase_UWashington_K562_mergedReplicates.5prime.bed
fi

# Calculate coverage
if [[ $TECH == CM ]]
	then
	echo "Calculating" $TECH " coverage in TSS set" $CAGE
	bedtools coverage -d \
	-a /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_CM.5prime.bed \
	-b $CAGEDIR/$CAGE.120bpSlop.bed \
	> $CAGEDIR/$CAGE.120bpSlop.CMcoverage.bed
elif [[ $TECH == IgG ]]
	then
	echo "Calculating" $TECH "coverage in TSS set" $CAGE
	bedtools coverage -d \
	-a /home/arendeiro/data/human/chipmentation/bed/IgG_K562_500k_CM.5prime.bed \
	-b $CAGEDIR/$CAGE.120bpSlop.bed \
	> $CAGEDIR/$CAGE.120bpSlop.IgGcoverage.bed
elif [[ $TECH == ChIP ]]
	then
	echo "Calculating" $TECH "coverage in TSS set" $CAGE
	bedtools coverage -d \
	-a /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_ChIP.5prime.bed \
	-b $CAGEDIR/$CAGE.120bpSlop.bed \
	> $CAGEDIR/$CAGE.120bpSlop.ChIPcoverage.bed
elif [[ $TECH == DNase ]]
	then
	echo "Calculating" $TECH "coverage in TSS set" $CAGE
	bedtools coverage -d \
	-a /home/arendeiro/data/human/chipmentation/bed/DNase_UWashington_K562_mergedReplicates.5prime.bed \
	-b $CAGEDIR/$CAGE.120bpSlop.bed \
	> $CAGEDIR/$CAGE.120bpSlop.DNasecoverage.bed
fi

date
