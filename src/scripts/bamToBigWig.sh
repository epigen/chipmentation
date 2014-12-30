#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Optional parameters
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --nodes=1
#SBATCH --job-name=bedGraphToBigWig
#SBATCH --output=/home/arendeiro/logs/bedGraphToBigWig_%j.out
date

SAMPLE=$1
SAMPLE_NAME=`basename ${SAMPLE%.bam}`

echo $SAMPLE
echo $SAMPLE_NAME

GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromlengths.txt
OUTPUT_FOLDER=/fhgfs/groups/lab_bock/public_html/arendeiro/chipmentation2014/bigWigTracks
URL=http://www.biomedical-sequencing.at/bocklab/arendeiro/chipmentation2014/bigWigTracks/
TRACKHUB=${OUTPUT_FOLDER}/tracks_hub.txt

# Track with full reads
echo "Making bigWig track for sample: " $SAMPLE_NAME
bamToBed -i $SAMPLE | \
genomeCoverageBed -i stdin -bg -g $GENOMESIZE > \
$OUTPUT_FOLDER/${SAMPLE_NAME}.cov

bedGraphToBigWig $OUTPUT_FOLDER/${SAMPLE_NAME}.cov $GENOMESIZE ${OUTPUT_FOLDER}/${SAMPLE}.bw
rm bedGraphToBigWig $OUTPUT_FOLDER/${SAMPLE_NAME}.cov

echo 'track type=bigWig name="'$SAMPLE_NAME'" description="'$SAMPLE_NAME'" visibility=3 bigDataUrl='$URL''$SAMPLE_NAME'.bw' >> $TRACKHUB

# Track with 5' position for ChIPmentation
if [[ $SAMPLE == *_CM.bam ]]
    then
    echo "Making 5 prime position bigWig track for sample: " $SAMPLE_NAME
    bamToBed -i $SAMPLE | \
    python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py | \
	genomeCoverageBed -i stdin -bg -g $GENOMESIZE > \
	$OUTPUT_FOLDER/${SAMPLE_NAME}.5prime.cov
    
    bedGraphToBigWig $OUTPUT_FOLDER/${SAMPLE_NAME}.5prime.cov $GENOMESIZE ${OUTPUT_FOLDER}/${SAMPLE}.5prime.bw
    rm bedGraphToBigWig $OUTPUT_FOLDER/${SAMPLE_NAME}.5prime.cov

    echo 'track type=bigWig name="'$SAMPLE_NAME'" description="'$SAMPLE_NAME'" visibility=3 bigDataUrl='$URL''${SAMPLE_NAME}'.5prime.bw' >> $TRACKHUB
fi

date
