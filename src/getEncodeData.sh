#!/usr/env bash

getBamFile () {
    SAMPLE=$1
    BAMFOLDER=/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/
    URL=hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/

    rsync -a -P rsync://${URL}/$SAMPLE.bam $BAMFOLDER
    rsync -a -P rsync://${URL}/$SAMPLE.bam.bai $BAMFOLDER
}
export -f getBamFile

mergeBams () {
    SAMPLE=$1
    BAMFOLDER=/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped

    samtools merge -f ${BAMFOLDER}/${SAMPLE}.bam ${BAMFOLDER}/${SAMPLE}Rep1.bam ${BAMFOLDER}/${SAMPLE}Rep2.bam
    samtools index ${BAMFOLDER}/${SAMPLE}.bam
}
export -f mergeBams

module load parallel
parallel getBamFile ::: {wgEncodeBroadHistoneK562H3k4me1StdAln,wgEncodeBroadHistoneK562H3k27me3StdAln,wgEncodeBroadHistoneK562H3k36me3StdAln,wgEncodeBroadHistoneK562H3k4me3StdAln}Rep{1,2}
parallel mergeBams ::: wgEncodeBroadHistoneK562H3k4me1StdAln wgEncodeBroadHistoneK562H3k27me3StdAln wgEncodeBroadHistoneK562H3k36me3StdAln wgEncodeBroadHistoneK562H3k4me3StdAln

# Input
getBamFile wgEncodeBroadHistoneK562ControlStdAlnRep1


# getPeakFile () {
#     SAMPLE=$1
#     BAMFOLDER=/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/
#     URL=hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/

#     rsync -a -P rsync://${URL}/$SAMPLE.bam $BAMFOLDER
#     rsync -a -P rsync://${URL}/$SAMPLE.bam.bai $BAMFOLDER
# }
# export -f getPeakFile

#parallel getPeakFile ::: {wgEncodeBroadHistoneK562H3k4me1StdAln,wgEncodeBroadHistoneK562H3k27me3StdAln,wgEncodeBroadHistoneK562H3k36me3StdAln,wgEncodeBroadHistoneK562H3k4me3StdAln}Rep{1,2}
