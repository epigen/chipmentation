# histones
wget -r http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/ \
--no-parent --no-directories \
--accept-regex wgEncodeUwHistoneK562H3k[4me3,27me3]Std.*


module load samtools

cd data/human/encode/chip-seq

# PU1
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep1.bam
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep2.bam

samtools sort wgEncodeHaibTfbsK562Pu1cPcr1xAlnRep1V2.bam wgEncodeHaibTfbsK562Pu1cPcr1xAlnRep1V2.sorted
samtools sort wgEncodeHaibTfbsK562Pu1cPcr1xAlnRep2V2.bam wgEncodeHaibTfbsK562Pu1cPcr1xAlnRep2V2.sorted

sambamba merge -t 12 wgEncodeHaibTfbsK562Pu1cPcr1xAln.bam wgEncodeHaibTfbsK562Pu1cPcr1xAlnRep1V2.sorted.bam wgEncodeHaibTfbsK562Pu1cPcr1xAlnRep2V2.sorted.bam


# CTCF
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC//wgEncodeHaibTfbs/wgEncodeHaibTfbsK562CtcfcPcr1xAlnRep1V2.bam
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC//wgEncodeHaibTfbs/wgEncodeHaibTfbsK562CtcfcPcr1xAlnRep2V2.bam

samtools sort wgEncodeHaibTfbsK562CtcfcPcr1xAlnRep1V2.bam wgEncodeHaibTfbsK562CtcfcPcr1xAlnRep1V2.sorted
samtools sort wgEncodeHaibTfbsK562CtcfcPcr1xAlnRep2V2.bam wgEncodeHaibTfbsK562CtcfcPcr1xAlnRep2V2.sorted

sambamba merge -t 12 wgEncodeHaibTfbsK562CtcfcPcr1xAln.bam wgEncodeHaibTfbsK562CtcfcPcr1xAlnRep1V2.sorted.bam wgEncodeHaibTfbsK562CtcfcPcr1xAlnRep2V2.sorted.bam
