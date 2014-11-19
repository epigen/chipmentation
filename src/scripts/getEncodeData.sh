wget -r http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwHistone/ \
--no-parent --no-directories \
--accept-regex wgEncodeUwHistoneK562H3k[4me3,27me3]Std.*


# PU1
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep1.bam
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep2.bam

sambamba sort -t 12 wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep1.bam > ~/tmp

sambamba merge -t 12 wgEncodeHaibTfbsK562Pu1Pcr1xAln.bam wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep1.bam wgEncodeHaibTfbsK562Pu1Pcr1xAlnRep2.bam
