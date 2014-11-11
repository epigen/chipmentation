
library(DNaseR)

chrmSizes = read.table("/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt",sep="\t")
names(chrmSizes) = c("chr", "size")

footprints(bam = "data/human/chipmentation/mapped/merged/PU1_K562_10mio_CM.bam",
 chrN = chrmSizes$chr, chrL = chrmSizes$size,
 p = 1e-4, width = c(6,40), N = 2e6)


python

import sys
import csv
import pandas as pd

sample = 'PU1_K562_10mio_CM'
projectDir = '/home/arendeiro/data/human/chipmentation/'

def parseBedCoverage(bedFile):
    """
    Reads tab-delimiter file with two rows describing the chromossomes and its lengths.
    Returns dictionary of chr:sizes.
    """
    prev_peak = ''
    peaks = {}
    with open(bedFile, "rb") as f:
        reader = csv.reader(f, delimiter="\t")
        for line, row in enumerate(reader):
            peak = str(row[3])
            bp = int(row[4])
            count = int(row[5])
            if peak != prev_peak:
                # new peak
                peaks[peak] = {}
                peaks[peak][bp] = count
                prev_peak = peak
            else:
                peaks[peak][bp] = count
    sortedPeaks = {}
    for key in sorted(peaks.iterkeys()):
        sortedPeaks[key] = peaks[key]
    return sortedPeaks

covCM = parseBedCoverage(projectDir + 'bed/' + sample + '_peak_coverage.GATAmotif.bed')
covCMDF = pd.DataFrame(covCM)
a = covCMDF.T
covCMDF.to_csv('/home/arendeiro/PU1_K562_10mio_CM_peaks.GATAmotif.csv', index = False, encoding='utf-8')



#
sh

cut -f 1,2,3,4 /home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.narrowPeak > tmp

bigWigAverageOverBed /home/arendeiro/reference/Homo_sapiens/hg38.phastCons7way.bw \
tmp \
/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.conservation.tab \
-bedOut=/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.conservation.bed


R
library(CENTIPEDE)

pu1 = read.csv("/home/arendeiro/PU1_K562_10mio_CM_peaks.csv")
pu1 = t(pu1)

gata = read.csv("/home/arendeiro/PU1_K562_10mio_CM_peaks.GATAmotif.csv")
gata = t(gata)
cons = read.table("/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.conservation.bed", sep = "\t", stringsAsFactors=FALSE)

score = peaks[,8]
# try different peak windows
# get conservation score for whole window
centFitPU1 <- fitCentipede(Xlist = list(DNase=as.matrix(pu1[,1850:2250])))
centFitGATA <- fitCentipede(Xlist = list(DNase=as.matrix(gata[,1850:2250])))
    #Y=cbind(cons[,5], score)

df = data.frame(postProb = centFitPU1$PostPr, sample = 'PU1 motifs')
df = rbind(df, data.frame(postProb = centFitGATA$PostPr, sample = 'GATA motifs'))

pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bp.PU1vsGATA.pdf")
ggplot(df, aes(sample, postProb)) +
    geom_boxplot() +
    theme_bw()
dev.off()


library(ggplot2)
pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bpAll.gata.pdf")
imageCutSites(gata[ , 1850:2250][order(centFitGATA$PostPr),])
dev.off()
pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.1000bpAll.gata.pdf")
imageCutSites(gata[ , 1500:2500][order(centFitGATA$PostPr),])
dev.off()
pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bp90.pdf")
imageCutSites(pu1[centFit$PostPr > 0.9 , 1850:2250])
dev.off()
pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.1000bp90.pdf")
imageCutSites(pu1[centFit$PostPr > 0.9 , 1500:2500])
dev.off()



pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.pdf")
par(mfrow = c(2,2))
plotProfile(centFit$LambdaParList[[1]],Mlen=21)
boxplot(centFit$PostPr, notch = TRUE, ylab = 'Posterior Prob')
dev.off()



# Export list of peaks
peaksPU1 = read.table("/home/arendeiro/data/human/chipmentation/bed/PU1_K562_10mio_CM.bed", sep = "\t", stringsAsFactors=FALSE)
colnames(peaksPU1) = c('chr','start','end','id')

peaksGATA = read.table("/home/arendeiro/data/human/chipmentation/bed/PU1_K562_10mio_CM.GATAmotif.bed", sep = "\t", stringsAsFactors=FALSE)
colnames(peaksGATA) = c('chr','start','end','id')



postPU1 = cbind(rownames(pu1), centFitPU1$PostPr)
colnames(postPU1) = c('id', 'posterior')
expPU1 = merge(peaksPU1, postPU1)
expPU1 = expPU1[, c('chr', 'start', 'end', 'id', 'posterior')]
expPU1$posterior = as.numeric(as.character(expPU1$posterior))

postGATA = cbind(rownames(gata), centFitGATA$PostPr)
colnames(postGATA) = c('id', 'posterior')
expGATA = merge(peaksGATA, postGATA)
expGATA = expGATA[, c('chr', 'start', 'end', 'id', 'posterior')]
expGATA$posterior = as.numeric(as.character(expGATA$posterior))



write.table(expPU1,
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# filter for posterior > 0.9
write.table(expPU1[expPU1$posterior > 0.9,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior>09.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# filter for posterior < 0.4
write.table(expPU1[expPU1$posterior < 0.4,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior<04.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')

write.table(expGATA,
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# filter for posterior > 0.9
write.table(expGATA[expGATA$posterior > 0.9,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif>09.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# filter for posterior < 0.4
write.table(expGATA[expGATA$posterior < 0.4,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif<04.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')



PROJECTDIR=/home/arendeiro/data/human/chipmentation
GATADIR=/home/arendeiro/data/human/encode/chip-seq

for GATA in wgEncodeAwgTfbsHaibK562Gata2sc267Pcr1xUniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gata1UcdUniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gata2UcdUniPk.narrowPeak wgEncodeAwgTfbsUchicagoK562Egata2UniPk.narrowPeak
do
echo $GATA
# Have footprint
## have peak
bedtools intersect -u \
-a "$PROJECTDIR/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif>09.bed" \
-b $GATADIR/$GATA \
| wc -l
## don't have peak (subtract this from above)
wc -l /home/arendeiro/data/human/chipmentation/bed/PU1_K562_10mio_CM.GATAmotif.bed

# Don't have footprint
## have peak
bedtools intersect -u \
-a "$PROJECTDIR/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif<04.bed" \
-b $GATADIR/$GATA \
| wc -l
## don't have peak (subtract this from above)
wc -l /home/arendeiro/data/human/chipmentation/bed/PU1_K562_10mio_CM.GATAmotif.bed
done

wgEncodeAwgTfbsHaibK562Gata2sc267Pcr1xUniPk.narrowPeak
521
10848
1943
9426
wgEncodeAwgTfbsSydhK562Gata1UcdUniPk.narrowPeak
182
11187
684
10685
wgEncodeAwgTfbsSydhK562Gata2UcdUniPk.narrowPeak
362
11007
1324
10045
wgEncodeAwgTfbsUchicagoK562Egata2UniPk.narrowPeak
344
11025
1363
10006


fisher_exact([[521,10848],[1943, 9426]])
(0.23299289712560178, 2.239359192110117e-213)
fisher_exact([[182,11187],[684, 10685]])
(0.25414184279267343, 1.3907346750027083e-71)
fisher_exact([[362,11007],[1324, 10045]])
(0.24951781577062879, 3.8063320109304112e-138)
fisher_exact([[344,11025],[1363, 10006]])
(0.22905748457367783, 5.7237439107939745e-154)