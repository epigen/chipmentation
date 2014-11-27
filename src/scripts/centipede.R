
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

#pu1 = read.csv("/home/arendeiro/PU1_K562_10mio_CM_peaks.csv")
#pu1 = t(pu1)

gata = read.csv("/home/arendeiro/PU1_K562_10mio_CM_peaks.GATAmotif.csv")
gata = t(gata)
#cons = read.table("/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.conservation.bed", sep = "\t", stringsAsFactors=FALSE)

#score = peaks[,8]
# try different peak windows
# get conservation score for whole window
#centFitPU1 <- fitCentipede(Xlist = list(DNase=as.matrix(pu1[,1850:2250])))
centFitGATA <- fitCentipede(Xlist = list(DNase=as.matrix(gata[,1850:2250])))
    #Y=cbind(cons[,5], score)

df = data.frame(postProb = centFitGATA$PostPr, sample = 'GATA motifs')

p = ggplot(df, aes(sample, postProb)) +
    geom_boxplot() +
    theme_bw()
ggsave(filename = "/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bp.PU1-GATA.boxplot.pdf",
    plot = p, height = 4, width = 1.25)

p = ggplot(df, aes(sample, postProb)) +
    geom_violin(alpha=0.5, color="gray") + 
    theme_bw()
ggsave(filename = "/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bp.PU1-GATA.violinplot.pdf",
    plot = p, height = 4, width = 1.25)

library(ggplot2)
# pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bpAll.pu1.pdf")
# imageCutSites(pu1[ , 1850:2250][order(centFitPU1$PostPr),])
# dev.off()
# pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.1000bpAll.pu1.pdf")
# imageCutSites(pu1[ , 1500:2500][order(centFitPU1$PostPr),])
# dev.off()
# pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bp90.pu1.pdf")
# imageCutSites(pu1[centFit$PostPr > 0.9 , 1850:2250])
# dev.off()
# pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.1000bp90.pu1.pdf")
# imageCutSites(pu1[centFit$PostPr > 0.9 , 1500:2500])
# dev.off()

pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bpAll.gata.pdf")
imageCutSites(gata[ , 1850:2250][order(centFitGATA$PostPr),])
dev.off()
pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.1000bpAll.gata.pdf")
imageCutSites(gata[ , 1500:2500][order(centFitGATA$PostPr),])
dev.off()
pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bp90.gata.pdf")
imageCutSites(gata[centFitGATA$PostPr > 0.9 , 1850:2250])
dev.off()
pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.1000bp90.gata.pdf")
imageCutSites(gata[centFitGATA$PostPr > 0.9 , 1500:2500])
dev.off()



# Export list of peaks
# peaksPU1 = read.table("/home/arendeiro/data/human/chipmentation/bed/PU1_K562_10mio_CM.bed", sep = "\t", stringsAsFactors=FALSE)
# colnames(peaksPU1) = c('chr','start','end','id')

peaksGATA = read.table("/home/arendeiro/data/human/chipmentation/bed/PU1_K562_10mio_CM.GATAmotif.bed", sep = "\t", stringsAsFactors=FALSE)
colnames(peaksGATA) = c('chr','start','end','id')

# postPU1 = cbind(rownames(pu1), centFitPU1$PostPr)
# colnames(postPU1) = c('id', 'posterior')
# expPU1 = merge(peaksPU1, postPU1)
# expPU1 = expPU1[, c('chr', 'start', 'end', 'id', 'posterior')]
# expPU1$posterior = as.numeric(as.character(expPU1$posterior))

postGATA = cbind(rownames(gata), centFitGATA$PostPr)
colnames(postGATA) = c('id', 'posterior')
expGATA = merge(peaksGATA, postGATA)
expGATA = expGATA[, c('chr', 'start', 'end', 'id', 'posterior')]
expGATA$posterior = as.numeric(as.character(expGATA$posterior))

# write.table(expPU1,
#     "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.bed",
#     quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# # filter for posterior > 0.9
# write.table(expPU1[expPU1$posterior > 0.9,],
#     "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior>09.bed",
#     quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# # filter for posterior < 0.4
# write.table(expPU1[expPU1$posterior < 0.4,],
#     "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior<04.bed",
#     quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')

write.table(expGATA,
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# filter for posterior > 0.99
write.table(expGATA[expGATA$posterior > 0.99,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif>099.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# filter for posterior > 0.9
write.table(expGATA[expGATA$posterior > 0.9,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif>09.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# filter for posterior > 0.8
write.table(expGATA[expGATA$posterior > 0.8,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif>08.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# filter for posterior > 0.5
write.table(expGATA[expGATA$posterior > 0.5,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif>05.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')


# filter for posterior < 0.8
write.table(expGATA[expGATA$posterior < 0.8,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif<08.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')

# filter for posterior < 0.4
write.table(expGATA[expGATA$posterior < 0.4,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif<04.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')
# filter for posterior < 0.4
write.table(expGATA[expGATA$posterior < 0.3,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif<03.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')

# filter for posterior < 0.4
write.table(expGATA[expGATA$posterior < 0.1,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif<01.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')




PROJECTDIR=/home/arendeiro/data/human/chipmentation
GATADIR=/home/arendeiro/data/human/encode/chip-seq

for T in 099 09 08
do
for GATA in wgEncodeAwgTfbsHaibK562Gata2sc267Pcr1xUniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gata1UcdUniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gata2UcdUniPk.narrowPeak wgEncodeAwgTfbsUchicagoK562Egata2UniPk.narrowPeak
do
# Have footprint
## have peak
yy=`bedtools intersect -u \
-a "$PROJECTDIR/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif>${T}.bed" \
-b $GATADIR/$GATA \
| wc -l`
## don't have peak
yn=`bedtools intersect -v \
-a "$PROJECTDIR/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif>${T}.bed" \
-b $GATADIR/$GATA \
| wc -l`

# Don't have footprint
## have peak
ny=`bedtools intersect -u \
-a "$PROJECTDIR/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif<08.bed" \
-b $GATADIR/$GATA \
| wc -l`
## don't have peak (subtract this from above)
nn=`bedtools intersect -v \
-a "$PROJECTDIR/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.GATAmotif<08.bed" \
-b $GATADIR/$GATA \
| wc -l`

echo "fisher_exact([[$yy, $yn], [$ny, $nn]])"

done
done

# < 0.8
# 0.99
fisher_exact([[371, 935], [2209, 6947]])
(1.2478543053091025, 0.00097949588576523258)
fisher_exact([[130, 1176], [753, 8403]])
(1.2336030029541698, 0.037943902478709074)
fisher_exact([[256, 1050], [1498, 7658]])
(1.2463907432131731, 0.0038385225223134039)
fisher_exact([[248, 1058], [1551, 7605]])
(1.1493529998939644, 0.071308157380166964)

# 0.9
fisher_exact([[521, 1350], [2209, 6947]])
(1.2136837516556847, 0.00080168127069898541)
fisher_exact([[182, 1689], [753, 8403]])
(1.2024890373379189, 0.036122506975226114)
fisher_exact([[362, 1509], [1498, 7658]])
(1.2263738441624397, 0.0020373398490278897)
fisher_exact([[344, 1527], [1551, 7605]])
(1.104604545644549, 0.13037068871019813)

# 0.8
fisher_exact([[613, 1600], [2209, 6947]])
(1.2048752263467633, 0.00054872924704832827)
fisher_exact([[213, 2000], [753, 8403]])
(1.1884721115537848, 0.037311978727615686)
fisher_exact([[418, 1795], [1498, 7658]])
(1.1904615624918646, 0.0048515516023193072)
fisher_exact([[398, 1815], [1551, 7605]])
(1.0752114071966368, 0.24496550187737376)

# < 0.4
# 0.99
fisher_exact([[371, 935], [1943, 5871]])
fisher_exact([[130, 1176], [684, 7130]])
fisher_exact([[256, 1050], [1324, 6490]])
fisher_exact([[248, 1058], [1363, 6451]])

# 0.9
fisher_exact([[521, 1350], [1943, 5871]])
fisher_exact([[182, 1689], [684, 7130]])
fisher_exact([[362, 1509], [1324, 6490]])
fisher_exact([[344, 1527], [1363, 6451]])

# 0.8
fisher_exact([[613, 1600], [1943, 5871]])
fisher_exact([[213, 2000], [684, 7130]])
fisher_exact([[418, 1795], [1324, 6490]])
fisher_exact([[398, 1815], [1363, 6451]])

# 0.5
fisher_exact([[809, 2374], [1943, 5871]])
fisher_exact([[258, 2925], [684, 7130]])
fisher_exact([[548, 2635], [1324, 6490]])
fisher_exact([[533, 2650], [1363, 6451]])



