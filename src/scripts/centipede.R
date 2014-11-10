
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

peaks = read.table("/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.narrowPeak", sep = "\t", stringsAsFactors=FALSE)
colnames(peaks) = c('chr','start','end','id')
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
pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.400bpAll.pdf")
imageCutSites(pu1[ , 1850:2250][order(centFit$PostPr),])
dev.off()
pdf("/home/arendeiro/projects/chipmentation/results/plots/PU1_K562_10mio_CM_peaks.footprintsCentipede.1000bpAll.pdf")
imageCutSites(pu1[ , 1500:2500][order(centFit$PostPr),])
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

post = cbind(rownames(pu1), centFit$PostPr)
colnames(post) = c('id', 'posterior')
exp = merge(peaks, post)
exp = exp[, c('chr', 'start', 'end', 'id', 'posterior')]
exp$posterior = as.numeric(as.character(exp$posterior))


write.table(exp,
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')


# filter for posterior > 0.9
write.table(exp[exp$posterior > 0.9,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior>09.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')

# filter for posterior < 0.4
write.table(exp[exp$posterior < 0.4,],
    "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.centipedePosterior<04.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep= '\t')

