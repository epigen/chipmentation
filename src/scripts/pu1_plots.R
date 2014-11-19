


cm = read.csv('pu1_CM.csv')
rownames(cm) = cm[,1]
colnames(cm) = seq(-2000, 1999)
cm = cm[,-1]
chip = read.csv('pu1_ChIP.csv')
rownames(chip) = chip[,1]
colnames(chip) = seq(-2000, 1999)
chip = chip[,-1]
igg = read.csv('pu1_IgG.csv')
rownames(igg) = igg[,1]
colnames(igg) = seq(-2000, 1999)
igg = igg[,-1]
dnase = read.csv('pu1_DNase.csv')
rownames(dnase) = dnase[,1]
colnames(dnase) = seq(-2000, 1999)
dnase = dnase[,-1]

encode = read.csv('PU1_K562_10mio_CM_peak_coverage_Encode.csv')
rownames(encode) = encode[,1]
colnames(encode) = seq(-2000, 1999)
encode = encode[,-1]

ChIPmentation = colMeans(cm) / (143376696 / 29824946.)
ChIP = colMeans(chip) / (30931217 / 29824946.)
IgG = colMeans(igg) / (29824946 / 29824946.)
DNase = colMeans(dnase) / (66663835 / 29824946.)
ChIP_Encode = colMeans(encode) / (34072563 / 29824946.)

colors = c("#3FAB35", "#4169E1", "#233D8C", "#696969", "#8B0000")

### Plot average profiles
library(ggplot2)
library("reshape2")
df = cbind(ChIPmentation, ChIP, ChIP_Encode, IgG, DNase)
df = melt(df)

p = ggplot(df, aes(Var1, value, colour = Var2)) +
	geom_line() + 
	facet_grid(Var2 ~ ., scales = "free") + 
	coord_cartesian(xlim = c(-1000, 1000)) +
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw()

ggsave(filename = "PU1_peaks_signal_2kb.pdf", plot = p, height = 3, width = 6)


p = ggplot(df, aes(Var1, value, colour = Var2)) +
	geom_line() + 
	facet_grid(Var2 ~ ., scales = "free") + 
	coord_cartesian(xlim = c(-200, 200)) +
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw()

ggsave(filename = "PU1_peaks_signal_400bp.pdf", plot = p, height = 3, width = 6)



### Plot heatmap of CM (footprint)
library(gplots)
require(made4)

window = seq(1600,2400)

r = sample(1:nrow(cm), 2000)
pdf("PU1_peaks_heatmap_400bp_CM.pdf")
heatplot(cm[r,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()

r = sample(1:nrow(chip), 2000)
pdf("PU1_peaks_heatmap_400bp_ChIP.pdf")
heatplot(chip[r,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()

# Plot it centered on motifs
pu1 = read.csv("/home/arendeiro/PU1_K562_10mio_CM_peaks.csv")
pu1 = t(pu1)

r = sample(1:nrow(cm), 2000)
pdf("PU1_peaks_heatmap_motifs_400bp_CM.pdf")
heatplot(pu1[r,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()





### GATA

gata = read.csv("/home/arendeiro/PU1_K562_10mio_CM_peaks.GATAmotif.csv")
gata = t(gata)

window = c(1600:2400)

library(CENTIPEDE)
centFitGATA <- fitCentipede(Xlist = list(DNase=as.matrix(gata[,1850:2250])))

g = gata[order(centFitGATA$PostPr, decreasing = TRUE), ]
r = sample(1:nrow(g), 2000)

pdf("PU1_peaks_heatmap_400bp_GATA_CM.pdf")
heatplot(g[r,window], dend = 'none', labRow = NA, , labCol = NA)
dev.off()

t = c(1:500, (nrow(g) - 500):nrow(g))
col = c(rep(1, 500), rep(2, 501))

pdf("PU1_peaks_heatmap_400bp_GATA_CM_500top_bottom.pdf")
heatplot(g[t,window], dend = 'none', labRow = NA, , labCol = NA, classvec = col)
dev.off()
