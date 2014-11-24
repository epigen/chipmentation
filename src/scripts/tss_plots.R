R

sample = 'H3K4me3_K562_500k_CM'

#sample = "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed"


names = c(paste("X.", c(2000:1), sep = ""), paste("X", c(0:2000), sep = ""))

cm = read.csv(paste(sample, "_tss_coverage.csv", sep = ''))
rownames(cm) = cm[,1]
cm = cm[,-1]
cm = cm[, names]
colnames(cm) = seq(-2000, 1999)

chip = read.csv(paste(sample, "_tss_coverage_ChIP.csv", sep = ''))
rownames(chip) = chip[,1]
chip = chip[,-1]
chip = chip[, names]
colnames(chip) = seq(-2000, 1999)
igg = read.csv(paste(sample, "_tss_coverage_IgG.csv", sep = ''))
rownames(igg) = igg[,1]
igg = igg[,-1]
igg = igg[, names]
colnames(igg) = seq(-2000, 1999)
dnase = read.csv(paste(sample, "_tss_coverage_DNase.csv", sep = ''))
rownames(dnase) = dnase[,1]
dnase = dnase[,-1]
dnase = dnase[, names]
colnames(dnase) = seq(-2000, 1999)

ChIPmentation = colMeans(cm) / (64054226 / 64054226.)
ChIP = colMeans(chip) / (76894657 / 64054226.)
IgG = colMeans(igg) / (66863728 / 64054226.)
DNase = colMeans(dnase) / (66663835 / 64054226.)

colors = c("#3FAB35", "#4169E1", "#696969", "#8B0000")

### Plot average profiles
library(ggplot2)
library("reshape2")
df = cbind(ChIPmentation, ChIP, IgG, DNase)
df = melt(df)

p = ggplot(df, aes(Var1, value, colour = Var2)) +
	#geom_line() +
	stat_smooth(method = "gam", formula = y ~ s(x, k = 240)) + 
	facet_wrap( ~ Var2) + # , ncol = 2, scales = "free"
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() #+
	#scale_size(range=c(0.0, 1.2), guide=FALSE)

ggsave(filename = "tss_signal.pdf", plot = p, height = 3, width = 6)



window = -1000:1000
a = df[df$Var1 %in% window, ]
p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 240), se = FALSE) + 
	facet_grid(Var2 ~ ., scales = "free") + 
	#coord_cartesian(xlim = c(-1000, 1000)) +
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw()

ggsave(filename = "tss_signal_1kb.pdf", plot = p, height = 3, width = 6)


window = -400:400
a = df[df$Var1 %in% window, ]

p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 550), se = FALSE) + 
	facet_grid(Var2 ~ ., scales = "free") + 
	#coord_cartesian(xlim = c(-200, 200)) +
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw()

ggsave(filename = "tss_signal_400bp.pdf", plot = p, height = 3, width = 6)



### Plot heatmap of CM
library(gplots)
require(made4)

pdf("tss_heatmap_CM.pdf")
heatplot(cm, dend = 'row', labRow = NA, , labCol = NA)
dev.off()

pdf("tss_heatmap_ChIP.pdf")
heatplot(chip, dend = 'row', labRow = NA, , labCol = NA)
dev.off()


window = seq(1000,3000)

pdf("tss_heatmap_1kb_CM.pdf")
heatplot(cm[,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()

pdf("tss_heatmap_1kb_ChIP.pdf")
heatplot(chip[,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()


window = seq(1600,2400)

pdf("tss_heatmap_400bp_CM.pdf")
heatplot(cm[,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()

pdf("tss_heatmap_400bp_ChIP.pdf")
heatplot(chip[,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()

