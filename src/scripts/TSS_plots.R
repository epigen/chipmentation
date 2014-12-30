# This plots average profile plots around TSSs


cageDir = "/fhgfs/groups/lab_bock/shared/data/cage_tss/"
cage = "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed"

cm = read.csv(paste(cageDir, cage, ".120bpSlop.CMcoverage.csv", sep = ''))
rownames(cm) = cm[,1]
cm = cm[,-1]
colnames(cm) = seq(-60, 60)

igg = read.csv(paste(cageDir, cage, ".120bpSlop.IgGcoverage.csv", sep = ''))
rownames(igg) = igg[,1]
igg = igg[,-1]
colnames(igg) = seq(-60, 60)

chip = read.csv(paste(cageDir, cage, ".120bpSlop.ChIPcoverage.csv", sep = ''))
rownames(chip) = chip[,1]
chip = chip[,-1]
colnames(chip) = seq(-60, 60)

dnase = read.csv(paste(cageDir, cage, ".120bpSlop.DNasecoverage.csv", sep = ''))
rownames(dnase) = dnase[,1]
dnase = dnase[,-1]
colnames(dnase) = seq(-60, 60)

# Normalize per library size
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
	stat_smooth(method = "gam", formula = y ~ s(x, k = 240), se = FALSE) + 
	facet_wrap( ~ Var2) + # , ncol = 2, scales = "free"
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank()) #+
	#scale_size(range=c(0.0, 1.2), guide=FALSE)

ggsave(filename = "tss_signal.pdf", plot = p, height = 2, width = 7)
ggsave(filename = "tss_signal_wide.pdf", plot = p, height = 2, width = 4)

p = ggplot(df, aes(Var1, value, colour = Var2)) +
	#geom_line() +
	stat_smooth(method = "gam", formula = y ~ s(x, k = 240), se = FALSE) + 
	#facet_wrap( ~ Var2) + # , ncol = 2, scales = "free"
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank()) #+
	#scale_size(range=c(0.0, 1.2), guide=FALSE)

ggsave(filename = "tss_signal_allinone_wide.pdf", plot = p, height = 2, width = 7)
ggsave(filename = "tss_signal_allinone.pdf", plot = p, height = 2, width = 4)

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
	theme_bw() +
	theme(legend.title=element_blank()) 

ggsave(filename = "tss_signal_1kb.pdf", plot = p, height = 3, width = 6)


window = -60:60
a = df[df$Var1 %in% window, ]

library(grid)
p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5) + 
	#facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	geom_segment(x = -28, xend = -28, y = 0.4, yend = 0.61, colour = "black", size = 0.8, arrow = arrow(length=unit(0.3,"cm"))) + 
	geom_segment(x = 0, xend = 0, y = 0.4, yend = 0.61, colour = "black", size = 0.8, arrow = arrow(length=unit(0.3,"cm"))) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_highlight_arrows.pdf", plot = p, height = 4, width = 4)
ggsave(filename = "tss_signal_120bp_wide_highlight_arrows.pdf", plot = p, height = 2.5, width = 6)


p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5) + 
	facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	#geom_segment(x = -35, xend = -35, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -20, xend = -20, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_facet_highlight.pdf", plot = p, height = 4, width = 4)
ggsave(filename = "tss_signal_120bp_facet_wide_highlight.pdf", plot = p, height = 2.5, width = 6)




# INDEPENDENT


df = cbind(ChIPmentation)
df = melt(df)
window = -60:60
a = df[df$Var1 %in% window, ]

p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5, colour = colors[1]) + 
	#facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	geom_segment(x = -35, xend = -35, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	geom_segment(x = -20, xend = -20, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_wide_CM.pdf", plot = p, height = 2, width = 6)


df = cbind(ChIP)
df = melt(df)
window = -60:60
a = df[df$Var1 %in% window, ]

p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5, colour = colors[2]) + 
	#facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	#geom_segment(x = -35, xend = -35, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -20, xend = -20, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_wide_ChIP.pdf", plot = p, height = 2, width = 6)



df = cbind(DNase)
df = melt(df)
window = -60:60
a = df[df$Var1 %in% window, ]

p = ggplot(a, aes(Var1, value)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5, colour = colors[4]) + 
	#facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	#geom_segment(x = -35, xend = -35, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -20, xend = -20, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_wide_DNase.pdf", plot = p, height = 2, width = 6)


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

