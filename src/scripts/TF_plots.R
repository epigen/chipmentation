
bedDir = "/home/arendeiro/data/human/chipmentation/bed/"
resultsDir = "/home/arendeiro/projects/chipmentation/results/"
plotsDir = "/home/arendeiro/projects/chipmentation/results/plots/"
sample = "PU1_K562_10mio_CM" # CTCF_K562_10mio_CM

cm = read.csv(paste0(bedDir, sample, "_peak_coverage_CM.csv"))
rownames(cm) = cm[,1]
colnames(cm) = seq(-2000, 1999)
cm = cm[,-1]
chip = read.csv(paste0(bedDir, sample, "_peak_coverage_ChIP.csv"))
rownames(chip) = chip[,1]
colnames(chip) = seq(-2000, 1999)
chip = chip[,-1]
igg = read.csv(paste0(bedDir, sample, "_peak_coverage_IgG.csv"))
rownames(igg) = igg[,1]
colnames(igg) = seq(-2000, 1999)
igg = igg[,-1]
dnase = read.csv(paste0(bedDir, sample, "_peak_coverage_DNase.csv"))
rownames(dnase) = dnase[,1]
colnames(dnase) = seq(-2000, 1999)
dnase = dnase[,-1]
encode = read.csv(paste0(bedDir, sample, "_peak_coverage_Encode.csv"))
rownames(encode) = encode[,1]
colnames(encode) = seq(-2000, 1999)
encode = encode[,-1]

# normalize by total library size
if (sample == "PU1_K562_10mio_CM") {
    ChIPmentation = colMeans(cm) / (143376696 / 29824946.)
    ChIP = colMeans(chip) / (30931217 / 29824946.)
    IgG = colMeans(igg) / (29824946 / 29824946.)
    DNase = colMeans(dnase) / (66663835 / 29824946.)
    ChIP_Encode = colMeans(encode) / (34072563 / 29824946.)
} else if (sample == "CTCF_K562_10mio_CM") {
    ChIPmentation = colMeans(cm) / (25369965 / 25369965.)
    ChIP = colMeans(chip) / (34063407 / 25369965.)
    IgG = colMeans(igg) / (29824946 / 25369965.)
    DNase = colMeans(dnase) / (66663835 / 25369965.)
    ChIP_Encode = colMeans(encode) / (55971147 / 25369965.)
}

colors = c("#3FAB35", "#4169E1", "black", "#696969", "#8B0000")

### Plot average profiles
library(ggplot2)
library("reshape2")
df = cbind(ChIPmentation, ChIP, ChIP_Encode, IgG, DNase)
df = melt(df)

p = ggplot(df, aes(Var1, value, colour = Var2)) +
    geom_line(size = 1.5, alpha = 0.9) + 
    facet_grid(Var2 ~ ., scales = "free") + 
    xlab("Distance to motif") +
    ylab("Tags") +
    scale_color_manual(values = colors) +
    theme_bw()

ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_4kb.strandSpecific.pdf"), plot = p, height = 6, width = 6)
ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_4kb.strandSpecific.wide.pdf"), plot = p, height = 6, width = 8)


p = ggplot(df, aes(Var1, value, colour = Var2)) +
    geom_line(size = 1.5, alpha = 0.9) + 
    facet_grid(Var2 ~ ., scales = "free") + 
    coord_cartesian(xlim = c(-1000, 1000)) +
    xlab("Distance to motif") +
    ylab("Tags") +
    scale_color_manual(values = colors) +
    theme_bw()

ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_2kb.strandSpecific.pdf"), plot = p, height = 6, width = 6)
ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_2kb.strandSpecific.wide.pdf"), plot = p, height = 6, width = 8)

p = ggplot(df, aes(Var1, value, colour = Var2)) +
    geom_line(size = 1.5, alpha = 0.9) + 
    facet_grid(Var2 ~ ., scales = "free") + 
    coord_cartesian(xlim = c(-500, 500)) +
    xlab("Distance to motif") +
    ylab("Tags") +
    scale_color_manual(values = colors) +
    theme_bw()

ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_1kb.strandSpecific.pdf"), plot = p, height = 6, width = 6)
ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_1kb.strandSpecific.wide.pdf"), plot = p, height = 6, width = 8)

p = ggplot(df, aes(Var1, value, colour = Var2)) +
    geom_line(size = 1.5, alpha = 0.9) + 
    facet_grid(Var2 ~ ., scales = "free") + 
    coord_cartesian(xlim = c(-200, 200)) +
    xlab("Distance to motif") +
    ylab("Tags") +
    scale_color_manual(values = colors) +
    theme_bw()

ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_400bp.strandSpecific.pdf"), plot = p, height = 6, width = 6)
ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_400bp.strandSpecific.wide.pdf"), plot = p, height = 6, width = 8)

p = ggplot(df, aes(Var1, value, colour = Var2)) +
    geom_line(size = 1.5, alpha = 0.9) + 
    facet_grid(Var2 ~ ., scales = "free") + 
    coord_cartesian(xlim = c(-100, 100)) +
    xlab("Distance to motif") +
    ylab("Tags") +
    scale_color_manual(values = colors) +
    theme_bw()

ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_200bp.strandSpecific.pdf"), plot = p, height = 6, width = 6)
ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_200bp.strandSpecific.wide.pdf"), plot = p, height = 6, width = 8)



# only CM
df2 = cbind(ChIPmentation[1000:3000])
df2 = melt(df2)
p = ggplot(df2, aes(Var1, value)) +
    geom_line(size = 1.5, color = colors[1]) + 
    xlab("Distance to motif") +
    ylab("Tags") +
    theme_bw()

ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_2kb_CMonly.strandSpecific.pdf"), plot = p, height = 6, width = 6)
ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_2kb_CMonly.strandSpecific.pdf"), plot = p, height = 6, width = 8)

df2 = cbind(ChIPmentation[1500:2500])
df2 = melt(df2)
p = ggplot(df2, aes(Var1, value)) +
    geom_line(size = 1.5, color = colors[1]) +
    xlab("Distance to motif") +
    ylab("Tags") +
    theme_bw()

ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_1kb_CMonly.strandSpecific.pdf"), plot = p, height = 6, width = 6)
ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_1kb_CMonly.strandSpecific.pdf"), plot = p, height = 6, width = 8)


df2 = cbind(ChIPmentation[1800:2200])
df2 = melt(df2)
p = ggplot(df2, aes(Var1, value)) +
    geom_line(size = 1.5, color = colors[1]) + 
    xlab("Distance to motif") +
    ylab("Tags") +
    theme_bw()

ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_400bp_CMonly.strandSpecific.pdf"), plot = p, height = 6, width = 6)
ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_400bp_CMonly.strandSpecific.pdf"), plot = p, height = 6, width = 8)


df2 = cbind(ChIPmentation[1900:2100])
df2 = melt(df2)
p = ggplot(df2, aes(Var1, value)) +
    geom_line(size = 1.5, color = colors[1]) + 
    xlab("Distance to motif") +
    ylab("Tags") +
    theme_bw()

ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_200bp_CMonly.strandSpecific.pdf"), plot = p, height = 6, width = 6)
ggsave(filename = paste0(plotsDir, sample, "_peaks_signal_200bp_CMonly.strandSpecific.pdf"), plot = p, height = 6, width = 8)


### Footprinting and heatmap of CM
motifScore = read.table(paste0(bedDir, sample, ".annotation.motifStrand.bed"), sep = "\t", header = TRUE)
colnames(motifScore) = c("peak", "strand", "score")

cm2 = cm
cm2$peak = gsub(x = rownames(cm2), pattern = "_[-,+]$", replacement = "")

cm2 = merge(cm2, motifScore[, c(1,3)], by = "peak")
cm2 = cm2[cm2$score != -1e+10, ] # remove peaks with no motif


# Footprint with Centipede
library(CENTIPEDE)
centFit <- fitCentipede(Xlist = list(DNase=as.matrix(cm2[,1850:2250])), Y=cbind(rep(1, dim(cm2)[1]), cm2$score))

# sort by Posterior Prob.
window = c(1800:2200)
clusterOrder = as.data.frame(cm2[order(centFit$PostPr, decreasing = TRUE), window])

# export CM data orederd by clustering to cdt (JTV)
# write cdt files out to Java Tree view
clusterOrder$X = rownames(clusterOrder)
clusterOrder$NAME = rownames(clusterOrder)
clusterOrder$GWEIGHT = 1
colnames(clusterOrder)[1:400] = paste("X", 1:400, sep = "")
clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:400, sep = ""))]
clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 400)), clusterOrder))
write.table(clusterOrder, paste0(resultsDir, sample, ".heatmap.cdt"), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

# make PostProb plots
library(ggplot2)
df = data.frame(postProb = centFit$PostPr, sample = paste0(sample, " motifs"))

p = ggplot(df, aes(sample, postProb)) +
    geom_boxplot() +
    theme_bw()
ggsave(filename = paste0(plotsDir, sample, "_peaks.footprintsCentipede.400bp.boxplot.pdf"),
    plot = p, height = 4, width = 2.5)

p = ggplot(df, aes(sample, postProb)) +
    geom_violin(alpha=0.5, color="gray") + 
    theme_bw()
ggsave(filename = paste0(plotsDir, sample, "_peaks.footprintsCentipede.400bp.violinplot.pdf"),
    plot = p, height = 4, width = 2.5)


#### Hierarchical clustering
cmNorm = t(apply(cm[,-ncol(cm)], 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
cmNorm = cmNorm[,-ncol(cmNorm)]

window = c(1800:2200)
r = sample(1:nrow(cmNorm), 13000)

# use ward method and inverse correlation as distance measurement
hclustfunc <- function(x) hclust(x, method="ward")
distfunc <- function(x) as.dist((1-cor(t(x))))
if (sample == "PU1_K562_10mio_CM") {
    d <- distfunc(as.matrix(cm[, window])) # distance
} else if (sample == "CTCF_K562_10mio_CM") {
    d <- distfunc(as.matrix(cm[r, window])) # distance - random 13000 rows
}
fit <- hclustfunc(d) # clustering and dendrogram

save(fit, file = paste0(resultsDir, sample, "_CM_normNoLog_hclust.R", sep = ""))

# export CM data ordered by clustering to cdt (JTV)
# write cdt files out to Java Tree view
clusterOrder = as.data.frame(cm[fit$order, window])
clusterOrder$X = rownames(clusterOrder)
clusterOrder$NAME = rownames(clusterOrder)
clusterOrder$GWEIGHT = 1
colnames(clusterOrder)[1:400] = paste("X", 1:400, sep = "")
clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:400, sep = ""))]
clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 400)), clusterOrder))
write.table(clusterOrder, paste0(resultsDir, sample, "_CM_normNoLog_kmclust.cdt"), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

