


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

cm = read.csv("/home/arendeiro/data/human/chipmentation/bed/PU1_K562_10mio_CM.peakCentered.peak_coverage.csv")

window = c(1600:2400)
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






##################### CLUSTERS #########################




# load data and normalize
cm = read.csv('pu1_CM.csv')
rownames(cm) = cm[,1]
colnames(cm) = seq(-2000, 1999)
cm = cm[,-1]

cmNorm = t(apply(cm[,-ncol(cm)], 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
cmNorm = cmNorm[,-ncol(cmNorm)]

window = c(1600:2400)
#r = sample(1:nrow(cm), 1000)

#dist2 <- function(x, ...)
#  as.dist(1-cor(t(x), method="pearson"))
hclust2 <- function(x, method="ward", ...)
  hclust(x, method=method, ...)

h = heatmap.2(log2(as.matrix(cmNorm[, window]) + 1),
    #distfun = dist2,
    hclust = hclust2,
    Colv = NA,
    trace = "none",
    labRow = NA, dendrogram=c("row"), scale = c('row'),
    symbreak=TRUE, col=bluered(256)
)

save(h, file = "PU1_peaks_heatmap_400bp_CM.R")

t = as.hclust(h$rowDendrogram)
memb <- cutree(t, h = 8)
nclus = max(unique(memb))

# order by cluster and plot again
D = data.frame(peak = names(memb), cluster = memb)
colnames(D) = c("peak", "cluster")
c = data.frame(cmNorm, peak = rownames(cmNorm))
q = merge(c, D, by = "peak")
rownames(cm) = cm[,1]
cm = cm[,-1]

color = redblue(max(unique(memb)))

h = heatmap.2(log2(as.matrix(q[, window]) + 1),
    RowSideColors = color[q$cluster],
    hclust = hclust2,
    Colv = NA,
    trace = "none",
    labRow = NA, dendrogram=c("row"), scale = c('row'),
    symbreak=TRUE#, col=bluered(256)
)


### Export with cluster annotation for all numbers of clusters

peaks = read.table("/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.narrowPeak")
peaks = peaks[ , 1:4 ]
colnames(peaks) = c("chrm", "start", "end", "peak")
df = merge(q, peaks, by = "peak")
df = df[, c("chrm", "start", "end", "peak", "cluster")]

write.table(df, "/fhgfs/groups/lab_bock/shared/data/pu1_clusters/all_clusters.bed", sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

for (clus in 1:nclus){
    write.table(df[df$cluster == clus , ], paste("/fhgfs/groups/lab_bock/shared/data/pu1_clusters/cluster_", clus, ".bed", sep = ''),
    sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
}



# kmeans
cm = read.csv('pu1_CM.csv')
rownames(cm) = cm[,1]
colnames(cm) = seq(-2000, 1999)
cm = cm[,-1]

km = kmeans(cm, centers = 6)
df = data.frame(peak = names(km$cluster), cluster = km$cluster)
cm$peak = rownames(cm)
d = merge(cm, df, by = "peak")

d = d[order(d$cluster) , ]

require(made4)
heatplot(d[,window], dend = 'none', labRow = NA, classvec = d$cluster)




#### INDEPENDENT
dir = "/fhgfs/groups/lab_bock/shared/data/pu1_clusters/"


cm = read.csv('pu1_CM.csv')
rownames(cm) = cm[,1]
colnames(cm) = seq(-2000, 1999)
cm = cm[,-1]

cmNorm = t(apply(cm[,-ncol(cm)], 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
cmNorm = cmNorm[,-ncol(cmNorm)]

window = c(1600:2400)

d <- as.dist(1-cor(t(log2(as.matrix(cmNorm[, window]) + 1)))); # distance
dd <- hclust(d, method="ward"); # clustering and dendrogram
save(dd, file = paste(dir, "hclust.R", sep = ""))

ddd <- as.dendrogram(dd)
dddd = heatmap.2(log2(as.matrix(cmNorm[, window]) + 1),
    Rowv = ddd,
    #distfun = dist2,
    hclust = hclust2,
    Colv = NA,
    trace = "none",
    labRow = NA, dendrogram=c("row"), scale = c('row'),
    symbreak=TRUE, col=bluered(256)
)


# Less groups
k = 7

c = cutree(dd, k = k)
cc = data.frame(peak = names(c), cluster = c)

pdf(paste(dir, k, "_clusters_dendrogram.pdf", sep = ""))
plot(dd, labels = FALSE, xlab = "")
rect.hclust(dd, h = height)
dev.off()


### Export with cluster annotation for all numbers of clusters
peaks = read.table("/home/arendeiro/data/human/chipmentation/bed/PU1_K562_10mio_CM.motif.bed")
colnames(peaks) = c("chrm", "start", "end", "peak")
peaks$start = peaks$start + 1999
peaks$end = peaks$end - 2000


df = merge(cc, peaks, by = "peak")
df = df[, c("chrm", "start", "end", "peak", "cluster")]

dir = "/fhgfs/groups/lab_bock/shared/data/pu1_clusters/"

write.table(df, paste(dir, as.character(k), "_clusters_all.bed", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

for (clus in 1:k){
    write.table(df[df$cluster == clus , ], paste(dir, as.character(k), "_clusters_", as.character(clus), ".bed", sep = ''),
    sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# More groups
k = 12

c = cutree(dd, height)
cc = data.frame(peak = names(c), cluster = c)

pdf(paste(dir, k, "_clusters_dendrogram.pdf", sep = ""))
plot(dd, labels = FALSE, xlab = "")
rect.hclust(dd, k = k)
dev.off()

df = merge(cc, peaks, by = "peak")
df = df[, c("chrm", "start", "end", "peak", "cluster")]

write.table(df, paste(dir, as.character(k), "_clusters_all.bed", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

for (clus in 1:k){
    write.table(df[df$cluster == clus , ], paste(dir, as.character(k), "_clusters_", as.character(clus), ".bed", sep = ''),
    sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
}








