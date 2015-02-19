# First paste from excel into text editor, save as tsv (tabs will be there already)

# read in tsv
d = read.table("motif_p-values.tsv", sep = "\t", header = TRUE)
# add names to the rows (motifs)
rownames(d) = d[,1]
# remove motif column, to have only numeric values
d = d[,-1]

# load plotting library
# if not available install with
# install.packages("")
library(made4)
# try to heatmap original values
heatplot(d)
# not enough differences to separate motifs


# do it with the logarithm
pdf("motif_p-values.heatmap.pdf")
heatplot(log2(d))
dev.off()


# change colors, get Z-score instead of value
pdf("motif_p-values.heatmap.red.pdf")
heatplot(log2(d), lowcol="darkred", highcol="white", cols.default = FALSE, dend = "row", dualScale=FALSE, scale = "row")
dev.off()

pdf("motif_p-values.heatmap.blue.pdf")
heatplot(log2(d), lowcol="darkblue", highcol="white", cols.default = FALSE, dend = "row", dualScale=FALSE, scale = "row")
dev.off()

pdf("motif_p-values.heatmap.black.pdf")
heatplot(log2(d), lowcol="black", highcol="white", cols.default = FALSE, dend = "row", dualScale=FALSE, scale = "row")
dev.off()


1 - cor(log2(d))
