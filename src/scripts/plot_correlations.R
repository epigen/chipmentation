# request more resources with salloc
# salloc --partition=develop --time=10:00:00 --job-name=R --nodes=1 --mem=20000 --nodelist=n002 srun R --no-save

require(lattice)
library(ggplot2)
library(reshape2)
library(LSD)

projectDir <- "/home/arendeiro/projects/chipmentation/"
dataDir <- "/home/arendeiro/data/human/chipmentation/"

counts <- read.csv(paste(dataDir, "counts_1kb_windows.nodups.TFs.tsv", sep = ""), sep= '\t')

# exclude samples
counts <- counts[ , !grepl(x = names(counts), pattern = "PBMC")]
counts <- counts[ , !grepl(x = names(counts), pattern = "IgG")]
counts <- counts[ , !grepl(x = names(counts), pattern = "Input")]
counts <- counts[ , !grepl(x = names(counts), pattern = "_R1")]
counts <- counts[ , !grepl(x = names(counts), pattern = "_R2")]
counts <- counts[ , !grepl(x = names(counts), pattern = "PU1")]
counts <- counts[ , !grepl(x = names(counts), pattern = "CTCF")]
counts <- counts[ , !grepl(x = names(counts), pattern = "cJUN")]
counts <- counts[ , !grepl(x = names(counts), pattern = "H3K4me3")]
counts <- counts[ , !grepl(x = names(counts), pattern = "H3K27me3")]
# normalize to size 
rawCounts <- counts

i = 1
for (sample in seq(1, length(counts))) {
    rawCounts[i] <- counts[ , sample] / sum(counts[ , sample]) * 1000000
    rawCounts[i] <- rawCounts[i][!is.na(rawCounts[i])]
    i = i + 1
}

# correlate
rawCor <- cor(rawCounts)
write.table(rawCor, paste(projectDir, "results/correlations_10kb_raw.tsv", sep = ""))

# scatterplot 
pdf(paste(projectDir, "results/plots/correlations_1kb_windows_scatter_mergedReplicates.pdf", sep = ""), height = 11, width = 7)
m = matrix(c(1,2,3,4,5,5,6,7) ,nrow = 4,ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(1,1,0.2,1))

par(oma = c(5,5,2,0) + 0.1,
    mar = c(0,0,0,0) + 0.1,
    pty = "s")

smoothScatter(log2(rawCounts$H3K4me3_K562_500k_CM), log2(rawCounts$H3K4me3_K562_500k_ChIP),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n')
text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$H3K4me3_K562_500k_CM, rawCounts$H3K4me3_K562_500k_ChIP), 3))),
    cex = 1.6
)
mtext(side = 2, "500.000 cells", line = 2)
mtext(side = 3, "H3K4me3")

smoothScatter(log2(rawCounts$H3K27me3_K562_500k_CM), log2(rawCounts$H3K27me3_K562_500k_ChIP),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n', yaxt = 'n')
text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$H3K27me3_K562_500k_CM, rawCounts$H3K27me3_K562_500k_ChIP), 3))),
    cex = 1.6
)
mtext(side = 3, "H3K27me3")

smoothScatter(log2(rawCounts$H3K4me3_K562_10k_CM), log2(rawCounts$H3K4me3_K562_500k_ChIP),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
mtext(side = 2, "10.000 cells", line = 2)
text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$H3K4me3_K562_10k_CM, rawCounts$H3K4me3_K562_500k_ChIP), 3))),
    cex = 1.6
)

smoothScatter(log2(rawCounts$H3K27me3_K562_10k_CM), log2(rawCounts$H3K27me3_K562_500k_ChIP),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, yaxt = 'n')
text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$H3K27me3_K562_10k_CM, rawCounts$H3K27me3_K562_500k_ChIP), 3))),
    cex = 1.6
)
mtext(side = 1, "ChIP", outer= FALSE, at = -6.2, cex = 1, line = 3)

plot.new()
#mtext(side = 1, "ChIP", line = 0, at = 0)

smoothScatter(log2(rawCounts$H3K4me3_K562_500k_CM), log2(rawCounts$H3K4me3_K562_10k_CM),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
mtext(side = 2, "500.000 cells", line = 2)
text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$H3K4me3_K562_500k_CM, rawCounts$H3K4me3_K562_10k_CM), 3))),
    cex = 1.6
)
mtext(side = 1, "10.000 cells", line = 2, cex.lab = 1.5)

smoothScatter(log2(rawCounts$H3K27me3_K562_500k_CM), log2(rawCounts$H3K27me3_K562_10k_CM),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, yaxt = 'n')
text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$H3K27me3_K562_500k_CM, rawCounts$H3K27me3_K562_10k_CM), 3))),
    cex = 1.6
)
mtext(side = 1, "10.000 cells", line = 2, cex.lab = 1.5)

title(xlab = "ChIPmentation",
      ylab = "ChIPmentation",
      outer = TRUE, cex.lab = 1.5)
dev.off()



# biological replicates scatterplot 
pdf(paste(projectDir, "results/plots/correlations_1kb_windows_scatter_biologicalReplicates.pdf", sep = ""))
par(mfrow = c(2,2),
    oma = c(5,6,2,0) + 0.1,
    mar = c(0,0,0,0) + 0.1)

smoothScatter(log2(rawCounts$H3K4me3_K562_500k_CM_CM11.1_R1), log2(rawCounts$H3K4me3_K562_500k_CM_CM12.1_R2),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n')
text(-3, max(log2(rawCounts$H3K4me3_K562_500k_CM_CM12.1_R2)), bquote(R^2 == .(round(cor(rawCounts$H3K4me3_K562_500k_CM_CM11.1_R1, rawCounts$H3K4me3_K562_500k_CM_CM12.1_R2), 3))))
mtext(side = 2, "50.000 cells", line = 2)
mtext(side = 3, "H3K4me3")

smoothScatter(log2(rawCounts$H3K27me3_K562_500k_CM_CM11.2_R1), log2(rawCounts$H3K27me3_K562_500k_CM_CM12.2_R2),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n', yaxt = 'n')
text(-3, max(log2(rawCounts$H3K27me3_K562_500k_CM_CM12.2_R2)), bquote(R^2 == .(round(cor(rawCounts$H3K27me3_K562_500k_CM_CM11.2_R1, rawCounts$H3K27me3_K562_500k_CM_CM12.2_R2), 3))))
mtext(side = 3, "H3K27me3")

smoothScatter(log2(rawCounts$H3K4me3_K562_10k_CM_CM11.5_R1), log2(rawCounts$H3K4me3_K562_10k_CM_CM12.5_R2),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
mtext(side = 2, "10.000 cells", line = 2)
text(-3, max(log2(rawCounts$H3K4me3_K562_10k_CM_CM12.5_R2)), bquote(R^2 == .(round(cor(rawCounts$H3K4me3_K562_10k_CM_CM11.5_R1, rawCounts$H3K4me3_K562_10k_CM_CM12.5_R2), 3))))

smoothScatter(log2(rawCounts$H3K27me3_K562_10k_CM_CM11.6_R1), log2(rawCounts$H3K27me3_K562_10k_CM_CM12.6_R2),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, yaxt = 'n')
text(-3, max(log2(rawCounts$H3K27me3_K562_10k_CM_CM12.6_R2)), bquote(R^2 == .(round(cor(rawCounts$H3K27me3_K562_10k_CM_CM11.6_R1, rawCounts$H3K27me3_K562_10k_CM_CM12.6_R2), 3))))

title(xlab = "Replicate 1",
      ylab = "Replicate 2",
      outer = TRUE, cex.lab = 1.5)
dev.off()



# PBMCs scatterplot 
pdf("correlations_1kb_windows_scatter_PBMCs_no100pg_axis.pdf")
par(mfrow = c(3,1))

smoothScatter(log2(rawCounts$H3K4me3_500k_PBMC_ChIP_ChIP8.8_R1), log2(rawCounts$H3K4me3_500k_PBMC_CM_next_100pg_ChIP8.8_R1),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
text(-2, max(log2(rawCounts$H3K4me3_500k_PBMC_CM_next_100pg_ChIP8.8_R1)), bquote(R^2 == .(round(cor(rawCounts$H3K4me3_500k_PBMC_ChIP_ChIP8.8_R1, rawCounts$H3K4me3_500k_PBMC_CM_next_100pg_ChIP8.8_R1), 3))))
mtext(side = 2, "ChIP", line = 2)
mtext(side = 3, "100pg")
#abline(lm(log2(rawCounts$H3K4me3_500k_PBMC_ChIP_ChIP8.8_R1) ~ log2(rawCounts$H3K4me3_500k_PBMC_CM_next_100pg_ChIP8.8_R1)))

fit <- glm(log2(rawCounts$H3K4me3_500k_PBMC_ChIP_ChIP8.8_R1), log2(rawCounts$H3K4me3_500k_PBMC_CM_next_100pg_ChIP8.8_R1))
co <- coef(fit)
abline(fit, col="blue", lwd=2)

smoothScatter(log2(rawCounts$H3K4me3_500k_PBMC_ChIP_ChIP8.8_R1), log2(rawCounts$H3K4me3_500k_PBMC_CM_next_10pg_ChIP8.8_R1),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
text(-2, max(log2(rawCounts$H3K4me3_500k_PBMC_CM_next_10pg_ChIP8.8_R1)), bquote(R^2 == .(round(cor(rawCounts$H3K4me3_500k_PBMC_ChIP_ChIP8.8_R1, rawCounts$H3K4me3_500k_PBMC_CM_next_10pg_ChIP8.8_R1), 3))))
mtext(side = 2, "ChIP", line = 2)
mtext(side = 3, "10pg")

smoothScatter(log2(rawCounts$H3K4me3_500k_PBMC_ChIP_ChIP8.8_R1), log2(rawCounts$H3K4me3_500k_PBMC_CM_next_2pg_ChIP8.8_R1),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
mtext(side = 2, "ChIP", line = 2)
mtext(side = 3, "2pg")
text(-2, max(log2(rawCounts$H3K4me3_500k_PBMC_CM_next_2pg_ChIP8.8_R1)), bquote(R^2 == .(round(cor(rawCounts$H3K4me3_500k_PBMC_ChIP_ChIP8.8_R1, rawCounts$H3K4me3_500k_PBMC_CM_next_2pg_ChIP8.8_R1), 3))))

smoothScatter(log2(rawCounts$H3K4me3_500k_PBMC_CM_next_100pg_ChIP8.8_R1), log2(rawCounts$H3K4me3_500k_PBMC_CM_next_10pg_ChIP8.8_R1),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
mtext(side = 2, "100pg", line = 2)
mtext(side = 3, "10pg")
text(1, max(log2(rawCounts$H3K4me3_500k_PBMC_CM_next_10pg_ChIP8.8_R1)), bquote(R^2 == .(round(cor(rawCounts$H3K4me3_500k_PBMC_CM_next_100pg_ChIP8.8_R1, rawCounts$H3K4me3_500k_PBMC_CM_next_10pg_ChIP8.8_R1), 3))))

smoothScatter(log2(rawCounts$H3K4me3_500k_PBMC_CM_next_100pg_ChIP8.8_R1), log2(rawCounts$H3K4me3_500k_PBMC_CM_next_2pg_ChIP8.8_R1),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
mtext(side = 2, "100pg", line = 2)
mtext(side = 3, "2pg")
text(1, max(log2(rawCounts$H3K4me3_500k_PBMC_CM_next_2pg_ChIP8.8_R1)), bquote(R^2 == .(round(cor(rawCounts$H3K4me3_500k_PBMC_CM_next_100pg_ChIP8.8_R1, rawCounts$H3K4me3_500k_PBMC_CM_next_2pg_ChIP8.8_R1), 3))))

smoothScatter(log2(rawCounts$H3K4me3_500k_PBMC_CM_next_10pg_ChIP8.8_R1), log2(rawCounts$H3K4me3_500k_PBMC_CM_next_2pg_ChIP8.8_R1),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
mtext(side = 2, "10pg", line = 2)
mtext(side = 3, "2pg")
text(-1, max(log2(rawCounts$H3K4me3_500k_PBMC_CM_next_2pg_ChIP8.8_R1)), bquote(R^2 == .(round(cor(rawCounts$H3K4me3_500k_PBMC_CM_next_10pg_ChIP8.8_R1, rawCounts$H3K4me3_500k_PBMC_CM_next_2pg_ChIP8.8_R1), 3))))

dev.off()






# TFs
# scatterplot 
pdf(paste(projectDir, "results/plots/correlations_1kb_windows_scatter_mergedReplicates.noDups.log.TFs.pdf", sep = ""))
m = matrix(c(1,2,3,3) ,nrow = 2,ncol = 2,byrow = TRUE)
layout(mat = m, heights = c(0.4, 0.02))

par(oma = c(5,6,2,0) + 0.1,
    mar = c(0,0,0,0) + 0.1,
    pty = "s")

smoothScatter(log2(rawCounts$PU1_K562_10mio_CM), log2(rawCounts$PU1_K562_10mio_ChIP),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
text(par('usr')[1] + 1.5, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$PU1_K562_10mio_CM, rawCounts$PU1_K562_10mio_ChIP), 3)))
)

mtext(side = 2, "10^7 cells", line = 2)
mtext(side = 3, "PU1")

smoothScatter(log2(rawCounts$CTCF_K562_10mio_CM), log2(rawCounts$CTCF_K562_10mio_ChIP),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, yaxt = 'n')
text(par('usr')[1] + 1.8, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$CTCF_K562_10mio_CM, rawCounts$CTCF_K562_10mio_ChIP), 3)))
)
mtext(side = 3, "CTCF")

plot.new()
mtext(side = 1, "10^7 cells", line = 2)


title(xlab = "ChIP",
      ylab = "ChIPmentation",
      outer = TRUE, cex.lab = 1.5)
dev.off()


# biological replicates scatterplot 
pdf(paste(projectDir, "results/plots/correlations_1kb_windows_scatter_biologicalReplicates.TFs.pdf", sep = ""))
m = matrix(c(1,2,3,4,5,5) ,nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m, heights = c(0.4, 0.4, 0.02))

par(oma = c(5,6,2,0) + 0.1,
    mar = c(0,1,0,0) + 0.1,
    pty = "s")

# ChIP
smoothScatter(log2(rawCounts$PU1_K562_10mio_ChIP_CM15.5_R1), log2(rawCounts$PU1_K562_10mio_ChIP_CM16.4_R2),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n')
text(par('usr')[1] + 1.5, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$PU1_K562_10mio_ChIP_CM15.5_R1, rawCounts$PU1_K562_10mio_ChIP_CM16.4_R2), 3)))
)
mtext(side = 2, "ChIP", line = 3)
mtext(side = 3, "PU.1")

smoothScatter(log2(rawCounts$CTCF_K562_10mio_ChIP_CM15.6_R1), log2(rawCounts$CTCF_K562_10mio_ChIP_CM16.5_R2),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, xaxt = 'n', yaxt = 'n')
text(par('usr')[1] + 1.5, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$CTCF_K562_10mio_ChIP_CM15.6_R1, rawCounts$CTCF_K562_10mio_ChIP_CM16.5_R2), 3)))
)
mtext(side = 3, "CTCF")

# CM
smoothScatter(log2(rawCounts$PU1_K562_10mio_CM_CM15.1_R1), log2(rawCounts$PU1_K562_10mio_CM_CM16.1_R2),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE)
text(par('usr')[1] + 1.5, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$PU1_K562_10mio_CM_CM15.1_R1, rawCounts$PU1_K562_10mio_CM_CM16.1_R2), 3)))
)
mtext(side = 2, "ChIPmentation", line = 3)

smoothScatter(log2(rawCounts$CTCF_K562_10mio_CM_CM15.2_R1), log2(rawCounts$CTCF_K562_10mio_CM_CM16.2_R2),
    col = rgb(104,104,104,50 , maxColorValue = 255), pch = 16, nrpoints = 0, ann = FALSE, yaxt = 'n')
text(par('usr')[1] + 1.5, par('usr')[4] - 0.5,
    bquote(R^2 == .(round(cor(rawCounts$CTCF_K562_10mio_CM_CM15.2_R1, rawCounts$CTCF_K562_10mio_CM_CM16.2_R2), 3)))
)

plot.new()
mtext(side = 1, "10^7 cells", line = 2)

title(xlab = "Replicate 1",
      ylab = "Replicate 2",
      outer = TRUE, cex.lab = 2)
dev.off()







#### Scatterplots with all data between replicates
# exclude PBMCs
repl <- counts[ , !grepl(x = names(counts), pattern = "PBMC")]
# exclude IgG
repl <- repl[ , seq(1, length(rep), 2)]

# plot all in one sheet
png(paste(projectDir, "results/plots/correlations_1kb_windows_scatter.png", sep = ""), units="in", width = 28, height = 28, res=600)

par(mfrow = c((length(repl) - 1)  / 5, (length(repl) - 1) / 5))

for (sample in seq(1, length(rep) / 2, 2)){
    heatscatter(log(repl[ , sample]), log(repl[ , sample + 1]), xlab = names(repl)[sample], ylab = names(repl)[sample + 1])
}
dev.off()

#rawCor <- cor(rawCounts[seq(1, length(counts), 2)])
rawCor <- cor(rawCounts)
write.table(rawCor, paste(projectDir, "results/correlations_10kb_raw.tsv", sep = ""))

### Normalized correlations
normCounts <- counts[ , seq(1, length(counts), 2)]

i = 1
for (sample in seq(1, length(counts), 2)) {
    normCounts[i] <- ((counts[ , sample] + 1) / (counts[ , sample + 1] + 1)) * (sum(counts[ , sample + 1]) / sum(counts[ , sample]))
    i = i + 1
}
names(normCounts) <- names(counts)[seq(1, length(counts), 2)]

normCor <- cor(normCounts)
write.table(normCor, paste(projectDir, "results/correlations_10kb_norm.tsv", sep = ""))

#### Plot
rawCor <- melt(rawCor)
rawCor <- rawCor[order(rawCor$Var1),]
p <- qplot(x=Var1, y=Var2, data = rawCor, fill=value, geom="tile")
ggsave(filename = paste(projectDir, "results/plots/correlations_10kb_windows_raw.pdf", sep = ""), plot=p, width=15, height=15, units="in")

## No inputs
rawCor <- rawCor[!grepl(x = rawCor$Var1, pattern = "IgG") , ]
rawCor <- rawCor[!grepl(x = rawCor$Var2, pattern = "IgG") , ]
rawCor <- rawCor[!grepl(x = rawCor$Var1, pattern = "Input") , ]
rawCor <- rawCor[!grepl(x = rawCor$Var2, pattern = "Input") , ]
p <- qplot(x=Var1, y=Var2, data = rawCor, fill=value, geom="tile")
ggsave(filename = paste(projectDir, "results/plots/correlations_10kb_windows_raw_samples.pdf", sep = ""), plot=p, width=15, height=15, units="in")


## Normalized
normCor <- melt(normCor)
normCor <- normCor[order(normCor$Var1),]
p <- qplot(x=Var1, y=Var2, data = normCor, fill=value, geom="tile")
ggsave(filename = paste(projectDir, "results/plots/correlations_10kb_windows_norm.pdf", sep = ""), plot=p, width=15, height=15, units="in")
