# request more resources with salloc
# salloc --partition=develop --time=10:00:00 --job-name=R --nodes=1 --mem=20000 --nodelist=n001 srun R --no-save

require(lattice)
library(ggplot2)
library(reshape2)
library(LSD)

projectDir <- "/home/arendeiro/projects/chipmentation/"
dataDir <- "/home/arendeiro/data/human/chipmentation/"

counts <- read.table(paste(dataDir, "counts_10kb_windows.tsv", sep = ""), header = TRUE)


#### Scatterplots with all data between replicates
# exclude PBMCs
repl <- counts[ , !grepl(x = names(counts), pattern = "PBMC")]
# exclude IgG
repl <- repl[ , seq(1, length(rep), 2)]

# plot all in one sheet
png(paste(projectDir, "results/plots/correlations_10kb_windows_scatter.png", sep = ""), units="in", width = 28, height = 28, res=600)

par(mfrow = c((length(repl) - 1)  / 5, (length(repl) - 1) / 5))

for (sample in seq(1, length(rep) / 2, 2)){
	heatscatter(log(repl[ , sample]), log(repl[ , sample + 1]), xlab = names(repl)[sample], ylab = names(repl)[sample + 1])
}
dev.off()


### Raw correlations
rawCounts <- counts

i = 1
for (sample in seq(1, length(counts))) {
	rawCounts[i] <- counts[ , sample] / sum(counts[ , sample])
	i = i + 1
}

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
