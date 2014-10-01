require(lattice)

projectDir <- "/home/arendeiro/projects/chipmentation/"
dataDir <- "/home/arendeiro/data/human/chipmentation/"

counts <- read.table(paste(dataDir, "counts_1kb_windows.tsv", sep = ""), header = TRUE)

### Raw correlations
rawCounts <- counts[ , seq(1, length(counts), 2)]

i = 1
for (sample in seq(1, length(counts), 2)) {
	rawCounts[i] <- counts[ , sample] / sum(counts[ , sample])
	i = i + 1
}

rawCor <- cor(counts[seq(1, length(counts), 2)])
write.table(rawCor, paste(projectDir, "results/correlations_raw.tsv", sep = ""))

### Normalized correlations
normCounts <- counts[ , seq(1, length(counts), 2)]

i = 1
for (sample in seq(1, length(counts), 2)) {
	normCounts[i] <- ((counts[ , sample] + 1) / (counts[ , sample + 1] + 1)) * (sum(counts[ , sample + 1]) / sum(counts[ , sample]))
	i = i + 1
}
names(normCounts) <- names(counts)[seq(1, length(counts), 2)]

normCor <- cor(normCounts)
write.table(normCor, paste(projectDir, "results/correlations_norm.tsv", sep = ""))

# Plot
pdf(paste(projectDir, "results/plots/correlations_1kb_windows_raw.pdf", sep = ""))
levelplot(as.matrix(rawCor))
dev.off()

pdf(paste(projectDir, "results/plots/correlations_1kb_windows_norm.pdf", sep = ""))
levelplot(as.matrix(normCor))
dev.off()