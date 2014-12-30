library(ggplot2)
library(gplots)

dataDir <- "/home/arendeiro/data/human/chipmentation"
projectDir <- "/home/arendeiro/projects/chipmentation"
samples <- c("PU1_K562_10mio_CM_CM15-1", "PU1_K562_500k_CM_CM11-4")

for (sample in samples) {
	sample_file <- paste(dataDir, "/bed/", sample, "_peak.2kb_coverage.bed", sep = "")

	S <- read.table(sample_file, sep = "\t", header = FALSE)
	names(S) <- c("chr", "start", "end", "peak", "bp", "PU1", "IgG", "AT", "GC", "A", "C", "G", "T", "cons")

	PU1 <- data.frame(); IGG <- data.frame(); AT <- data.frame(); GC <- data.frame(); A <- data.frame(); C <- data.frame(); G <- data.frame(); T <- data.frame(); CONS <- data.frame();

	# extract bp info
	peaks <- c()
	peakNum <- 0
	for (i in 1:nrow(S)) {
		peakName <- paste("peak_", gsub(".*peak_", "", as.character(S$peak[i])), sep = "")
		
		if (peakName %in% unique(peaks)) {
			# same peak
			j <- 6
			for (variable in c(PU1, IGG, AT, GC, A, C, G, T, CONS)) {
				variable[peakNum, bp] <- S[, j] # find how to assign to the actual variables
				j <- j + 1
			}
			bp <- bp + 1
		} else {
			# new peak
			peaks <- append(peaks, peakName)
			peakNum <- peakNum + 1
			bp <- 1
			j <- 6
			for (variable in c(PU1, IGG, AT, GC, A, C, G, T, CONS)) {
				variable[peakNum, bp] <- S[, j] # find how to assign to the actual variables
				j <- j + 1
			}
			bp <- bp + 1
		}
	}
	# add peak names
	for (variable in c(PU1, IGG, AT, GC, A, C, G, T, CONS)) {
		rownames(variable) <- peaks
	}

	# save objects
	save(PU1, file = paste(projectDir, "/results/", sample, "_peak.2kb_coverage.Rdata", sep = ""))
	save(IGG, file = paste(projectDir, "/results/", sample, "_peak.2kb_coverage.Rdata", sep = ""))
	# ...

	# get average for each bp over all peaks
	avePU1 <- colMeans(PU1)

	# plot enrichment of PU1 averages over bps
	p <- ggplot(aes(x = -2000:1999)) +
		geom_line(aes(y = colMeans(PU1), colour = "PU1")) + 
		geom_line(aes(y = colMeans(IGG), colour = "IgG")) + 
		geom_line(aes(y = colMeans(A), colour = "A")) + 
		geom_line(aes(y = colMeans(T), colour = "T")) + 
		geom_line(aes(y = colMeans(C), colour = "C")) + 
		geom_line(aes(y = colMeans(G), colour = "G")) + 
		geom_line(aes(y = colMeans(CONS), colour = "Conservation")) + 
		xlab("distance to motif") +
		theme_bw() #+ xlim(-500, 500)
	
	ggsave(filename = paste(projectDir, "/results/plots/", sample, "_peak_PU1coverage_average.pdf", sep = ""), plot = p, height = 5, width = 7)


	# Subset peaks by amount of enrichment
	peakAvePU1 <- rowMeans(PU1)

	# suset by threshold
	sub <- subset(peakAvePU1, peakAvePU1 > 10)

	# subset by quantiles
	q <- quantile(peakAvePU1)
	q1 <- subset(peakAvePU1, peakAvePU1 < q[2])
	q2 <- subset(peakAvePU1, peakAvePU1 > q[2] & peakAvePU1 <= q[3])
	q3 <- subset(peakAvePU1, peakAvePU1 > q[3] & peakAvePU1 <= q[4])
	q4 <- subset(peakAvePU1, peakAvePU1 > q[4])

	aveq1 <- colMeans(q1)
	aveq2 <- colMeans(q2)
	aveq3 <- colMeans(q3)
	aveq4 <- colMeans(q4)

	# plot enrichment of PU1 averages over bps per quantile
	p <- ggplot(aes(x = -2000:1999)) +
		geom_line(aes(y = colMeans(q1), colour = "q1")) + 
		geom_line(aes(y = colMeans(q2), colour = "q2")) + 
		geom_line(aes(y = colMeans(q3), colour = "q3")) + 
		geom_line(aes(y = colMeans(q4), colour = "q4"))
		xlab("distance to motif") +
		theme_bw() #+ xlim(-500, 500)
	
	ggsave(filename = paste(projectDir, "/results/plots/", sample, "_peak_PU1coverage_average_quantile.pdf", sep = ""), plot = p, height = 5, width = 7)


	### Cluster
	# cluster peaks based on own coverage
	d <- dist(t(PU1), method = "euclidean") # distance matrix
	fit <- hclust(d, method = "ward")
	# display dendogram
	pdf(paste(projectDir, "/results/", sample, "_peak.2kb_coverage.dendogram.pdf", sep = ""))
	plot(fit)
	groups <- cutree(fit, k=3)
	rect.hclust(fit, k=3, border="red")
	dev.off()

	pdf(paste(projectDir, "/results/", sample, "_peak.2kb_coverage.cluster.pdf", sep = ""))
	heatmap.2(log2(t(PU1) + 1), Colv = FALSE, key = TRUE, symkey = FALSE,
		density.info = "none", trace = "none",,
		labRow = NA, labCol = NA)
	dev.off()
}





