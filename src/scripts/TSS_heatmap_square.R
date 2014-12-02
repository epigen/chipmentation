
R
library(ggplot2)
library(gplots)
require(made4)

plotsDir = "/home/arendeiro/projects/chipmentation/results/plots/"
resultsDir = "/home/arendeiro/projects/chipmentation/results/"
cageDir = "/fhgfs/groups/lab_bock/shared/data/cage_tss/"
cagePeaks = c(
	"hg19.cage_peak_coord_robust",
	"hg19.cage_peak_coord_robust.TATA_Annotated.TATA",
	"hg19.cage_peak_coord_robust.TATA_Annotated.TATA-less",
	"hg19.cage_peak_coord_robust.TATA_Annotated.CpG",
	"hg19.cage_peak_coord_robust.TATA_Annotated.CpG-less",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA-less",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.CpG",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.CpG-less"
	)

signals = c(
	"CM",
	"IgG",
	"ChIP",
	"DNase"
	)

cage = 6
selclus = 5


for (signal in signals){
	print(paste("Starting with", signal))
	# read in
	df = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.", signal, "coverage.csv", sep = ""), header = TRUE)
	colnames(df) <- c("peak", seq(-60, 60))
	rownames(df) = df[,1]
	df = df[,-1]

	# normalize
	dfNorm = t(apply(df, 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
	dfNorm = dfNorm[complete.cases(dfNorm),]

	# cluster
	kmfit = kmeans(dfNorm, centers = selclus, iter.max = 10000, nstart = 25)
	save(kmfit, file = paste(resultsDir, cagePeaks[cage], "_", signal, "_kmfit.R", sep = ""))

	dfClust = as.data.frame(cbind(dfNorm, rownames(dfNorm), kmfit$cluster))
	colnames(dfClust) = c(colnames(dfNorm), "peak", "cluster")
	clusterOrder = dfClust[order(dfClust[, c("cluster")]),]
	clusterOrder = as.data.frame(clusterOrder[, -c(ncol(clusterOrder) - 1, ncol(clusterOrder))])

	# pdf(paste(plotsDir, cagePeaks[cage], "_Norm", selclus, ".", sig, "coverage.pdf", sep = ""))
	# heatplot(log2(as.matrix(clusterOrder[,-ncol(clusterOrder)]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder[, ncol(clusterOrder)], dualScale= FALSE, scale = "col")
	# dev.off()

	# write cdt files out to Java Tree view
	clusterOrder$X = rownames(clusterOrder)
	clusterOrder$NAME = rownames(clusterOrder)
	clusterOrder$GWEIGHT = 1
	colnames(clusterOrder)[1:121] = paste("X", as.character(1:121), sep = "")
	clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:121, sep = ""))]
	clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 121)), clusterOrder))
	write.table(clusterOrder, paste(cageDir, cagePeaks[cage], ".120bpSlop.", "clusteredBy_", signal, ".", signal, "coverage.cdt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

	# order all per clustering
	for (sig in signals){		
		if (sig != signal){
			print(paste("Ordering", sig, "by", signal))
			df2 = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.", sig, "coverage.csv", sep = ""), header = TRUE)
			colnames(df2) <- c("peak", seq(-60, 60))
			rownames(df2) = df2[,1]
			df2 = df2[,-1]
			df2Norm = t(apply(df2, 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
			df2Norm = df2Norm[complete.cases(df2Norm),]

			dfClust2 = as.data.frame(cbind(df2Norm, rownames(df2Norm)))
			colnames(dfClust2) = c(colnames(df2Norm), "peak")

			peaks = as.data.frame(rownames(df2Norm))
			colnames(peaks) = "peak"

			ids = merge(dfClust, peaks, by = "peak")
			ids = ids[, c("peak", "cluster")]

			dfClust2 = merge(dfClust2, ids, by = "peak")
			clusterOrder = dfClust2[order(dfClust2[, c("cluster")]),]
			clusterOrder = as.data.frame(clusterOrder[, -c(1, ncol(clusterOrder))])

			# pdf(paste(plotsDir, cagePeaks[cage], "_Norm", selclus, ".", sig, "coverage.pdf", sep = ""))
			# heatplot(log2(as.matrix(clusterOrder[,-ncol(clusterOrder)]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder[, ncol(clusterOrder)], dualScale= FALSE, scale = "col")
			# dev.off()

			# write cdt files out to Java Tree view
			clusterOrder$X = rownames(clusterOrder)
			clusterOrder$NAME = rownames(clusterOrder)
			clusterOrder$GWEIGHT = 1
			colnames(clusterOrder)[1:121] = paste("X", as.character(1:121), sep = "")
			clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:121, sep = ""))]
			clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 121)), clusterOrder))
			write.table(clusterOrder, paste(cageDir, cagePeaks[cage], ".120bpSlop.", "clusteredBy_", signal, ".", sig, "coverage.cdt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

		} 
	}
}




#### % TATA - CpG
library("reshape2")
library("ggplot2")

plotsDir = "/home/arendeiro/projects/chipmentation/results/plots/"
resultsDir = "/home/arendeiro/projects/chipmentation/results/"
cageDir = "/fhgfs/groups/lab_bock/shared/data/cage_tss/"
signals = c(
	"CM",
	"IgG",
	"ChIP",
	"DNase"
	)

cage = 6
selclus = 5


annot = read.table("/home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.bed", sep = "\t", header = FALSE)

colnames(annot)[4] = "peak"
colnames(annot)[11] = "tata"
colnames(annot)[12] = "cpg"

annot = annot[, c("peak", "tata", "cpg")]

l = list()

for (signal in signals){
	print(signal)

	load(paste(resultsDir, cagePeaks[cage], "_", signal, "_kmfit.R", sep = ""))

	df = data.frame(kmfit$cluster, names(kmfit$cluster))
	colnames(df) = c("cluster", "peak")

	M = merge(annot, df, by = "peak")

	l[[signal]] = list()

	for (clus in 1:5){
		d = M[M$cluster == clus , ]

		t = length(which(d$tata == "TATA"))
		tl = length(which(d$tata == "TATA-less"))
		c = length(which(d$cpg == "CpG"))
		cl = length(which(d$cpg == "CpG-less"))
		tp = round(t / (t + tl) * 100, 2)
		tlp = round(tl / (t + tl) * 100, 2)
		cp = round(c / (c + cl) * 100, 2)
		clp = round(cl / (c + cl) * 100, 2)

		DT = data.frame(clus, tp, tlp)
		DC = data.frame(clus, cp, clp)
		
		l[[signal]][[clus]] = c(tp, tlp, cp, clp)

		DT = melt(DT, "clus")

		pie <- ggplot(DT, aes(x = as.factor(1), y = value, fill = as.factor(variable))) +
			geom_bar(stat="identity") + 
			facet_grid(. ~ clus) +
			geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), label = paste(value, "%", sep = "")), size = 5) + 
			coord_polar(theta = "y") +
			scale_fill_manual(values = c("#F50052", "#00B8B8")) +
			theme_bw() + 
			theme(axis.ticks = element_blank(), 
				axis.text.y = element_blank(),
				axis.text.x = element_blank()
			)
		ggsave(filename = paste(plotsDir, cagePeaks[cage], ".120bpSlop.", "clusteredBy_", signal, ".cluster", clus, "TATAPercent.pdf", sep = ""), plot = pie, height = 4, width = 4)

		DC = melt(DC, "clus")

		pie <- ggplot(DC, aes(x = as.factor(1), y = value, fill = as.factor(variable))) +
			geom_bar(stat="identity") + 
			facet_grid(. ~ clus) +
			geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), label = paste(value, "%", sep = "")), size = 5) + 
			coord_polar(theta = "y") +
			scale_fill_manual(values = c("#002AB3", "#B38900")) +
			theme_bw() + 
			theme(axis.ticks = element_blank(), 
				axis.text.y = element_blank(),
				axis.text.x = element_blank()
			)

		ggsave(filename = paste(plotsDir, cagePeaks[cage], ".120bpSlop.", "clusteredBy_", signal, ".cluster", clus, "CpGPercent.pdf", sep = ""), plot = pie, height = 4, width = 4)

		# pie(l[[signal]][[clus]][1:2], labels = c("TATA", "TATA-less"))
		# pie(l[[signal]][[clus]][3:4], labels = c("CpG", "CpG-less"))
	}
}






#### % TATA - CpG on TATA only tsss

cage = 7
selclus = 3



load(paste("/home/arendeiro/projects/chipmentation/results/", cagePeaks[cage], "_kmfit.R"))

df = data.frame(kmfit[[3]]$cluster, names(kmfit[[3]]$cluster))
colnames(df) = c("cluster", "peak")

M = merge(annot, df, by = "peak")


l = list()
for (clus in 1:3){
	d = M[M$cluster == clus , ]

	t = length(which(d$tata == "TATA"))
	tl = length(which(d$tata == "TATA-less"))
	c = length(which(d$cpg == "CpG"))
	cl = length(which(d$cpg == "CpG-less"))
	tp = round(t / (t + tl) * 100, 2)
	tlp = round(tl / (t + tl) * 100, 2)
	cp = round(c / (c + cl) * 100, 2)
	clp = round(cl / (c + cl) * 100, 2)

	DT = data.frame(clus, tp, tlp)
	DC = data.frame(clus, cp, clp)
	
	l[[clus]] = c(tp, tlp, cp, clp)

	DT = melt(DT, "clus")

	pie <- ggplot(DT, aes(x = as.factor(1), y = value, fill = as.factor(variable))) +
		geom_bar(stat="identity") + 
		facet_grid(. ~ clus) +
		geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), label = paste(value, "%", sep = "")), size = 5) + 
		coord_polar(theta = "y") +
		scale_fill_manual(values = c("#F50052", "#00B8B8")) +
		theme_bw() + 
		theme(axis.ticks = element_blank(), 
			axis.text.y = element_blank(),
			axis.text.x = element_blank()
		)
	ggsave(filename = paste(plotsDir, cagePeaks[cage], ".120bpSlop.", "clusteredBy_", signal, ".cluster", clus, "TATAPercent.pdf", sep = ""), plot = pie, height = 4, width = 4)

	DC = melt(DC, "clus")

	pie <- ggplot(DC, aes(x = as.factor(1), y = value, fill = as.factor(variable))) +
		geom_bar(stat="identity") + 
		facet_grid(. ~ clus) +
		geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), label = paste(value, "%", sep = "")), size = 5) + 
		coord_polar(theta = "y") +
		scale_fill_manual(values = c("#002AB3", "#B38900")) +
		theme_bw() + 
		theme(axis.ticks = element_blank(), 
			axis.text.y = element_blank(),
			axis.text.x = element_blank()
		)

	ggsave(filename = paste(plotsDir, cagePeaks[cage], ".120bpSlop.", "clusteredBy_", signal, ".cluster", clus, "CpGPercent.pdf", sep = ""), plot = pie, height = 4, width = 4)

	# pie(l[[signal]][[clus]][1:2], labels = c("TATA", "TATA-less"))
	# pie(l[[signal]][[clus]][3:4], labels = c("CpG", "CpG-less"))
}






#### LOCATION ENRICHMENT
library("reshape2")
library("ggplot2")
library("GenomicRanges")
library("simpleCache")
library("LOLA")

options(SHARE.DIR="/fhgfs/groups/lab_bock/nsheffield/share/")
source(paste0(getOption("SHARE.DIR"), "initDefault.R"))

utility("funcEnrichment.R")
utility("funcGenomeLocations.R")

loadLocationEnrichmentDatabases()


plotsDir = "/home/arendeiro/projects/chipmentation/results/plots/"
resultsDir = "/home/arendeiro/projects/chipmentation/results/"
cageDir = "/fhgfs/groups/lab_bock/shared/data/cage_tss/"
cagePeaks = c(
	"hg19.cage_peak_coord_robust",
	"hg19.cage_peak_coord_robust.TATA_Annotated.TATA",
	"hg19.cage_peak_coord_robust.TATA_Annotated.TATA-less",
	"hg19.cage_peak_coord_robust.TATA_Annotated.CpG",
	"hg19.cage_peak_coord_robust.TATA_Annotated.CpG-less",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA-less",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.CpG",
	"hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.CpG-less"
	)

signals = c(
	"CM",
	"IgG",
	"ChIP",
	"DNase"
	)

cage = 6
selclus = 5

# load annotation
annot = read.table("/home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.bed", sep = "\t", header = FALSE)
annot = annot[, 1:6]
colnames(annot) = c("chr", "start", "end", "peak", "score", "strand")


# load annotation
annot = read.table("/home/arendeiro/reference/Homo_sapiens/hg19.cage_peak_coord_robust.bed", sep = "\t", header = FALSE)
annot = annot[, 1:6]
colnames(annot) = c("chr", "start", "end", "peak", "score", "strand")


for (signal in signals){
	print(paste0("Doing signal", signal))
	# load clusters
	load(paste(resultsDir, cagePeaks[cage], "_", signal, "_kmfit.R", sep = ""))

	df = data.frame(kmfit$cluster, names(kmfit$cluster))
	colnames(df) = c("cluster", "peak")

	M = merge(annot, df, by = "peak")

	s = (M$start + M$end) / 2
	M$end = ((M$start + M$end) / 2) + 1
	M$start = s

	universe = makeGRangesFromDataFrame(M[, c("chr", "start", "end", "peak", "score", "strand")], keep.extra.columns = TRUE)
	universe = promoters(universe, upstream = 1000, downstream = 100)

	for (clus in 1:5){
		print(paste0("Doing signal", signal, ", cluster", clus))
		d = M[M$cluster == clus , c("chr", "start", "end", "peak", "score", "strand")]
		testSet = makeGRangesFromDataFrame(d, keep.extra.columns = TRUE)
		testSet = promoters(testSet, upstream = 1000, downstream = 100)

		# check appropriatness
		checkUniverseAppropriateness(testSet, universe)
		# test enrichment
		locResults = locationEnrichment(testSet, universe)

		#Output results:
		writeCombinedEnrichment(locResults, outFolder = paste0(resultsDir, "locationResults_1100bp_robustUniverse/locationResults_", signal, "_cluster", clus), includeSplits = TRUE)
	}
}

