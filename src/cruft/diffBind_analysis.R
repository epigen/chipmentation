library(DiffBind)

projectDir <- "/home/arendeiro/projects/chipmentation/"
dataDir <- "/home/arendeiro/data/human/chipmentation/"

samples = read.csv(paste(projectDir, "samples_peaks.txt", sep = ""), header = FALSE, sep = "\t")

sampleInfo <- data.frame()

for (sample in samples[,1]) {
	sample = as.character(sample)
	S = strsplit(sample, "_")[[1]]
	control = ifelse(grepl("Encode", x = sample),
		paste("Input", S[2], S[3], S[4], S[5], S[6], sep = "_"),
		paste("IgG", S[2], S[3], S[4], S[5], S[6], sep = "_")
	)
	
	df <- data.frame(SampleID = sample,
						#Tissue = S[2], # ignore cell line and add N. of cells
						Tissue = ifelse(grepl("500", x = S[3]),
							"k500",
							"k10"
						),
						Factor = S[1],
						Condition = S[4],
						Treatment = "None",
						Replicate = S[6],
						bamReads = paste(dataDir, "bed/", sample, ".bam", sep = ""),
						ControlID = control,
						bamControl = paste(dataDir, "bed/", control, ".bam", sep = ""),
						Peaks = ifelse(S[1] != "H3K27me3",
							paste(dataDir, "peaks/", sample, ".narrowPeak", sep = ""),
							paste(dataDir, "peaks/", sample, ".broadPeak", sep = "")
						),						
						PeakCaller = ifelse(S[1] != "H3K27me3",
							"narrow",
							"broad"
						),
						Technique = S[4],
						noCells = S[3],
						Experiment = S[5]
	)
	sampleInfo <- rbind(sampleInfo, df)
}

write.csv(sampleInfo, paste(projectDir, "diffBind_sample_info.csv", sep = ""))

#sampleInfo <- sampleInfo[1:4,]
### Goal: derive consensus peaksets based on occupancy

# Get list of samples
# and determine similarities of peak occupancy based on intervals
pdf(paste(projectDir, "results/plots/", "diffBind_occupancy.pdf", sep = ""))
peaks = dba(sampleSheet=sampleInfo)
dev.off()

pdf(paste(projectDir, "results/plots/", "diffBind_occupancy2.pdf", sep = ""))
plot(peaks)
dev.off()

# Check for global overlap between peak sets (THIS IS NOT differential binding)
olap.rate = dba.overlap(peaks, mode = DBA_OLAP_RATE)

pdf(paste(projectDir, "results/plots/", "diffBind_overlapRates.pdf", sep = ""))
plot(olap.rate, type= 'b', ylab = '# peaks', xlab = 'Overlap at least this many peaksets')
dev.off()

# Can also inspect within a subset (e.g. TF and/or TISSUE)
#dba.overlap(peaks, peaks$masks$Technique & peaks$masks$Factor, mode = DBA_OLAP_RATE)

# Derive consensus peaksets for each tissue
# and with overlap threshold adjusted for each tissue (as opposed to globally)
allPeaks = dba.peakset(
	peaks,
	consensus = c(DBA_CONDITION, DBA_TISSUE, DBA_FACTOR),
	minOverlap = 2 # 0.66 requires that 2thirds of the peaks are common within a overlap group
)

# get only consensus peaks
consensusPeaks = dba(allPeaks, mask = allPeaks$masks$Consensus, bCorPlot = FALSE)

# plot venn
pdf(paste(projectDir, "results/plots/", "diffBind_overlapTissues.pdf", sep = ""))
dba.plotVenn(allPeaks, allPeaks$masks$Consensus)
dev.off()

# Get consensus peaks between replicates
ChIP_H3K4me3_k500 = dba(consensusPeaks, mask = consensusPeaks$masks$ChIP, bCorPlot = FALSE)
CM_H3K4me3_k500 = dba(consensusPeaks, mask = consensusPeaks$masks$CM & consensusPeaks$masks$k500, bCorPlot = FALSE)
CM_H3K4me3_k10 = dba(consensusPeaks, mask = consensusPeaks$masks$CM & consensusPeaks$masks$k10, bCorPlot = FALSE)

# Save consensus peaks
write.table(x = as.data.frame(ChIP_H3K4me3_k500),
			file = paste(dataDir, "peaks/", "consensus_peaks_H3K4me3_K562_ChIP_500k.bed", sep = ""),
			sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
			)
write.table(x = as.data.frame(CM_H3K4me3_k500),
			file = paste(dataDir, "peaks/", "consensus_peaks_H3K4me3_K562_CM_500k.bed", sep = ""),
			sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
			)
write.table(x = as.data.frame(CM_H3K4me3_k10),
			file = paste(dataDir, "peaks/", "consensus_peaks_H3K4me3_K562_CM_10k.bed", sep = ""),
			sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
			)

# PCA 
dba.plotPCA(allPeaks, DBA_TISSUE, label = DBA_CONDITION)