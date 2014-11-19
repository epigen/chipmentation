

sample_file = '/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed'
#sample_file = '/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed'

df = read.table(sample_file, sep = "\t", header = FALSE)
names(df) = c("chr", "start", "end", "peak", "strand", "rep1", "rep2", "rep3", "bp", "count")

m = matrix(NA, nrow = nrow(df), ncol = 121)

peaks = unique(df$peak)

prevPeak = ''
for (i in 1:nrow(df)) {
	peak = df$peak[i]

	peakNum = which(peaks == peak)

	if (peak != prevPeak) {
		# new Peak
		bp = 1
		m[peakNum , bp] = df$count[i]
		bp = bp + 1
	} else {
		m[peakNum , bp] = df$count[i]
		bp = bp + 1
	}
}
rownames(m) = peaks

# add peak names
for (variable in c(PU1, IGG, AT, GC, A, C, G, T, CONS)) {
	rownames(variable) = peaks
}