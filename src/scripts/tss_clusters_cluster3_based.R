sample = '/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.bed'
covCM = parseBedCoverage(sample)
unique = pd.DataFrame(covCM).T
pos = pd.DataFrame([unique.ix[row] for row in range(1, len(unique)) if "+" in unique.index[row]])
neg = pd.DataFrame([unique.ix[row][::-1] for row in range(1, len(unique)) if "-" in unique.index[row]])
uniqueRev = pos.append(neg)
uniqueRev.columns = ["X" + str(i) for i in uniqueRev.columns]
uniqueRev.to_csv("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv", sep = '\t')

sed -i '1s/^/X/' /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv
cluster3 -f /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv \
-g 4 -k 5 -r 100 -ng

/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage_K_G5.cdt



R
# load read data
dataDir = '/fhgfs/groups/lab_bock/shared/data/cage_tss/'
df = read.table(paste(dataDir, "hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage_K_G5.cdt", sep = ""),
    header = TRUE,
    sep = "\t")
df = df[-1,-1]
colnames(df) <- c("peak", "weight", seq(-60, 60))
rownames(df) = df[,1]
df = df[,-2]

# heatmap just to check sanity
require(made4)
heatplot(log2(as.matrix(df[ , -1])), dend = 'none', labRow = NA, scale = "none", dualScale = FALSE)


# clustering data
clusters = read.table(paste(dataDir, "hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage_K_G5.kgg", sep = ""),
    header = TRUE,
    sep = "\t")
colnames(clusters) = c("peak", "cluster")

df = merge(df, clusters, by = "peak")

# get genomic locations
tmp = do.call(rbind, strsplit(as.character(df$peak), split= ','))[,1]
strand = do.call(rbind, strsplit(as.character(df$peak), split= ','))[,2]
chrm = do.call(rbind, strsplit(as.character(tmp), split= ':'))[,1]
tmp = do.call(rbind, strsplit(as.character(tmp), split= ':'))[,2]
s = as.numeric(do.call(rbind, strsplit(as.character(tmp), split= '\\..'))[,1])
e = as.numeric(do.call(rbind, strsplit(as.character(tmp), split= '\\..'))[,2])
peak = as.character(df$peak)

start = as.numeric(s + round((e - s) / 2) - 60)
end = as.numeric(60 + s + round((e - s) / 2))
d = data.frame(cbind(chrm, start, end, peak, strand))

df = merge(d, df, by = "peak")

# Write Out
k = 5
for (clust in 0:(k - 1)) {
    a = df[df$cluster == clust, c("chrm", "start", "end", "peak", "cluster", "strand")]

    write.table(a, paste(dataDir, "/cage.H3K4me3_K562_500k_CM.unique.K5", "_cluster_", clust, ".bed", sep = ""),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

