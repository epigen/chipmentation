sample = '/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed'
covCM = parseBedCoverage(sample)
uniq = pd.DataFrame(covCM).T
pos = pd.DataFrame([uniq.ix[row] for row in range(1, len(uniq)) if "+" in uniq.index[row]])
neg = pd.DataFrame([uniq.ix[row][::-1] for row in range(1, len(uniq)) if "-" in uniq.index[row]])
uniqueRev = pos.append(neg)
uniqueRev.columns = ["X" + str(i) for i in uniqueRev.columns]
uniqueRev.to_csv("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv", sep = '\t')


sed -i '1s/^/X/' /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv
sed -i '1s/^/X/' /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv

mkdir -p /fhgfs/groups/lab_bock/shared/data/cage_tss/10000it /fhgfs/groups/lab_bock/shared/data/cage_tss/1000it /fhgfs/groups/lab_bock/shared/data/cage_tss/250it
cp /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv /fhgfs/groups/lab_bock/shared/data/cage_tss/10000it/
cp /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv /fhgfs/groups/lab_bock/shared/data/cage_tss/10000it/
cp /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv /fhgfs/groups/lab_bock/shared/data/cage_tss/1000it/
cp /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv /fhgfs/groups/lab_bock/shared/data/cage_tss/1000it/
cp /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv /fhgfs/groups/lab_bock/shared/data/cage_tss/250it/
cp /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv /fhgfs/groups/lab_bock/shared/data/cage_tss/250it/

cluster3 -f /fhgfs/groups/lab_bock/shared/data/cage_tss/10000it/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv -g 4 -k 5 -r 10000 -ng
cluster3 -f /fhgfs/groups/lab_bock/shared/data/cage_tss/10000it/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv -g 4 -k 5 -r 10000 -ng
cluster3 -f /fhgfs/groups/lab_bock/shared/data/cage_tss/1000it/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv -g 4 -k 5 -r 1000 -ng
cluster3 -f /fhgfs/groups/lab_bock/shared/data/cage_tss/1000it/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv -g 4 -k 5 -r 1000 -ng
cluster3 -f /fhgfs/groups/lab_bock/shared/data/cage_tss/250it/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv -g 4 -k 5 -r 100 -l
cluster3 -f hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv -g 4 -k 5 -r 100 -l


cluster3 -f  /fhgfs/groups/lab_bock/shared/data/cage_tss/GRCh37_hg19_refSeq.tss.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv -g 4 -k 5


DIR=/fhgfs/groups/lab_bock/shared/data/cage_tss/1000it
python projects/chipmentation/src/lib/count_peak_change_cluster.py \
$DIR/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage_K_G5.1.kgg \
$DIR/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage_K_G5.2.kgg

DIR=/fhgfs/groups/lab_bock/shared/data/cage_tss/250it
python projects/chipmentation/src/lib/count_peak_change_cluster.py \
$DIR/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage_K_G5.1.kgg \
$DIR/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage_K_G5.2.kgg

# hierarkical
cluster3 -f /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv -g 4 -m m -ng

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
heatplot(log2(as.matrix(df[ , -1]) + 1), dend = 'none', labRow = NA, scale = "none", dualScale = FALSE)


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
d = data.frame(cbind(chrm, start, end, peak, strand), stringsAsFactors = FALSE)

df = merge(d, df, by = "peak")

# Write Out
# export bed file of all
a = df[, c("chrm", "start", "end", "peak", "cluster", "strand")]
write.table(a, paste(dataDir, "/cage.H3K4me3_K562_500k_CM.unique.K5_all.bed", sep = ""),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# for each cluster
k = 5
xlims = list(c(23,33),c(25,32),c(-34,-24), c(-5,5), c(-7,8)) # locations of enrichment

for (clust in 0:(k - 1)) {
    # export bed file of cluster on full -60:60 window
    a = df[df$cluster == clust, c("chrm", "start", "end", "peak", "cluster", "strand")]
    write.table(a, paste(dataDir, "/cage.H3K4me3_K562_500k_CM.unique.K5", "_cluster_", clust, ".bed", sep = ""),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

    # export bed file of cluster focused on enrichment/depletion location
    x1 = xlims[[clust + 1]][1]
    x2 = xlims[[clust + 1]][2]

    a = df[df$cluster == clust, c("chrm", "start", "end", "peak", "cluster", "strand")]
    a$start = as.numeric(a$start) + 60 + x1
    a$end = as.numeric(a$end) - 60 + x2
    write.table(a, paste(dataDir, "/cage.H3K4me3_K562_500k_CM.unique.K5", "_cluster.focused_", clust, ".bed", sep = ""),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


    # export bed file of all focused on enrichment/depletion location
    a = df[, c("chrm", "start", "end", "peak", "cluster", "strand")]

    a$start = as.numeric(a$start) + 60 + x1
    a$end = as.numeric(a$end) - 60 + x2

    write.table(a, paste(dataDir, "/cage.H3K4me3_K562_500k_CM.unique.K5", "_all.focused_", clust, ".bed", sep = ""),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

}


### GET CAGE PEAK INTENSITY
#expr = read.table("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.bed", sep = "\t", header = FALSE)
expr = read.table("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.tsv", sep = "\t", header = TRUE)
expr = expr[, c(1, (ncol(expr) - 2):ncol(expr))]
colnames(expr) = c("peak", "rep1", "rep2", "rep3")
expr$int = log2(apply(expr[,-1], 1, mean) + 1)

# annotate peaks with expression 
M = merge(df, expr, by = 'peak')
M$cluster = factor(M$cluster)

# plot cage signal per cluster
library(ggplot2)

give.n <- function(x){
  return(c(y = median(x) + 2, label = length(x)))
}

p <- ggplot(M, aes(cluster, int)) +
  geom_boxplot() +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
  ylab("Log2 CAGE peak intensity") +
  theme_bw()
ggsave(filename = "~/projects/chipmentation/results/plots/TSS_unique_clusters_noNorm_expression.pdf", plot = p, height = 8, width = 12)




#### Nucleotide composition
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa

# get fasta
tail -n +2 /fhgfs/groups/lab_bock/shared/data/cage_tss/cage.H3K4me3_K562_500k_CM.unique.K5_all.bed \
| bedtools getfasta -name -fi $GENOMEREF -bed stdin -fo /fhgfs/groups/lab_bock/shared/data/cage_tss/cage.H3K4me3_K562_500k_CM.unique.K5_all.fa

# load sequences into R 
fasta = read.table("/fhgfs/groups/lab_bock/shared/data/cage_tss/cage.H3K4me3_K562_500k_CM.unique.K5_all.fa", sep = '\t')
d = data.frame(peak = fasta[seq(1, nrow(fasta), 2), ], fasta = fasta[seq(2, nrow(fasta) + 1, 2), ])
d$peak = gsub(x = d$peak, pattern = ">", replacement = "")

# merge with cluster annotation
D = merge(D, d, by = "peak")


# plot nucleotide compostition 
library(ape)

fastaCol = 7

# for all
pdf("/home/arendeiro/projects/chipmentation/results/plots/TSS_unique_noNorm_nucleotide_composition.pdf") #/home/arendeiro/projects/chipmentation/results/plots/H3K4me3_tss_k-means_clustering_averageProfile.pdf")
# order by cluster
d = D[order(D$cluster),]
m = do.call(rbind, strsplit(tolower(as.character(d[, fastaCol ])), split= ''))
colnames(m) = seq(-60, 59)
image.DNAbin(as.DNAbin(m), legend = FALSE)
dev.off()


## per cluster at site of enrichment
pdf("/home/arendeiro/projects/chipmentation/results/plots/TSS_unique_noNorm_clusters_nucleotide_composition.pdf") #/home/arendeiro/projects/chipmentation/results/plots/H3K4me3_tss_k-means_clustering_averageProfile.pdf")
colors = rainbow(kmax)
m = matrix(c(1:kmax), nrow = kmax / 4, ncol = 4, byrow = TRUE)
layout(mat = m, heights = c(rep(0.4, kmax / 4), 0.2))

par(oma = c(5,5,2,0) + 0.1,
    mar = c(1,0,0,0) + 0.1,
    pty = "s")

for (clus in 1:kmax){
    m = do.call(rbind, strsplit(tolower(as.character(D[D$cluster == clus , fastaCol ])), split= ''))
    colnames(m) = seq(-60, 59)
    image.DNAbin(as.DNAbin(m), legend = FALSE)
    #image.DNAbin(as.DNAbin(m), c("a", "t"), "dark blue", legend = FALSE)
}

dev.off()

