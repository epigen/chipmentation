### CAGE

### get data
# new robust set - all tissues
# this link exactly: http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_coord_robust.bed.gz
/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.bed

# new robust set - K562
# follow this link exactly: http://fantom.gsc.riken.jp/5/tet/#!/search/hg19.cage_peak_tpm_ann_decoded.osc.txt.gz?c=0&c=559&c=560&c=561
/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.tsv



### Filter out non expressed genes (Encode set only)
# tsss with at least one 0 in one of the three replicates
awk -v OFS='\t' '{ if ($8 != 0 || $9 != 0 || $10 != 0) print $1,$8,$9,$10 }' /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.tsv \
| python /home/arendeiro/projects/chipmentation/src/lib/split_string_to_bed_positions.py \
> /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.bed


### Intersect both cage sets with H3K4me3 500k CM
# bedtools intersect -wa \
# -a /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.bed \
# -b /fhgfs/groups/lab_bock/shared/data/peaks/H3K4me3_K562_500k_CM_peaks/H3K4me3_K562_500k_CM_peaks.narrowPeak \
# > /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.bed

# bedtools intersect -wa \
# -a /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.bed \
# -b /fhgfs/groups/lab_bock/shared/data/peaks/H3K4me3_K562_500k_CM_peaks/H3K4me3_K562_500k_CM_peaks.narrowPeak \
# > /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.bed

# get closest tss to peak (first one) - select only the one that actually intersects the peak (distance == 0)
bedtools closest -d -t first \
-a /fhgfs/groups/lab_bock/shared/data/peaks/H3K4me3_K562_500k_CM_peaks/H3K4me3_K562_500k_CM_peaks.narrowPeak \
-b /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.bed \
| cut -f 11,12,13,14,15,16,17,18,19,20 \
| awk -v OFS='\t' '{ if ($10 == 0) print $1,$2,$3,$4,$5,$6}' \
> /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.bed

# normal gene TSSs
bedtools intersect -wa \
-a reference/Homo_sapiens/GRCh37_hg19_refSeq.tss.bed \
-b /fhgfs/groups/lab_bock/shared/data/peaks/H3K4me3_K562_500k_CM_peaks/H3K4me3_K562_500k_CM_peaks.narrowPeak \
> /fhgfs/groups/lab_bock/shared/data/cage_tss/GRCh37_hg19_refSeq.tss.H3K4me3_K562_500k_CM_peaks.bed


### Count reads on -60 to 60 windows around tsss
# Get reads to bed format and 5' position of reads
bedtools bamtobed -i /home/arendeiro/data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam \
| python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_CM.bed

# Get middle position of CAGE peak
# Make windows around TSSs
# Count reads in this windows
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt

for SET in hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique GRCh37_hg19_refSeq.tss.H3K4me3_K562_500k_CM_peaks
do
    SET=GRCh37_hg19_refSeq.tss.H3K4me3_K562_500k_CM_peaks
    cat /fhgfs/groups/lab_bock/shared/data/cage_tss/${SET}.bed \
    | python /home/arendeiro/projects/chipmentation/src/lib/get_peak_center.py \
    | bedtools slop -i stdin -g $GENOMESIZE -b 60 \
    | bedtools coverage -d \
    -a /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_CM.bed \
    -b stdin \
    > /fhgfs/groups/lab_bock/shared/data/cage_tss/${SET}.120bpCoverage.bed
done

# get same columns in both files
cut -f 1,2,3,4,6,10,11 /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed \
> tmp
mv tmp /fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed


cut -f 1,2,3,6,13,14 /fhgfs/groups/lab_bock/shared/data/cage_tss/GRCh37_hg19_refSeq.tss.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed \
> tmp
mv tmp /fhgfs/groups/lab_bock/shared/data/cage_tss/GRCh37_hg19_refSeq.tss.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed

# Parse coverage file


import sys
import csv
import pandas as pd

def parseBedCoverage(bedFile):
    """
    Reads tab-delimiter file with two rows describing the chromossomes and its lengths.
    Returns dictionary of chr:sizes.
    """
    prev_peak = ''
    peaks = {}
    with open(bedFile, "rb") as f:
        reader = csv.reader(f, delimiter="\t")
        for line, row in enumerate(reader):
            chrm = str(row[0])
            start = str(row[1])
            end = str(row[2])
            peak = str(row[3])
            strand = str(row[5])
            bp = int(row[9])
            count = int(row[10])
            
            if peak != prev_peak:
                # new peak
                peaks[peak] = {}
                peaks[peak][bp] = count
                prev_peak = peak
            else:
                peaks[peak][bp] = count
    sortedPeaks = {}
    for key in sorted(peaks.iterkeys()):
        sortedPeaks[key] = peaks[key]
    return sortedPeaks

covCM = parseBedCoverage('/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed')
all = pd.DataFrame(covCM).T
all.to_csv("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv")


covCM = parseBedCoverage('/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed')
k562 = pd.DataFrame(covCM).T
k562.to_csv("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv")

sample = '/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.bed'
covCM = parseBedCoverage(sample)
unique = pd.DataFrame(covCM).T


pos = pd.DataFrame([unique.ix[row] for row in range(1, len(unique)) if "+" in unique.index[row]])
neg = pd.DataFrame([unique.ix[row][::-1] for row in range(1, len(unique)) if "-" in unique.index[row]])
uniqueRev = pos.append(neg)
uniqueRev.columns = ["X" + str(i) for i in uniqueRev.columns]
uniqueRev.to_csv("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv", sep = '\t')


refseqCM = parseBedCoverage("/fhgfs/groups/lab_bock/shared/data/cage_tss/GRCh37_hg19_refSeq.tss.H3K4me3_K562_500k_CM_peaks.120bpCoverage.bed")
refseq = pd.DataFrame(refseqCM).T
refseq.to_csv("/fhgfs/groups/lab_bock/shared/data/cage_tss/GRCh37_hg19_refSeq.tss.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv")

R
library(gplots)
require(made4)

# read data
set.seed(2)

#df = read.csv('/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv')
#df = read.csv("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv")
#df = read.csv("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.H3K4me3_K562_500k_CM_peaks.unique.120bpCoverage.tsv")
df = read.csv("/fhgfs/groups/lab_bock/shared/data/cage_tss/GRCh37_hg19_refSeq.tss.H3K4me3_K562_500k_CM_peaks.120bpCoverage.tsv")

colnames(df) <- c("peak", seq(-60, 60))
rownames(df) = df[,1]
df = df [,-1]

# align TSSs
minus = df[grepl('-', rownames(df)),]
plus = df[!grepl('-', rownames(df)),]
minusRev = as.data.frame(t(apply(minus, 1, rev)))
colnames(minusRev) = colnames(plus)

dfInv = rbind(plus, minusRev)

dfNorm = dfInv
s = sample(1:nrow(dfNorm), 2000)
heatplot(as.matrix(dfNorm[s,]), dend = 'row', labRow = NA, scale = 'none')


# feature scaling
#dfNorm = t(apply(df, 1, function(x)(x - mean(x)) / sd(x))) # scale to [-1:1]
dfNorm = t(apply(dfInv, 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
dfNorm = dfNorm[complete.cases(dfNorm),]


# kmeans
kmax = 16 # the maximum number of clusters
totwss = rep(0, kmax) # will be filled with total sum of within group sum squares
kmfit = list() # create and empty list

for (k in 1:kmax){
    print(k)
    kclus = kmeans(dfNorm, centers = k, iter.max = 1000)
    totwss[k] = kclus$tot.withinss
    kmfit[[k]] = kclus

    dfClust <- cbind(dfNorm, kmfit[[k]]$cluster)
    colnames(dfClust) = c(colnames(dfNorm), "cluster")

    # order rows by cluster
    clusterOrder = dfClust[order(dfClust[,ncol(dfClust)]),]

    # select 2000 random rows still keeping order
    s = sample(1:nrow(clusterOrder), 2000)
    rnd = clusterOrder[s,]
    rnd = rnd[order(rnd[,ncol(rnd)]),]

    # plot heatmap
    # pdf(paste("/home/arendeiro/projects/chipmentation/results/plots/TSS_cluster_", k, ".pdf", sep = ""))
    # heatplot(as.matrix(clusterOrder[,-ncol(clusterOrder)]), dend = 'none', labRow = NA, classvec = clusterOrder[, ncol(clusterOrder)])
    # dev.off()
    pdf(paste("/home/arendeiro/projects/chipmentation/results/plots/TSS_unique_cluster_noNorm_subset_", k, ".pdf", sep = ""))
    heatplot(as.matrix(rnd[,-ncol(rnd)]), dend = 'none', labRow = NA, classvec = rnd[, ncol(rnd)])
    dev.off()
}

pdf("/home/arendeiro/projects/chipmentation/results/plots/TSS_unique_clusters_noNorm.pdf") #/home/arendeiro/projects/chipmentation/results/plots/H3K4me3_tss_k-means_clustering_averageProfile.pdf")
colors = rainbow(kmax)
m = matrix(c(1:kmax, rep((kmax + 1), 4)), nrow = (kmax + 4) / 4, ncol = 4, byrow = TRUE)
layout(mat = m, heights = c(rep(0.4, kmax / 4), 0.2))

for (k in 1:kmax){
    nclus = k
    clustered = cbind(dfNorm, kmfit[[k]]$cluster)
    colnames(clustered) <- c(seq(-60, 60), "cluster")
    for (clus in 1:nclus){
        if (clus == 1){
            par(mar = c(2,2,1,1))
            plot(x = colnames(clustered)[-ncol(clustered)],
                y = colMeans(clustered[clustered[,ncol(clustered)] == clus , -ncol(clustered)]),
                type = 'l',
                main = paste(k, 'clusters'),
                xlab = NA,
                ylab = NA,
                col = colors[clus],
                ylim = c(-1, 1)
                )
        } else {
            lines(x = as.numeric(colnames(clustered)[-ncol(clustered)]),
                y = colMeans(clustered[clustered[,ncol(clustered)] == clus , -ncol(clustered)]),
                col = colors[clus],
                type = 'l'
                )
        }
    }
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top", inset = 0, legend = 1:kmax, col = colors, bg = "white", lwd = 2, ncol = 8,horiz = FALSE)

dev.off()

# dist.pear <- function(x) as.dist(1-cor(t(x)))
# hclust.ave <- function(x) hclust(x, method="average")

# heatmap.2(as.matrix(rnd[,-ncol(rnd)]),
#     RowSideColors = bluered(length(rnd[, ncol(rnd)])), 
#     Colv = NA,
#     Rowv = NA,
#     trace = "none",
#     labRow = NA, dendrogram=c("none"), scale = c('none'),
#     symbreak=TRUE, col=bluered(256),
#     distfun=dist.pear, hclustfun=hclust.ave
# )

# calculate the adjusted R-squared
n = nrow(dfNorm)
rsq = 1-(totwss*(n-1))/(totwss[1]*(n-seq(1,kmax)))

# calculate aic
kmeansAIC = function(fit){
    m = ncol(fit$centers)
    n = length(fit$cluster)
    k = nrow(fit$centers)
    D = fit$tot.withinss
    return(D + 2*m*k)
}

aic = sapply(kmfit, kmeansAIC)

pdf("/home/arendeiro/projects/chipmentation/results/plots/TSS_unique_clusters_noNorm_R2-AIC.pdf")
require("sfsmisc")
mult.fig(2)
plot(seq(1,kmax), rsq, xlab = "Number of clusters", ylab = "Adjusted R-squared", pch=20, cex=2)
plot(seq(1,kmax), aic, xlab = "Number of clusters", ylab = "AIC", pch=20, cex=2)
dev.off()




### GET CAGE PEAK INTENSITY
expr = read.table("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.bed", sep = "\t", header = FALSE)
colnames(expr) = c("chrm", "start", "end", "peak", "strand", "rep1", "rep2", "rep3")

int = cbind(expr, log2(rowMeans(expr[,6:8])))
colnames(int) = c(colnames(expr), "int")

# annotate peaks with clusters 
selclus = 7
dfClust <- cbind(dfNorm, kmfit[[selclus]]$cluster)
dfClust <- as.data.frame(cbind(rownames(dfClust), dfClust))
colnames(dfClust) = c("peak", colnames(dfNorm), "cluster")

M = merge(int, dfClust, by = 'peak')
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


### Export with cluster annotation for all numbers of clusters
nclus = selclus
D = M[, c("chrm", "start", "end", "peak", "strand", "cluster")]

# reconstitute 120bp window
center = round( (D$start + D$end) / 2)
D$start = center - 60
D$end = center + 60

write.table(D, paste("/fhgfs/groups/lab_bock/shared/data/cage_tss/tss_unique_peaks_noNorm_cage_", nclus, "clusters.tsv",sep = ''), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

for (clus in 1:nclus){
    write.table(D[D$cluster == clus , ], paste("/fhgfs/groups/lab_bock/shared/data/cage_tss/tss_unique_peaks_noNorm_cage_", nclus, "clusters_cluster", clus, ".tsv", sep = ''),
    sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE) # 1
}



#exit R
OVERLAPS=/home/arendeiro/projects/chipmentation/results/tss_clusters_noNorm_tf_overlap.txt
rm $OVERLAPS
for tf in wgEncodeAwgTfbsSydhK562TbpIggmusUniPk.narrowPeak wgEncodeSydhTfbsK562TbpIggmusPk.narrowPeak wgEncodeAwgTfbsSydhK562Yy1UcdUniPk.narrowPeak wgEncodeSydhTfbsK562Yy1UcdPk.narrowPeak wgEncodeAwgTfbsHaibK562Taf1V0416101UniPk.narrowPeak wgEncodeAwgTfbsHaibK562Taf7sc101167V0416101UniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gtf2f1ab28179IggrabUniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gtf2bUniPk.narrowPeak K562_TFIIF_narrowPeak.bed K562_TFIIB_narrowPeak.bed
do
    for nclus in {1..7}
    do
        echo $tf $nclus
        y=`bedtools intersect -u \
        -a /fhgfs/groups/lab_bock/shared/data/cage_tss/tss_unique_peaks_noNorm_cage_7clusters_cluster${nclus}.tsv \
        -b /home/arendeiro/data/human/encode/chip-seq/${tf} \
        | wc -l`
        n=`bedtools intersect -v \
        -a /fhgfs/groups/lab_bock/shared/data/cage_tss/tss_unique_peaks_noNorm_cage_7clusters_cluster${nclus}.tsv \
        -b /home/arendeiro/data/human/encode/chip-seq/${tf} \
        | wc -l`
        r=`echo "$y $n" | awk '{printf "%.2f \n", ($1+1)/($2+1)}'`
        echo $tf $nclus $y $n $r >> $OVERLAPS
    done
done

# Plot overlaps
R
library(ggplot2)
overlaps = read.table("/home/arendeiro/projects/chipmentation/results/tss_clusters_noNorm_tf_overlap.txt", sep = ' ', header = FALSE)
colnames(overlaps) = c("tf", "cluster", "yes", "no", "ratio")

p <- ggplot(overlaps, aes(cluster, ratio)) +
geom_point() +
facet_grid(tf ~ . , scales = "free") +
theme_bw()
ggsave(filename = "~/projects/chipmentation/results/plots/TSS_unique_clusters_noNorm_overlap_with_TFs.pdf", plot = p)



#### Nucleotide composition
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa

# get fasta
tail -n +2 /fhgfs/groups/lab_bock/shared/data/cage_tss/tss_unique_peaks_noNorm_cage_7clusters.tsv \
| bedtools getfasta -name -fi $GENOMEREF -bed stdin -fo /fhgfs/groups/lab_bock/shared/data/cage_tss/tss_unique_peaks_noNorm_cage_7clusters.fa

# load sequences into R 
fasta = read.table("/fhgfs/groups/lab_bock/shared/data/cage_tss/tss_unique_peaks_noNorm_cage_7clusters.fa", sep = '\t')
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








### OLD STUUFFF
plot(colnames(df), colMeans(df), type = 'l')

plot(colnames(dfNorm), colMeans(dfNorm), type = 'l')


write.table(dfNorm, "/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.120bpCoverage.oriented.tsv",
    col.names = TRUE, row.names = TRUE, quote = FALSE, sep = '\t')



cluster3 -f hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.120bpCoverage.oriented.tsv -k 8 -g 7 -m a
cluster3 -f hg19.cage_peak_coord_robust.Encode.K562.expressed.H3K4me3_K562_500k_CM_peaks.120bpCoverage.oriented.tsv -g 7 -s -x 8 -y 8
