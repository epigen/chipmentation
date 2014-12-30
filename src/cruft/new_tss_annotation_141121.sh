

# strip header from one file
sed '/^#/ d' /home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.osc \
| tail -n +2 > /home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.bed

#
R
cageDir = "/fhgfs/groups/lab_bock/shared/data/cage_tss/"
### Intersect robust set with TATA/CpG annotation
robust = read.table(paste(cageDir, "hg19.cage_peak_coord_robust.bed", sep = ""), sep = "\t", header = FALSE)
colnames(robust)[1:4] = c("chrm", "start", "end", "peak")

# get center of cage peak
robust$start = round((robust$end + robust$start) / 2)
robust$end = robust$start + 1

annotation = read.table("/home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.bed", sep = "\t", header = FALSE)
colnames(robust)[4] = "peak"


annotation = annotation[, c(4, 11, 12)]
colnames(annotation) = c("peak", "TATA", "CpG")
M = merge(robust, annotation, by = "peak")
M = M[, c(2,3,4,5,1,6:11)]

write.table(M, paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(M[M$TATA == "TATA" , ], paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.TATA.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(M[M$TATA == "TATA-less" , ], paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.TATA-less.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(M[M$CpG == "CpG" , ], paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.CpG.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(M[M$CpG == "CpG-less" , ], paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.CpG-less.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


##### Intersect with K562 expressed
K562 = read.table(paste(cageDir, "hg19.cage_peak_coord_robust.Encode.K562.tsv", sep = ""), sep = "\t", header = TRUE)
colnames(K562Exp)[1] = "peak"

# filter expressed
K562$mean = apply(K562[, c(8:10)], 1, mean)
K562Exp = K562[K562$mean > 0.1]

M = merge(M, K562Exp, by = "peak")
M = M[, c(c(1:11), ncol(M))]
M = M[, c(2,3,4,5,1,6:12)]

write.table(M, paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(M[M$TATA == "TATA" , ], paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(M[M$TATA == "TATA-less" , ], paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA-less.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(M[M$CpG == "CpG" , ], paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.CpG.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(M[M$CpG == "CpG-less" , ], paste(cageDir, "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.CpG-less.bed", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# exit R


##### Calculate coverage for each annotation set
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt

declare -a arr=("hg19.cage_peak_coord_robust"
    "hg19.cage_peak_coord_robust.TATA_Annotated.TATA"
    "hg19.cage_peak_coord_robust.TATA_Annotated.TATA-less"
    "hg19.cage_peak_coord_robust.TATA_Annotated.CpG"
    "hg19.cage_peak_coord_robust.TATA_Annotated.CpG-less"
    "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed"
    "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA"
    "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA-less"
    "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.CpG"
    "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.CpG-less"
    )

# Get window around CAGE peak
for CAGE in "${arr[@]}"
do
    bedtools slop -b 60 \
    -i $CAGEDIR/$CAGE.bed \
    -g $GENOMESIZE \
    > $CAGEDIR/$CAGE.120bpSlop.bed
done


# Calculate coverage
## produce bed files from the 3 techniques
bedtools bamtobed -i /home/arendeiro/data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam \
| python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_CM.bed

bedtools bamtobed -i /home/arendeiro/data/human/chipmentation/mapped/merged/IgG_K562_500k_CM.bam \
| python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> /home/arendeiro/data/human/chipmentation/bed/IgG_K562_500k_CM.bed

bedtools bamtobed -i /home/arendeiro/data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_ChIP.bam \
| python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_ChIP.bed

bedtools bamtobed -i /home/arendeiro/data/human/chipmentation/mapped/merged/DNase_UWashington_K562_mergedReplicates.bam \
| python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
> /home/arendeiro/data/human/chipmentation/bed/DNase_UWashington_K562_mergedReplicates.bed

# calculate coverage
for CAGE in "${arr[@]}"
do
    for TECH in IgG # ChIP DNase CM 
    do
        sbatch /home/arendeiro/projects/chipmentation/src/scripts/cage_tss_coverage_job.sh $CAGE $TECH
        # bedtools coverage -d -a \
        # /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_CM.bed \
        # -b $CAGEDIR/$CAGE.120bpSlop.bed \
        # > $CAGEDIR/$CAGE.120bpSlop.CMcoverage.bed
    done
done

# parse coverage to a csv (also orients TSSs!)
for CAGE in "${arr[@]}"
do
    for TECH in IgG # ChIP DNase CM 
    do
        sbatch /home/arendeiro/projects/chipmentation/src/scripts/cage_tss_coverage-pythonParse_job.sh \
        $CAGEDIR/$CAGE.120bpSlop.${TECH}coverage.bed
        # python /home/arendeiro/projects/chipmentation/src/lib/parseBedCoverage.py $CAGEDIR/$CAGE.120bpSlop.CMcoverage.bed
    done
done



##### START EXPLORATORY ANALYSIS

R
library(ggplot2)
library(gplots)
require(made4)

plotsDir = "/home/arendeiro/projects/chipmentation/results/plots/"
cageDir = "/fhgfs/groups/lab_bock/shared/data/cage_tss/"
cagePeaks = c("hg19.cage_peak_coord_robust",
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


### Plot average signal for various cage annotations
colors = rainbow(length(cagePeaks))
for (cage in 2:length(cagePeaks)){
    df = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.CMcoverage.csv", sep = ""), header = TRUE)
    colnames(df) <- c("peak", seq(-60, 60))
    rownames(df) = df[,1]
    df = df [,-1]

    if (grepl(".TATA", cagePeaks[cage])) {
            tata = "TATA"
    } else if (grepl(".TATA-less", cagePeaks[cage])) {
        tata = "TATA-less"
    } else {tata = "All"}
    
    if (grepl(".CpG", cagePeaks[cage])) {
        cpg = "CpG"
    } else if (grepl(".CpG-less", cagePeaks[cage])) {
        cpg = "CpG-less"
    } else {cpg = "All"}

    if (grepl("K562_expressed", cagePeaks[cage])) {
        k562 = "K562"
    } else {k562 = "All"}

    if (cage == 2){
        df2 = data.frame(set = cagePeaks[cage],
        tata = tata,
        cpg = cpg,
        k562 = k562,
        bp = -60:60,
        counts = colMeans(df)
        )
        pdf("various_cage_annotations.pdf")
        plot(df2$bp, df2$counts, type = 'l', col = colors[cage], ylim = c(0, 2))
    } else {
        df2 = data.frame(set = cagePeaks[cage],
            tata = tata,
            cpg = cpg,
            k562 = k562,
            bp = -60:60,
            counts = colMeans(df)
            )
        lines(df2$bp, df2$counts, col = colors[cage])
    }
}
legend(30, 2, 2:length(cagePeaks), col = colors[2:length(cagePeaks)], lty=c(1,1))
dev.off()





########### CLUSTER PEAKS
cage = 6

df = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.CMcoverage.csv", sep = ""), header = TRUE)
colnames(df) <- c("peak", seq(-60, 60))
rownames(df) = df[,1]
df = df [,-1]

# feature scaling
#dfNorm = t(apply(df, 1, function(x)(x - mean(x)) / sd(x))) # scale to [-1:1]
dfNorm = t(apply(df, 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
dfNorm = dfNorm[complete.cases(dfNorm),]

pdf(paste(cageDir, cagePeaks[cage], "_Norm_hierarchical_clustering.pdf", sep = ""))
heatplot(dfNorm, dend = 'row', labRow = NA)
dev.off()


# kmeans clustering
kmax = 8 # the maximum number of clusters
totwss = rep(0, kmax) # will be filled with total sum of within group sum squares
kmfit = list() # create and empty list

for (k in 1:kmax){
    print(k)
    kclus = kmeans(dfNorm, centers = k, iter.max = 10000)
    totwss[k] = kclus$tot.withinss
    kmfit[[k]] = kclus
}

#save(kmfit, file = paste("/home/arendeiro/projects/chipmentation/results/", cagePeaks[cage], "_kmfit.R"))

load(paste("/home/arendeiro/projects/chipmentation/results/", cagePeaks[cage], "_kmfit.R"))

# calculate the adjusted R-squared and AIC
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

# plot
require("sfsmisc")
pdf(paste(plotsDir, cagePeaks[cage], "_Norm_R2-AIC.pdf", sep = ""))
    mult.fig(2)
    plot(seq(1,kmax), rsq, xlab = "Number of clusters", ylab = "Adjusted R-squared", pch=20, cex=2)
    plot(seq(1,kmax), aic, xlab = "Number of clusters", ylab = "AIC", pch=20, cex=2)
dev.off()



# heatmap of clusters
for (k in 1:kmax){
    print(k)
    dfClust = cbind(dfNorm, kmfit[[k]]$cluster)
    colnames(dfClust) = c(colnames(dfNorm), "cluster")
    clusterOrder = dfClust[order(dfClust[,ncol(dfClust)]),]


    # write cdt files out to Java Tree view
    clusterOrder = as.data.frame(clusterOrder[, -c(ncol(clusterOrder))])
    clusterOrder$X = rownames(clusterOrder)
    clusterOrder$NAME = rownames(clusterOrder)
    clusterOrder$GWEIGHT = 1
    colnames(clusterOrder)[1:121] = paste("X", 1:121, sep = "")
    clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:121, sep = ""))]
    clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 121)), clusterOrder))
    write.table(clusterOrder, paste(cageDir, cagePeaks[cage], ".120bpSlop.CMcoverage.cdt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

    # Plot all
    pdf(paste(plotsDir, cagePeaks[cage], "_Norm", k, ".pdf", sep = ""))
    heatplot(log2(as.matrix(clusterOrder[,-ncol(clusterOrder)]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder[, ncol(clusterOrder)], dualScale= FALSE, scale = "col")
    dev.off()

    # Plot 2000 random rows 
    s = sample(1:nrow(clusterOrder), 2000)
    rnd = clusterOrder[s,]
    rnd = rnd[order(rnd[,ncol(rnd)]),]
    pdf(paste(plotsDir, cagePeaks[cage], "_Norm", k, ".pdf", sep = ""))
    heatplot(as.matrix(rnd[,-ncol(rnd)]), dend = 'none', labRow = NA, classvec = rnd[, ncol(rnd)])
    dev.off()   
}



# Plot average profile of clusters
pdf(paste(plotsDir, cagePeaks[cage], "_Norm_k-means_clustering_averageProfile.pdf", sep = "")) #/home/arendeiro/projects/chipmentation/results/plots/H3K4me3_tss_k-means_clustering_averageProfile.pdf")
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
                ylim = c(0, 0.8)
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


######## Pick appropriate k - Analyse those clusters ##########
selclus = 5

### Plot pretty average profile for each cluster
# annotate peaks with clusters 
dfClust = as.data.frame(cbind(dfNorm, kmfit[[selclus]]$cluster))
dfClust = as.data.frame(cbind(rownames(dfClust), dfClust))
colnames(dfClust) = c("peak", colnames(dfNorm), "cluster")

# get mean per cluster
df = rbind(colMeans(dfClust[dfClust$cluster == 1, -c(1, ncol(dfClust))]),
    colMeans(dfClust[dfClust$cluster == 2, -c(1, ncol(dfClust))]),
    colMeans(dfClust[dfClust$cluster == 3, -c(1, ncol(dfClust))]),
    colMeans(dfClust[dfClust$cluster == 4, -c(1, ncol(dfClust))]),
    colMeans(dfClust[dfClust$cluster == 5, -c(1, ncol(dfClust))])
)
df = as.data.frame(df)
df$cluster = 1:5
df = melt(df, id.vars="cluster")
colnames(df) = c("cluster", "bp", "count")

library(ggplot2)
library("reshape2")

p = ggplot(df, aes(x = bp, y = count, group = cluster)) +
    geom_line() +
    #stat_smooth(method = "gam", formula = y ~ s(x, k = 60)) + 
    facet_wrap( ~ cluster , ncol = 1, scales = "free") +
    xlab("Distance to peak") +
    ylab("Tags") +
    #coord_cartesian(xlim = c(-20, 20)) +
    #scale_x_continuous(breaks = seq(-60, 60, 10), labels = seq(-60, 60, 10)) + 
    #scale_color_manual(values = colors) +
    theme_bw() #+
    #scale_size(range=c(0.0, 1.2), guide=FALSE)

ggsave(filename = paste(plotsDir, cagePeaks[cage], "_Norm_k-means_clusters_averageProfile.pdf", sep = ""), plot = p, height = 8, width = 24)



# heatmap ChIP, DNase signal on selected clusters

# IgG

dfIgG = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.IgGcoverage.csv", sep = ""), header = TRUE)
colnames(dfIgG) <- c("peak", seq(-60, 60))
rownames(dfIgG) = dfIgG[,1]

# merge CM rows with these
d = data.frame(peak = names(kmfit[[selclus]]$cluster), cluster = kmfit[[selclus]]$cluster)
dfIgG = merge(dfIgG, d, by = "peak")
clusterOrder = dfIgG[order(dfIgG$cluster),]

pdf(paste(plotsDir, cagePeaks[cage], "_Norm", selclus, ".IgGcoverage.pdf", sep = ""))
heatplot(log2(as.matrix(clusterOrder[,-c(1, ncol(clusterOrder))]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder$cluster, dualScale= FALSE, scale = "col")
dev.off()

# write cdt files out to Java Tree view
clusterOrder = as.data.frame(clusterOrder[, -c(1, ncol(clusterOrder))])
clusterOrder$X = rownames(clusterOrder)
clusterOrder$NAME = rownames(clusterOrder)
clusterOrder$GWEIGHT = 1
colnames(clusterOrder)[1:121] = paste("X", 1:121, sep = "")
clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:121, sep = ""))]
clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 121)), clusterOrder))
write.table(clusterOrder, paste(cageDir, cagePeaks[cage], ".120bpSlop.IgGcoverage.cdt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

# ChIP


dfChIP = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.ChIPcoverage.csv", sep = ""), header = TRUE)
colnames(dfChIP) <- c("peak", seq(-60, 60))
rownames(dfChIP) = dfChIP[,1]

# merge ChIP rows with these
d = data.frame(peak = names(kmfit[[selclus]]$cluster), cluster = kmfit[[selclus]]$cluster)
dfChIP = merge(dfChIP, d, by = "peak")
clusterOrder = dfChIP[order(dfChIP$cluster),]

pdf(paste(plotsDir, cagePeaks[cage], "_Norm", selclus, ".ChIPcoverage.pdf", sep = ""))
heatplot(log2(as.matrix(clusterOrder[,-c(1, ncol(clusterOrder))]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder$cluster, dualScale= FALSE, scale = "col")
dev.off()


# write cdt files out to Java Tree view
clusterOrder = as.data.frame(clusterOrder[, -c(1, ncol(clusterOrder))])
clusterOrder$X = rownames(clusterOrder)
clusterOrder$NAME = rownames(clusterOrder)
clusterOrder$GWEIGHT = 1
colnames(clusterOrder)[1:121] = paste("X", 1:121, sep = "")
clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:121, sep = ""))]
clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 121)), clusterOrder))
write.table(clusterOrder, paste(cageDir, cagePeaks[cage], ".120bpSlop.ChIPcoverage.cdt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

# DNase


dfDNase = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.DNasecoverage.csv", sep = ""), header = TRUE)
colnames(dfDNase) <- c("peak", seq(-60, 60))
rownames(dfDNase) = dfDNase[,1]

# merge ChIP rows with these
d = data.frame(peak = names(kmfit[[selclus]]$cluster), cluster = kmfit[[selclus]]$cluster)
dfDNase = merge(dfDNase, d, by = "peak")
clusterOrder = dfDNase[order(dfDNase$cluster),]

pdf(paste(plotsDir, cagePeaks[cage], "_Norm", selclus, ".DNasecoverage.pdf", sep = ""))
heatplot(log2(as.matrix(clusterOrder[,-ncol(clusterOrder)]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder$cluster, dualScale= FALSE, scale = "col")
dev.off()

# write cdt files out to Java Tree view
clusterOrder = as.data.frame(clusterOrder[, -c(1, ncol(clusterOrder))])
clusterOrder$X = rownames(clusterOrder)
clusterOrder$NAME = rownames(clusterOrder)
clusterOrder$GWEIGHT = 1
colnames(clusterOrder)[1:121] = paste("X", 1:121, sep = "")
clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:121, sep = ""))]
clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 121)), clusterOrder))
write.table(clusterOrder, paste(cageDir, cagePeaks[cage], ".120bpSlop.DNasecoverage.cdt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all


### Correlate with expression

# load expression
expr = read.table(paste(cageDir, cagePeaks[cage], ".bed", sep = ""), header = FALSE)
colnames(expr) = c("chrm", "start", "end", "1", "peak", "strand", "2", "3", "4", "5", "6", "expr")

M = merge(expr, dfClust, by = 'peak')
M$cluster = factor(M$cluster)

# plot cage signal per cluster
library(ggplot2)

give.n <- function(x){
  return(c(y = median(x) + 1, label = length(x)))
}

p <- ggplot(M, aes(cluster, log2(expr))) +
    geom_boxplot() +
    stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
    ylab("Log2 CAGE peak intensity") +
    theme_bw()

ggsave(filename = paste(plotsDir, cagePeaks[cage], "_Norm_expression.pdf", sep = ""),
    plot = p, height = 5, width = 5)


p <- ggplot(M, aes(cluster, log2(expr))) +
    geom_violin(alpha=0.5, color="gray")
    #geom_jitter(alpha=0.5, aes(color=cluster), position = position_jitter(width = 0.1))+
    #coord_flip()
ggsave(filename = paste(plotsDir, cagePeaks[cage], "_Norm_expression.violinPlot.pdf", sep = ""),
    plot = p, height = 4, width = 4)

#################  OVERLAP WITH BASAL TFs
### Export with cluster annotation for all numbers of clusters
D = M[, c("chrm", "start", "end", "peak", "cluster", "strand")]

#### reconstitute 120bp window
# center on CAGE peaks
center = round( (D$start + D$end) / 2)
D$start = center - 60
D$end = center + 60

# write out
write.table(D, paste(cageDir, cagePeaks[cage], "_Norm_", nclus, "_clusters.tsv",sep = ''), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

for (clus in 1:selclus){
    write.table(D[D$cluster == clus , ], paste(cageDir, cagePeaks[cage], "_Norm_", nclus, "_clusters_cluster_", clus, ".tsv", sep = ''),
    sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE) # 1
}


#exit R
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss
ANNOTATION=hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed
OVERLAPS=/home/arendeiro/projects/chipmentation/results/tss_clusters_Norm_tf_overlap.txt
rm $OVERLAPS
for tf in wgEncodeAwgTfbsSydhK562TbpIggmusUniPk.narrowPeak wgEncodeSydhTfbsK562TbpIggmusPk.narrowPeak wgEncodeAwgTfbsSydhK562Yy1UcdUniPk.narrowPeak wgEncodeSydhTfbsK562Yy1UcdPk.narrowPeak wgEncodeAwgTfbsHaibK562Taf1V0416101UniPk.narrowPeak wgEncodeAwgTfbsHaibK562Taf7sc101167V0416101UniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gtf2f1ab28179IggrabUniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gtf2bUniPk.narrowPeak K562_TFIIF_narrowPeak.bed K562_TFIIB_narrowPeak.bed
do
    for nclus in {1..5}
    do
        echo $tf $nclus
        y=`bedtools intersect -u \
        -a $CAGEDIR/${ANNOTATION}_Norm_5_clusters_cluster_${nclus}.tsv \
        -b /home/arendeiro/data/human/encode/chip-seq/${tf} \
        | wc -l`
        n=`bedtools intersect -v \
        -a $CAGEDIR/${ANNOTATION}_Norm_5_clusters_cluster_${nclus}.tsv \
        -b /home/arendeiro/data/human/encode/chip-seq/${tf} \
        | wc -l`
        r=`echo "$y $n" | awk '{printf "%.2f \n", ($1+1)/($2+1)}'`
        echo $tf $nclus $y $n $r >> $OVERLAPS
    done
done

# Plot overlaps
R
library(ggplot2)
overlaps = read.table("/home/arendeiro/projects/chipmentation/results/tss_clusters_Norm_tf_overlap.txt", sep = ' ', header = FALSE)
colnames(overlaps) = c("tf", "cluster", "yes", "no", "ratio")

p <- ggplot(overlaps, aes(cluster, ratio)) +
    geom_point() +
    facet_wrap( ~ tf, scales = "free", ncol = 2) +
    theme_bw()

ggsave(filename = paste(plotsDir, cagePeaks[cage], "_Norm_", selclus, "_clusters_overlap_with_TFs.pdf", sep = ''), plot = p)


#### NUCLEOTIDE COMPOSITION
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss
ANNOTATION=hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa

# get fasta with strand information!!!!
tail -n +2 $CAGEDIR/${ANNOTATION}_Norm_5_clusters.tsv \
| bedtools getfasta -s -name -fi $GENOMEREF -bed stdin -fo $CAGEDIR/${ANNOTATION}_Norm_5_clusters.fa

# load sequences into R 
fasta = read.table(paste(cageDir, cagePeaks[cage], "_Norm_", selclus, "_clusters.fa",sep = ''), sep = '\t')
d = data.frame(peak = fasta[seq(1, nrow(fasta), 2), ], fasta = fasta[seq(2, nrow(fasta) + 1, 2), ])
d$peak = gsub(x = d$peak, pattern = ">", replacement = "")

# merge with cluster annotation
D = merge(D, d, by = "peak")


# plot nucleotide compostition 
library(ape)

fastaCol = 7

# for all
pdf(paste(plotsDir, cagePeaks[cage], "_Norm_", selclus, "_all_nucleotide_composition.pdf", sep = ""))
# order by cluster
d = D[order(D$cluster, decreasing = TRUE),]
m = do.call(rbind, strsplit(tolower(as.character(d[, fastaCol ])), split= ''))
colnames(m) = seq(-60, 59)
image.DNAbin(as.DNAbin(m), legend = FALSE)
dev.off()


## per cluster at site of enrichment
pdf(paste(plotsDir, cagePeaks[cage], "_Norm_", selclus, "_clusters_nucleotide_composition.pdf", sep = ""))
#colors = rainbow(selclus)
m = matrix(c(1:6), nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = m, heights = c(rep(0.4, 2)))

par(oma = c(5,5,2,0) + 0.1,
    mar = c(1,0,0,0) + 0.1,
    pty = "s")

for (clus in 1:selclus){
    m = do.call(rbind, strsplit(tolower(as.character(D[D$cluster == clus , fastaCol ])), split= ''))
    colnames(m) = seq(-60, 59)
    image.DNAbin(as.DNAbin(m), legend = FALSE)
    #image.DNAbin(as.DNAbin(m), c("a", "t"), "dark blue", legend = FALSE)
}

dev.off()




### GENE ONTOLOGY ## GO of clusters
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss
ANNOTATION=hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed

bedtools closest -d -t first \
-a $CAGEDIR/${ANNOTATION}.bed \
-b /home/arendeiro/reference/Homo_sapiens/GRCh37_hg19_refSeq.tss.bed \
| awk 'BEGIN {FS='\t'} {if ($24 < 300) print }' \
| cut -f 5,16 > $CAGEDIR/${ANNOTATION}.refseqIDs.txt

# load gene info
genes = read.table(paste(cageDir, cagePeaks[cage], ".refseqIDs.txt", sep = ""), sep = '\t', header = FALSE)
colnames(genes) = c("peak", "gene")

# merge with cluster info
genes = merge(M, genes)


# get relation between refseq and ensembl ids
library(biomaRt)
library(org.Hs.eg.db)
library(GOstats)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
#listAttributes(mart)[grep('RefSeq', listAttributes(mart)$description),]

IDs = getBM(attributes = c("refseq_mrna", "refseq_ncrna", 'entrezgene'), values = genes$gene, mart = ensembl) # or (mart = mart)
save(IDs, file = "/home/arendeiro/reference/Homo_sapiens/refseq_entrez_IDs.R")
IDs = load("/home/arendeiro/reference/Homo_sapiens/refseq_entrez_IDs.R")

# annotate with ensembl ids
genesNR = data.frame(genes[grep('NR', as.character(genes$gene)),])
colnames(genesNR)[ncol(genesNR)] = 'refseq_ncrna'
genesE = merge(genesNR, IDs)
genesNM = data.frame(genes[grep('NM', as.character(genes$gene)),])
colnames(genesNM)[ncol(genesNM)] = 'refseq_mrna'
genesEE = merge(genesNM, IDs)
genesE = as.data.frame(rbind(genesE, genesEE))

# do the test
BPtestResults = list()
for (k in 1:selclus){
    print(k)
    selected = unique(genesE[kmfit[[selclus]]$cluster == k , ncol(genesE)])
    param = new("GOHyperGParams", geneIds=selected,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology="BP",pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp = hyperGTest(param)
    BPtestResults[[k]] = hyp
    #sumTable <- summary(hyp)
}

MFtestResults = list()
for (k in 1:selclus){
    print(k)
    selected = unique(genesE[kmfit[[selclus]]$cluster == k , ncol(genesE)])
    param = new("GOHyperGParams", geneIds=selected,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology="MF",pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp = hyperGTest(param)
    MFtestResults[[k]] = hyp
    #sumTable <- summary(hyp)
}

CCtestResults = list()
for (k in 1:selclus){
    print(k)
    selected = unique(genesE[kmfit[[selclus]]$cluster == k , ncol(genesE)])
    param = new("GOHyperGParams", geneIds=selected,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology="CC",pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp = hyperGTest(param)
    CCtestResults[[k]] = hyp
    #sumTable <- summary(hyp)
}


# export
for (k in 1:selclus){
    print(k)
    write.table(summary(BPtestResults[[k]]), paste(plotsDir, "../", cagePeaks[cage], "_Norm_", selclus, "_clusters_cluster", k, ".GO-BP.tsv",sep = ''),
        sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(summary(MFtestResults[[k]]), paste(plotsDir, "../", cagePeaks[cage], "_Norm_", selclus, "_clusters_cluster", k, ".GO-MF.tsv",sep = ''),
        sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(summary(CCtestResults[[k]]), paste(plotsDir, "../", cagePeaks[cage], "_Norm_", selclus, "_clusters_cluster", k, ".GO-CC.tsv",sep = ''),
        sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
}






















#### TATA Analysis

cage = 7

df = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.CMcoverage.csv", sep = ""), header = TRUE)
    colnames(df) <- c("peak", seq(-60, 60))
    rownames(df) = df[,1]
    df = df [,-1]

# feature scaling
#dfNorm = t(apply(df, 1, function(x)(x - mean(x)) / sd(x))) # scale to [-1:1]
dfNorm = t(apply(df, 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
dfNorm = dfNorm[complete.cases(dfNorm),]

# pdf(paste(cageDir, cagePeaks[cage], "_Norm_hierarchical_clustering.pdf", sep = ""))
# heatplot(dfNorm, dend = 'row', labRow = NA)
# dev.off()


# kmeans clustering
kmax = 8 # the maximum number of clusters
totwss = rep(0, kmax) # will be filled with total sum of within group sum squares
kmfit = list() # create and empty list

for (k in 1:kmax){
    print(k)
    kclus = kmeans(dfNorm, centers = k, iter.max = 10000)
    totwss[k] = kclus$tot.withinss
    kmfit[[k]] = kclus
}

save(kmfit, file = paste("/home/arendeiro/projects/chipmentation/results/", cagePeaks[cage], "_kmfit.R"))

load(paste("/home/arendeiro/projects/chipmentation/results/", cagePeaks[cage], "_kmfit.R"))

# calculate the adjusted R-squared and AIC
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

# plot
require("sfsmisc")
pdf(paste(plotsDir, cagePeaks[cage], "_Norm_R2-AIC.pdf", sep = ""))
    mult.fig(2)
    plot(seq(1,kmax), rsq, xlab = "Number of clusters", ylab = "Adjusted R-squared", pch=20, cex=2)
    plot(seq(1,kmax), aic, xlab = "Number of clusters", ylab = "AIC", pch=20, cex=2)
dev.off()



# heatmap of clusters
for (k in 1:kmax){
    print(k)
    dfClust = cbind(dfNorm, kmfit[[k]]$cluster)
    colnames(dfClust) = c(colnames(dfNorm), "cluster")
    clusterOrder = dfClust[order(dfClust[,ncol(dfClust)]),]

    # Plot all
    clusterOrder = dfClust[order(dfClust[,ncol(dfClust)]),]
    pdf(paste(plotsDir, cagePeaks[cage], "_Norm", k, ".pdf", sep = ""))
    heatplot(log2(as.matrix(clusterOrder[,-ncol(clusterOrder)]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder[, ncol(clusterOrder)], dualScale= FALSE, scale = "col")
    dev.off()

    # Plot 2000 random rows 
    # s = sample(1:nrow(clusterOrder), 2000)
    # rnd = clusterOrder[s,]
    # rnd = rnd[order(rnd[,ncol(rnd)]),]
    # pdf(paste(plotsDir, cagePeaks[cage], "_Norm", k, ".pdf", sep = ""))
    # heatplot(as.matrix(rnd[,-ncol(rnd)]), dend = 'none', labRow = NA, classvec = rnd[, ncol(rnd)])
    # dev.off()   
}



# Plot average profile of clusters
pdf(paste(plotsDir, cagePeaks[cage], "_Norm_k-means_clustering_averageProfile.pdf", sep = "")) #/home/arendeiro/projects/chipmentation/results/plots/H3K4me3_tss_k-means_clustering_averageProfile.pdf")
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
                ylim = c(0, 0.8)
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


######## Pick appropriate k - Analyse those clusters ##########
selclus = 3

### Plot pretty average profile for each cluster
# annotate peaks with clusters 
dfClust <- as.data.frame(cbind(dfNorm, kmfit[[selclus]]$cluster))
dfClust <- as.data.frame(cbind(rownames(dfClust), dfClust))
colnames(dfClust) = c("peak", colnames(dfNorm), "cluster")

# get mean per cluster
df = rbind(colMeans(dfClust[dfClust$cluster == 1, -c(1, ncol(dfClust))]),
    colMeans(dfClust[dfClust$cluster == 2, -c(1, ncol(dfClust))]),
    colMeans(dfClust[dfClust$cluster == 3, -c(1, ncol(dfClust))])
)
df = as.data.frame(df)
df$cluster = 1:selclus

library(reshape2)
df = melt(df, id.vars="cluster")
colnames(df) = c("cluster", "bp", "count")

library(ggplot2)
library("reshape2")

p = ggplot(df, aes(x = bp, y = count, group = cluster)) +
    geom_line() +
    #stat_smooth(method = "gam", formula = y ~ s(x, k = 60)) + 
    facet_wrap( ~ cluster , ncol = 1, scales = "free") +
    xlab("Distance to peak") +
    ylab("Tags") +
    #coord_cartesian(xlim = c(-20, 20)) +
    #scale_x_continuous(breaks = seq(-60, 60, 10), labels = seq(-60, 60, 10)) + 
    #scale_color_manual(values = colors) +
    theme_bw() #+
    #scale_size(range=c(0.0, 1.2), guide=FALSE)

ggsave(filename = paste(plotsDir, cagePeaks[cage], "_Norm_k-means_clusters_averageProfile.pdf", sep = ""), plot = p, height = 8, width = 24)



# heatmap ChIP, DNase signal on selected clusters

# CM
dfClust = cbind(dfNorm, kmfit[[selclus]]$cluster)
colnames(dfClust) = c(colnames(dfNorm), "cluster")
clusterOrder = dfClust[order(dfClust[,ncol(dfClust)]),]

pdf(paste(plotsDir, cagePeaks[cage], "_Norm", selclus, ".CMcoverage.pdf", sep = ""))
heatplot(log2(as.matrix(clusterOrder[,-ncol(clusterOrder)]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder[, ncol(clusterOrder)], dualScale= FALSE, scale = "col")
dev.off()

# write cdt files out to Java Tree view
clusterOrder = as.data.frame(clusterOrder[, -c(ncol(clusterOrder))])
clusterOrder$X = rownames(clusterOrder)
clusterOrder$NAME = rownames(clusterOrder)
clusterOrder$GWEIGHT = 1
colnames(clusterOrder)[1:121] = paste("X", 1:121, sep = "")
clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:121, sep = ""))]
clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 121)), clusterOrder))
write.table(clusterOrder, paste(cageDir, cagePeaks[cage], ".120bpSlop.CMcoverage.cdt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all


# CHIP

dfChIP = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.ChIPcoverage.csv", sep = ""), header = TRUE)
colnames(dfChIP) <- c("peak", seq(-60, 60))
rownames(dfChIP) = dfChIP[,1]
dfChIP = dfChIP[,-1]
dfChIP = t(apply(dfChIP, 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
dfClust = as.data.frame(cbind(dfChIP, kmfit[[selclus]]$cluster))
dfClust = dfClust[complete.cases(dfChIP),]
colnames(dfClust) = c(colnames(dfChIP), "cluster")
clusterOrder = dfClust[order(dfClust[,ncol(dfClust)]),]

pdf(paste(plotsDir, cagePeaks[cage], "_Norm", selclus, ".ChIPcoverage.pdf", sep = ""))
heatplot(log2(as.matrix(clusterOrder[,-ncol(clusterOrder)]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder[, ncol(clusterOrder)], dualScale= FALSE, scale = "col")
dev.off()

# write cdt files out to Java Tree view
clusterOrder = as.data.frame(clusterOrder[, -c(ncol(clusterOrder))])
clusterOrder$X = rownames(clusterOrder)
clusterOrder$NAME = rownames(clusterOrder)
clusterOrder$GWEIGHT = 1
colnames(clusterOrder)[1:121] = paste("X", 1:121, sep = "")
clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:121, sep = ""))]
clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 121)), clusterOrder))
write.table(clusterOrder, paste(cageDir, cagePeaks[cage], ".120bpSlop.ChIPcoverage.cdt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

# DNASE

dfDNase = read.csv(paste(cageDir, cagePeaks[cage], ".120bpSlop.DNasecoverage.csv", sep = ""), header = TRUE)
colnames(dfDNase) <- c("peak", seq(-60, 60))
rownames(dfDNase) = dfDNase[,1]
dfDNase = dfDNase[,-1]
dfDNase = t(apply(dfDNase, 1, function(x)((x - min(x)) / (max(x) - min(x))))) # scale to [0:1]
dfClust = as.data.frame(cbind(dfDNase, kmfit[[selclus]]$cluster))
dfClust = dfClust[complete.cases(dfDNase),]
colnames(dfClust) = c(colnames(dfDNase), "cluster")
clusterOrder = dfClust[order(dfClust[,ncol(dfClust)]),]

pdf(paste(plotsDir, cagePeaks[cage], "_Norm", selclus, ".DNasecoverage.pdf", sep = ""))
heatplot(log2(as.matrix(clusterOrder[,-ncol(clusterOrder)]) + 1), dend = 'none', labRow = NA, classvec = clusterOrder[, ncol(clusterOrder)], dualScale= FALSE, scale = "col")
dev.off()

# write cdt files out to Java Tree view
clusterOrder = as.data.frame(clusterOrder[, -c(ncol(clusterOrder))])
clusterOrder$X = rownames(clusterOrder)
clusterOrder$NAME = rownames(clusterOrder)
clusterOrder$GWEIGHT = 1
colnames(clusterOrder)[1:121] = paste("X", 1:121, sep = "")
clusterOrder = clusterOrder[, c("X", "NAME", "GWEIGHT", paste("X", 1:121, sep = ""))]
clusterOrder = as.data.frame(rbind(c("EWEIGHT", "", "", rep(1, 121)), clusterOrder))
write.table(clusterOrder, paste(cageDir, cagePeaks[cage], ".120bpSlop.DNasecoverage.cdt", sep = ""), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all


### Correlate with expression

# load expression
expr = read.table(paste(cageDir, cagePeaks[cage], ".bed", sep = ""), header = FALSE)
colnames(expr) = c("chrm", "start", "end", "1", "peak", "strand", "2", "3", "4", "5", "6", "expr")

dfClust = as.data.frame(cbind(dfNorm, kmfit[[selclus]]$cluster))
colnames(dfClust)[ncol(dfClust)] = "cluster"
dfClust$peak = rownames(dfClust)
M = merge(expr, dfClust, by = 'peak')
M$cluster = factor(M$cluster)

# plot cage signal per cluster
library(ggplot2)

give.n <- function(x){
  return(c(y = median(x) + 1, label = length(x)))
}

p <- ggplot(M, aes(cluster, log2(expr))) +
    geom_boxplot() +
    stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
    ylab("Log2 CAGE peak intensity") +
    theme_bw()

ggsave(filename = paste(plotsDir, cagePeaks[cage], "_Norm_expression.pdf", sep = ""),
    plot = p, height = 5, width = 5)


p <- ggplot(M, aes(cluster, log2(expr))) +
    geom_violin(alpha=0.5, color="gray")
    #geom_jitter(alpha=0.5, aes(color=cluster), position = position_jitter(width = 0.1))+
    #coord_flip()
ggsave(filename = paste(plotsDir, cagePeaks[cage], "_Norm_expression.violinPlot.pdf", sep = ""),
    plot = p, height = 4, width = 4)

#################  OVERLAP WITH BASAL TFs
### Export with cluster annotation for all numbers of clusters
D = M[, c("chrm", "start", "end", "peak", "cluster", "strand")]

#### reconstitute 120bp window
# center on CAGE peaks
center = round( (D$start + D$end) / 2)
D$start = center - 60
D$end = center + 60

# write out
write.table(D, paste(cageDir, cagePeaks[cage], "_Norm_", nclus, "_clusters.tsv",sep = ''), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

for (clus in 1:selclus){
    write.table(D[D$cluster == clus , ], paste(cageDir, cagePeaks[cage], "_Norm_", nclus, "_clusters_cluster_", clus, ".tsv", sep = ''),
    sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE) # 1
}


#exit R
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss
ANNOTATION=hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA
OVERLAPS=/home/arendeiro/projects/chipmentation/results/tss_clusters_Norm_tf_overlap_TATA.txt
rm $OVERLAPS
for tf in wgEncodeAwgTfbsSydhK562TbpIggmusUniPk.narrowPeak wgEncodeSydhTfbsK562TbpIggmusPk.narrowPeak wgEncodeAwgTfbsSydhK562Yy1UcdUniPk.narrowPeak wgEncodeSydhTfbsK562Yy1UcdPk.narrowPeak wgEncodeAwgTfbsHaibK562Taf1V0416101UniPk.narrowPeak wgEncodeAwgTfbsHaibK562Taf7sc101167V0416101UniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gtf2f1ab28179IggrabUniPk.narrowPeak wgEncodeAwgTfbsSydhK562Gtf2bUniPk.narrowPeak K562_TFIIF_narrowPeak.bed K562_TFIIB_narrowPeak.bed
do
    for nclus in {1..3}
    do
        echo $tf $nclus
        y=`bedtools intersect -u \
        -a $CAGEDIR/${ANNOTATION}_Norm_3_clusters_cluster_${nclus}.tsv \
        -b /home/arendeiro/data/human/encode/chip-seq/${tf} \
        | wc -l`
        n=`bedtools intersect -v \
        -a $CAGEDIR/${ANNOTATION}_Norm_3_clusters_cluster_${nclus}.tsv \
        -b /home/arendeiro/data/human/encode/chip-seq/${tf} \
        | wc -l`
        r=`echo "$y $n" | awk '{printf "%.2f \n", ($1+1)/($2+1)}'`
        echo $tf $nclus $y $n $r >> $OVERLAPS
    done
done

# Plot overlaps
R
library(ggplot2)
overlaps = read.table("/home/arendeiro/projects/chipmentation/results/tss_clusters_Norm_tf_overlap_TATA.txt", sep = ' ', header = FALSE)
colnames(overlaps) = c("tf", "cluster", "yes", "no", "ratio")

p <- ggplot(overlaps, aes(cluster, ratio)) +
    geom_point() +
    facet_wrap( ~ tf, scales = "free", ncol = 2) +
    theme_bw()

ggsave(filename = paste(plotsDir, cagePeaks[cage], "_Norm_", selclus, "_clusters_overlap_with_TFs.pdf", sep = ''), plot = p)


#### NUCLEOTIDE COMPOSITION
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss
ANNOTATION=hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa

# get fasta with strand information!!!!
tail -n +2 $CAGEDIR/${ANNOTATION}_Norm_5_clusters.tsv \
| bedtools getfasta -s -name -fi $GENOMEREF -bed stdin -fo $CAGEDIR/${ANNOTATION}_Norm_3_clusters.fa

# load sequences into R 
fasta = read.table(paste(cageDir, cagePeaks[cage], "_Norm_", selclus, "_clusters.fa",sep = ''), sep = '\t')
d = data.frame(peak = fasta[seq(1, nrow(fasta), 2), ], fasta = fasta[seq(2, nrow(fasta) + 1, 2), ])
d$peak = gsub(x = d$peak, pattern = ">", replacement = "")

# merge with cluster annotation
D = merge(D, d, by = "peak")


# plot nucleotide compostition 
library(ape)

fastaCol = 7

# for all
pdf(paste(plotsDir, cagePeaks[cage], "_Norm_", selclus, "_all_nucleotide_composition.pdf", sep = ""))
# order by cluster
d = D[order(D$cluster, decreasing = TRUE),]
m = do.call(rbind, strsplit(tolower(as.character(d[, fastaCol ])), split= ''))
colnames(m) = seq(-60, 59)
image.DNAbin(as.DNAbin(m), legend = FALSE)
dev.off()


## per cluster at site of enrichment
pdf(paste(plotsDir, cagePeaks[cage], "_Norm_", selclus, "_clusters_nucleotide_composition.pdf", sep = ""))
#colors = rainbow(selclus)
m = matrix(c(1:6), nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = m, heights = c(rep(0.4, 2)))

par(oma = c(5,5,2,0) + 0.1,
    mar = c(1,0,0,0) + 0.1,
    pty = "s")

for (clus in 1:selclus){
    m = do.call(rbind, strsplit(tolower(as.character(D[D$cluster == clus , fastaCol ])), split= ''))
    colnames(m) = seq(-60, 59)
    image.DNAbin(as.DNAbin(m), legend = FALSE)
    #image.DNAbin(as.DNAbin(m), c("a", "t"), "dark blue", legend = FALSE)
}

dev.off()



### GENE ONTOLOGY ## GO of clusters
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss
ANNOTATION=hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA

bedtools closest -d -t first \
-a $CAGEDIR/${ANNOTATION}.bed \
-b /home/arendeiro/reference/Homo_sapiens/GRCh37_hg19_refSeq.tss.bed \
| awk 'BEGIN {FS='\t'} {if ($24 < 300) print }' \
| cut -f 5,16 > $CAGEDIR/${ANNOTATION}.refseqIDs.txt

# load gene info
genes = read.table(paste(cageDir, cagePeaks[cage], ".refseqIDs.txt", sep = ""), sep = '\t', header = FALSE)
colnames(genes) = c("peak", "gene")

# merge with cluster info
genes = merge(M, genes)

# get relation between refseq and ensembl ids
library(biomaRt)
library(org.Hs.eg.db)
library(GOstats)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
#listAttributes(mart)[grep('RefSeq', listAttributes(mart)$description),]

# IDs = getBM(attributes = c("refseq_mrna", "refseq_ncrna", 'entrezgene'), values = genes$gene, mart = ensembl) # or (mart = mart)
# save(IDs, file = "/home/arendeiro/reference/Homo_sapiens/refseq_entrez_IDs.R")
IDs = load("/home/arendeiro/reference/Homo_sapiens/refseq_entrez_IDs.R")

# annotate with ensembl ids
genesNR = data.frame(genes[grep('NR', as.character(genes$gene)),])
colnames(genesNR)[ncol(genesNR)] = 'refseq_ncrna'
genesE = merge(genesNR, IDs)
genesNM = data.frame(genes[grep('NM', as.character(genes$gene)),])
colnames(genesNM)[ncol(genesNM)] = 'refseq_mrna'
genesEE = merge(genesNM, IDs)
genesE = as.data.frame(rbind(genesE, genesEE))

# do the test
BPtestResults = list()
for (k in 1:selclus){
    print(k)
    selected = unique(genesE[kmfit[[selclus]]$cluster == k , ncol(genesE)])
    param = new("GOHyperGParams", geneIds=selected,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology="BP",pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp = hyperGTest(param)
    BPtestResults[[k]] = hyp
    #sumTable <- summary(hyp)
}

MFtestResults = list()
for (k in 1:selclus){
    print(k)
    selected = unique(genesE[kmfit[[selclus]]$cluster == k , ncol(genesE)])
    param = new("GOHyperGParams", geneIds=selected,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology="MF",pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp = hyperGTest(param)
    MFtestResults[[k]] = hyp
    #sumTable <- summary(hyp)
}

CCtestResults = list()
for (k in 1:selclus){
    print(k)
    selected = unique(genesE[kmfit[[selclus]]$cluster == k , ncol(genesE)])
    param = new("GOHyperGParams", geneIds=selected,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology="CC",pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp = hyperGTest(param)
    CCtestResults[[k]] = hyp
    #sumTable <- summary(hyp)
}

# export
for (k in 1:selclus){
    print(k)
    write.table(summary(BPtestResults[[k]]), paste(plotsDir, "../", cagePeaks[cage], "_Norm_", selclus, "_clusters_cluster", k, ".GO-BP.tsv",sep = ''),
        sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(summary(MFtestResults[[k]]), paste(plotsDir, "../", cagePeaks[cage], "_Norm_", selclus, "_clusters_cluster", k, ".GO-MF.tsv",sep = ''),
        sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
    write.table(summary(CCtestResults[[k]]), paste(plotsDir, "../", cagePeaks[cage], "_Norm_", selclus, "_clusters_cluster", k, ".GO-CC.tsv",sep = ''),
        sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
}


