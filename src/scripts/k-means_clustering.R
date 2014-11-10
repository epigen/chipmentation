library(gplots)
require(made4)

# read data
#df = read.table("ann_inH3K4me3peak_K562_DPIcluster_hg19_20120116.permissive_set_ALL.bed_ghist.txt",
#df = read.table("ann_inH3K4me3peak_K562_DPIcluster_hg19_20120116.permissive_set_ALL.bed_ghist.txt",
df = read.table("/fhgfs/groups/lab_bock/shared/data/ChIPmentation_Freeze/homer/annotatePeaks/ann_inH3K4me3peak_K562_DPIcluster_hg19_20120116.permissive_set_ALL.bed_ghist.txt",
    sep='\t', header=TRUE)

colnames(df) <- c("peak", seq(-75, 75))
rownames(df) = df[,1]
df = df [,-1]

# align TSSs
#minus = df[grepl('-', rownames(df)),]
#plus = df[!grepl('-', rownames(df)),]
#minusRev = as.data.frame(t(apply(minus, 1, rev)))
#colnames(minusRev) = colnames(plus)

#dfInv = rbind(plus, minusRev)

# feature scaling
#dfNorm = t(apply(df, 1, function(x)(x - mean(x) / sd(x)))) # scale to [-1:1]
dfNorm = t(apply(df, 1, function(x)(x - min(x)) / (max(x) - min(x)))) # scale to [0:1]
dfNorm = dfNorm[complete.cases(dfNorm),]


# plot heatmap
dist.pear <- function(x) as.dist(1-cor(t(x)))
hclust.ave <- function(x) hclust(x, method="average")

s = sample(1:nrow(dfNorm), 200)
heatplot(as.matrix(dfNorm[s,]), dend = 'row', labRow = NA, labCol = NA)

heatmap.2(as.matrix(dfNorm[s,]),
    Colv = NA,
    trace = "none",
    labRow = NA, labCol = NA, dendrogram=c("row"), scale = c('row'), symbreak=TRUE, col=bluered(256),
    distfun=dist.pear, hclustfun=hclust.ave
)

#### K-MEANS CLUSTERING

# Rule of thumb for number of k
# k = sqrt(nrow(dfNorm)/2)

kmax = 16 # the maximum number of clusters we will examine; you can change this
totwss = rep(0, kmax) # will be filled with total sum of within group sum squares
kmfit = list() # create and empty list
for (i in 1:kmax){
    kclus = kmeans(dfNorm, centers = i, iter.max = 20)
    totwss[i] = kclus$tot.withinss
    kmfit[[i]] = kclus

    print(i)#, kclus$tot.withinss)
}

####################################################################
# calculate the adjusted R-squared
# and plot it
####################################################################
n = nrow(dfNorm)
rsq = 1-(totwss*(n-1))/(totwss[1]*(n-seq(1,kmax)))

kmeansAIC = function(fit){
    m = ncol(fit$centers)
    n = length(fit$cluster)
    k = nrow(fit$centers)
    D = fit$tot.withinss
    return(D + 2*m*k)
}

aic = sapply(kmfit, kmeansAIC)


pdf("/home/arendeiro/projects/chipmentation/results/plots/H3K4me3_tss_k-means_clustering.pdf")
require("sfsmisc")
mult.fig(2)
plot(seq(1,kmax), rsq, xlab = "Number of clusters", ylab = "Adjusted R-squared", pch=20, cex=2)
plot(seq(1,kmax), aic, xlab = "Number of clusters", ylab = "AIC", pch=20, cex=2)
dev.off()

####################################################################
# try to locate the "elbow" point
####################################################################
v = -diff(aic)
nv = length(v)
fom = v[1:(nv-1)]/v[2:nv]
nclus = which.max(fom)+1
cat("The apparent number of clusters is: ",nclus,"\n")
#points(nclus,aic[nclus],col=2,pch=20,cex=2)


### Export with cluster annotation for all numbers of clusters
for (k in 1:kmax){
    nclus = k
    clustered = cbind(rownames(dfNorm), dfNorm, kmfit[[nclus]]$cluster)
    colnames(clustered) <- c("peak", seq(-75, 75), "cluster")

    write.table(clustered, paste("tss_peaks_cage_", nclus, "clusters.tsv",sep = ''), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # all

    for (clus in 1:nclus){
        write.table(clustered[clustered$cluster == clus , ], paste("tss_peaks_cage_", nclus, "clusters_cluster", clus, ".tsv", sep = ''),
        sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE) # 1
    }
}


# Plot clusters
pdf("clusters_16.pdf")#/home/arendeiro/projects/chipmentation/results/plots/H3K4me3_tss_k-means_clustering_averageProfile.pdf")
colors = rainbow(kmax)
m <- matrix(c(seq(1,16), rep(17,4)),nrow = 5,ncol = 4,byrow = TRUE)
layout(mat = m,heights = c(rep(0.4, 8), 0.2))

for (k in 1:kmax){
    nclus = k
    clustered = cbind(dfNorm, kmfit[[k]]$cluster)
    colnames(clustered) <- c(seq(-75, 75), "cluster")
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
                ylim = c(0, 280)
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
legend(x = "top", inset = 0, legend = 1:kmax, col = colors, bg = "white", lwd = 2, ncol = 4,horiz = FALSE)

dev.off()


# Plot heatmap with clusters
s = sample(1:nrow(dfNorm), 200)
heatplot(as.matrix(dfNorm[s,]), dend = 'row', labRow = NA, labCol = NA, classvec = kmfit[[16]]$cluster[s])




## GO of clusters
library(biomaRt)
library(org.Hs.eg.db)
library(topGO)

# export genomic locations of peaks
chr = apply(df, 1, function(x) strsplit(as.character(x[1]), ':')[[1]][1])
start = apply(df, 1, function(x) strsplit(strsplit(as.character(x[1]), ':')[[1]][2], '..', fixed=TRUE)[[1]][1])
end = apply(df, 1, function(x) strsplit(strsplit(strsplit(as.character(x[1]), ':')[[1]][2], '..', fixed=TRUE)[[1]][2], ',')[[1]][1])
strand = apply(df, 1, function(x) strsplit(strsplit(as.character(x[1]), ':')[[1]][2], ',')[[1]][2])

write.table(cbind(chr,start,end,strand), 
    'ann_inH3K4me3peak_K562_DPIcluster_hg19_20120116.permissive_set_ALL.bed_ghist.bed', 
    sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

#bedtools closest -t first \
#-a ann_inH3K4me3peak_K562_DPIcluster_hg19_20120116.permissive_set_ALL.bed_ghist.bed \
#-b reference/Homo_sapiens/GRCh37_hg19_refSeq.tss.bed \
#| cut -f 8 \
#>  ann_inH3K4me3peak_K562_DPIcluster_hg19_20120116.permissive_set_ALL.bed_ghist.genes.bed

genes = read.table("ann_inH3K4me3peak_K562_DPIcluster_hg19_20120116.permissive_set_ALL.bed_ghist.genes.bed", sep = '\t', header = FALSE)

# get relation between refseq and ensembl ids
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
refseq <- genes[,1]
#listAttributes(mart)[grep('RefSeq', listAttributes(mart)$description),]
ids = getBM(attributes = c("refseq_mrna", "refseq_ncrna", 'ensembl_gene_id'), values = genes[,1], mart = mart)

# annotate with ensembl ids
genesNR = data.frame(genes[grep('NR', as.character(genes[,1])),])
colnames(genesNR) = 'refseq_ncrna'
genesE = merge(genesNR, ids)$ensembl_gene_id
genesNM = data.frame(genes[grep('NM', as.character(genes[,1])),])
colnames(genesNM) = 'refseq_mrna'
genesE = c(genesE, merge(genesNM, ids)$ensembl_gene_id)

allGenes = rep(0, length(genesE))
names(allGenes) <- genesE

kmax = 16 # use case with 16 clusters
testResults = list()
for (k in 1:kmax){
    genesCluster = allGenes[kmfit[[kmax]]$cluster == k]
    GOdata <- new("topGOdata", ontology = "BP", allGenes = genesCluster, geneSel = function(p) p < 1e-2, 
        description = "Test", annot = annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")

    resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    testResults[[k]] = GenTable(GOdata, classicFisher = resultFisher, topNodes = 100)
}

