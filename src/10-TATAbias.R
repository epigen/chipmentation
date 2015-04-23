#GOAL: Check for bias in tagmentation. We noticed signal enrichment at TATA boxes; is it just because of sequence bias toward TATA?
#Strategy: Check the signal at TATAs in boxes and compare to the signal at TATAs outside of boxes.

project.init("chipmentation", "projects/chipmentation")
library(rtracklayer)

### 1. Read in CAGE peaks and identify all TATAs.
loadBSgenome("hg19")

bed = fread("../data/hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed.TATA.bed", sep="\t")
setnames(bed, c("chr", "start", "stop", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"))

seqLength = 50
bedGR= promoters(dtToGr(bed, strand="strand"), upstream=seqLength, downstream=0)

bedGR[strand(bedGR)=="+"]=shift(bedGR[strand(bedGR)=="+"],-5)
bedGR[strand(bedGR)=="-"]=shift(bedGR[strand(bedGR)=="-"],5)

export.bed(bedGR, "bed/promoters.bed")

seq = extractSequences(bedGR, mask=FALSE)
seq[strand(bedGR)=="-"] = reverseComplement(seq)[strand(bedGR)=="-"]

v = vmatchPattern("TATA", seq, max.mismatch=1); 

s= unlist(lapply(v, start))
sl = lapply(v, start)
sExpanded = as.vector(sapply(s, function(x) { x+c(0:3); }))
par(mfrow=c(2,1))
barplot(table(s))
barplot(table(sExpanded))
table(sapply(v, length))
boxTATAs = as.integer(names(which( table(s)/length(s) > .01)))
as.integer(names(which.max(table(sExpanded)))) #center:
boxTATAs = as.integer(names(which.max(table(sExpanded)))) + c(-5:5)
boxTATAs
#nonboxTATAs = (1:300)[-boxTATAs]
#lapply(sl, function(x) { x %in% boxTATAs; })

# Build a DT of tata coordinates in each seq
DT=data.table("regionID"=numeric(0), "TATAcoord"=numeric(0), "origin"=logical(0))
for (i in 1:length(sl)) {
	if(length(sl[[i]]) == 0) next;
	DT=rbind(DT, 	cbind(regionID=i, TATAcoord=sl[[i]]+1, origin=TRUE))
	DT=rbind(DT, 	cbind(regionID=i, TATAcoord=sl[[i]]+2, origin=FALSE))
	DT=rbind(DT, 	cbind(regionID=i, TATAcoord=sl[[i]]+3, origin=FALSE))
	DT=rbind(DT, 	cbind(regionID=i, TATAcoord=sl[[i]]+4, origin=FALSE))
}
#restrict to within sequence

max(DT$TATAcoord)
DT = DT [TATAcoord < seqLength,]
DT
#annotate TATAs as in box or not.
DT[, isBox:=TATAcoord %in% boxTATAs]

pdf("fig/TATAlocations.pdf")
par(mfrow=c(1,1))
#plot(as.matrix(DT[,list(TATAcoord, regionID)]), pch=16, cex=.3, main="TATA locations")
x = DT[origin==1, TATAcoord]
y = DT[origin==1, regionID]
r = cbind(rep(4, length(x)),.05)
co = DT[origin==1, isBox]
plot(as.matrix(DT[origin==1,list(TATAcoord, regionID)]), pch=16, cex=.3, main="TATA locations", type="n")
symbols(x=x, y=y, rectangles=r, add=TRUE, inches=FALSE)
abline(v=max(boxTATAs), col="red", lwd=4)
abline(v=min(boxTATAs), col="red", lwd=4)
dev.off()

### 2. Load up the signal data

# bwSlurp("/fhgfs/groups/lab_bock/arendeiro/share/chipmentationBedFiles/H3K4me3_K562_500k_ChIP.5prime.bw", chr="chr6", start=36561943, end=36562142)

# Read bamfiles
bam = bamSlurp("/fhgfs/groups/lab_bock/arendeiro/share/chipmentationBamFiles/H3K4me3_K562_500k_CM.bam", bedGR)
#bamSlurpPlot(bam[1:5], bedGR[1:5])

#Convert bam data into vectors
bamVectors = bamSlurpVectors(bam, bedGR)
bamVec = lapply(bamVectors, function(x) { as.numeric(rep(as.vector(runValue(x)), as.vector(runLength(x)))) })

#combine vectors into a matrix
bamMat = do.call(rbind, bamVec)
class(bamMat); dim(bamMat); #apply(bamMat, 2, sum)
#Take a look at aggregate signal:

pdf("fig/line_signal.pdf")
plot(apply(bamMat > 0, 2, sum), type="l", main="Total Signal")
dev.off()

### 3. Combine signal data with TATA location information.

#Extract signal values at TATA spots
DT[, signalValue:=bamMat[as.matrix(DT[,list(regionID, TATAcoord)])]]
DT[, signalValue:=capData(signalValue)]

pdf("fig/vio_TATAsignal.pdf")
compareSignals(list(box=DT[isBox==TRUE, signalValue], no_box=DT[isBox==FALSE, signalValue]))
dev.off()

DTsummary = DT[,.N, by=list(isBox, signalValue)]
DTsummary[, totalCount:=sum(N), by=list(isBox)]
DTsummary[, proportion := N/totalCount]
DTsummary
DT[, list(hit=sum(signalValue > 0), count=.N), by=isBox]
DT[, list(hit=sum(signalValue > 0)/.N), by=isBox]
fisher.test(as.matrix(DT[, list(hit=sum(signalValue > 0), count=.N), by=isBox][, list(hit, count)]))

#ggplot(DT, aes(x=signalValue, group=isBox)) + geom_histogram() + theme_classic() +facet_grid(isBox ~ .)
#ggplot(DT, aes(x=signalValue, fill=isBox)) + geom_histogram(alpha=0.8, position="dodge", binwidth=1)+ theme_classic() + scale_x_discrete()

pdf("fig/bar_TATAsignal.pdf")
DTsummary
ggplot(DTsummary, aes(y=proportion, x=factor(signalValue), fill=isBox)) + geom_bar(color="black", stat="identity", position="dodge") + scale_x_discrete() + theme_classic() + theme(aspect.ratio=1)
dev.off()

pdf("fig/bar_TATA_binary.pdf", height=3);
par(mfrow=c(2,1))
labeledBarplot(as.matrix(c(.217, 1-.217)), "No box", col=c("yellow", "darkgray"));
labeledBarplot(as.matrix(c(.351, 1-.351)), "TATA box", col=c("yellow", "darkgray"));
dev.off()


