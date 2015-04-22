project.init("chipmentation", "projects/chipmentation")
#psa = loadPSA()
svloadenv(loadPSA())

# Look for combined samples

# Divide the CAGE peaks into TATA and Cpg Island groups:
tss = loadCageTSS();
tss
tssGroupIds = split(1:nrow(tss), tss$group)
names(tssGroupIds)
tssGroupIds$TATA.CpG

msa = psa[biologicalReplicate==0 | technicalReplicate==0,]


psa

tss

pdf(paste0(resDir, "tss400_raw.pdf"), width=16)
for (i in 1:nrow(msa)) {
	message(i, ": ", msa[i, sampleName])
	simpleCache(paste0(msa[i, sampleName]), cacheSubDir="tss400bp", reload=FALSE, assignToVariable="mat")

	par(mfrow=c(2,4))
	for (type in names(tssGroupIds)) {
	plot(-199:200, apply(mat[tssGroupIds[[type]],], 2, sum), type="l", main = paste0(msa[i, sampleName]), xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	plot(-50:50, apply(mat[tssGroupIds[[type]],150:250], 2, sum), type="l", xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	} # for
}
dev.off()


pdf(paste0(resDir, "tss400_cap.pdf"), width=16)
for (i in 1:nrow(msa)) {
	message(i, ": ", msa[i, sampleName])
	simpleCache(paste0(msa[i, sampleName], "_cap"), cacheSubDir="tss400bp_capped", reload=FALSE, assignToVariable="mat")

	par(mfrow=c(2,4))
	for (type in names(tssGroupIds)) {
	plot(-199:200, apply(mat[tssGroupIds[[type]],], 2, sum), type="l", main = paste0(msa[i, sampleName]), xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	plot(-50:50, apply(mat[tssGroupIds[[type]],150:250], 2, sum), type="l", xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	} # for
}
dev.off()


pdf(paste0(resDir, "tss400_bin.pdf"), width=16)
for (i in 1:nrow(msa)) {
	message(i, ": ", msa[i, sampleName])
	simpleCache(paste0(msa[i, sampleName], "_bin"), cacheSubDir="tss400bp_bin", reload=TRUE, assignToVariable="mat")
	par(mfrow=c(2,4))
	for (type in names(tssGroupIds)) {
	plot(-199:200, apply(mat[tssGroupIds[[type]],], 2, sum), type="l", main = paste0(msa[i, sampleName]), xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	plot(-50:50, apply(mat[tssGroupIds[[type]],150:250], 2, sum), type="l", xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	} # for
}
dev.off()




# This is the way I was loading signals from the bam files directly;
# I now switched to using bigwig files, which is faster.
tss = loadCageTSS();
tss
#bedGR = dtToGr(tss, "V1", "V2", "V3", strand="V6", splitFactor="group")
#bedGR=GRangesList(bedGR)

seqLength = 50
bedGR= promoters(bedGR, upstream=100, downstream=100)
bedGR
lapply(bedGR, length)

bam = bamSlurp(psa[biologicalReplicate==0 & technicalReplicate==0 & ip=="",filePath], bedGR[["TATA.CpG-less"]])

v = bamSlurpVectors(bam, bedGR[["TATA.CpG-less"]])
# Convert from Rle to numeric vectors
bamVec = lapply(v, function(x) { as.numeric(rep(as.vector(runValue(x)), as.vector(runLength(x)))) })

bamMat = do.call(rbind, bamVec)
bamMat

plot(apply(bamMat > 0, 2, sum), type="l", main="Total Signal")


bamFile = psa[biologicalReplicate==0 & technicalReplicate==0 & ip=="",filePath]

wholeChrom <- GRangesForUCSCGenome(genome, chrom=chr)
param = ScanBamParam(which = wholeChrom, what=scanBamWhat())
x = rep(NA, length(singleChromGR))
tryCatch( {
sb = scanBam(bamFile, param=param);
sb

bedGRL = split(bedGR[["TATA.CpG-less"]], seqnames(bedGR[["TATA.CpG-less"]]))	

bamSlurpList(bamFile, bedGR[["TATA.CpG-less"]])












