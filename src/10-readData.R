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

i=1

pdf(paste0(resDir, "tss400_raw.pdf"), width=16)
for (i in 1:nrow(SV$msa)) {
	message(i, ": ", SV$msa[i, sampleName])
	simpleCache(paste0(SV$msa[i, sampleName]), cacheSubDir="tss400bp", reload=FALSE, assignToVariable="mat")

	par(mfrow=c(2,4))
	for (type in names(tssGroupIds)) {
	plot(-199:200, apply(mat[tssGroupIds[[type]],], 2, sum), type="l", main = paste0(SV$msa[i, sampleName]), xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	plot(-50:50, apply(mat[tssGroupIds[[type]],150:250], 2, sum), type="l", xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	} # for
}
dev.off()


pdf(paste0(resDir, "tss400_cap.pdf"), width=16)
for (i in 1:nrow(SV$msa)) {
	message(i, ": ", SV$msa[i, sampleName])
	simpleCache(paste0(SV$msa[i, sampleName], "_cap"), cacheSubDir="tss400bp_capped", reload=FALSE, assignToVariable="mat")

	par(mfrow=c(2,4))
	for (type in names(tssGroupIds)) {
	plot(-199:200, apply(mat[tssGroupIds[[type]],], 2, sum), type="l", main = paste0(SV$msa[i, sampleName]), xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	plot(-50:50, apply(mat[tssGroupIds[[type]],150:250], 2, sum), type="l", xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	} # for
}
dev.off()


pdf(paste0(resDir, "tss400_bin.pdf"), width=16)
for (i in 1:nrow(SV$msa)) {
	message(i, ": ", SV$msa[i, sampleName])
	simpleCache(paste0(SV$msa[i, sampleName], "_bin"), cacheSubDir="tss400bp_bin", reload=TRUE, assignToVariable="mat")
	par(mfrow=c(2,4))
	for (type in names(tssGroupIds)) {
	plot(-199:200, apply(mat[tssGroupIds[[type]],], 2, sum), type="l", main = paste0(SV$msa[i, sampleName]), xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	plot(-50:50, apply(mat[tssGroupIds[[type]],150:250], 2, sum), type="l", xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	} # for
}
dev.off()



# Now, can I cluster the ones that do or do not contribute to the pattern?
# Or, predict TSS?

# A good example is: K562_10M_CM_H3K27AC_nan_nan_0_0_hg19
i = 8
SV$msa[8,sampleName]
simpleCache(paste0(SV$msa[i, sampleName], "_bin"), cacheSubDir="tss400bp_bin", reload=TRUE, assignToVariable="mat")
mat
mod = apply(mat[tssGroupIds[["TATA-less.CpG-less"]],150:250], 2, sum)

subMat = mat[tssGroupIds[["TATA-less.CpG-less"]],150:250]

plot(mod, type="l")
scalemod = scale(mod)[,1]
plot(scalemod, type="l")
sum(scalemod * subMat[1,])
plot(scalemod, type="l")
for(i in showme) {
	lines(subMat[i,], type='l', col=i)
}
zscore(mod)

c(4,5,6) * 
matrix(c(2,2,2,4,4,4), ncol=3)
#sweep is right. 
sts = sweep(subMat, MARGIN=2, scalemod, `*`)

dim(sts)
sum(sts[1,])
sum(subMat[1,] * scalemod)

stsSum=apply(sts, 1, sum)
hitCount = apply(subMat, 1, sum)
score = stsSum/hitCount
showme = which(score > 2)
i = 118
plot(density(na.omit(score)))
abline(v=0)

table(hitCount)






#COMBINE ALL CM DATA:


SV$msa[technique=="CM" & substr(ip, 1, 2) == "H3" & numberCells != "1K",]
cmIds = SV$msa[technique=="CM" & substr(ip, 1, 2) == "H3" & numberCells != "1K", which=TRUE]

SV$msa[technique=="CM",]
cmIds = SV$msa[technique=="CM", which=TRUE]

allCM = matrix(0, nrow=184827, ncol=400)
for (i in cmIds) {
	message(i, ": ", SV$msa[i, sampleName])
	simpleCache(paste0(SV$msa[i, sampleName], "_bin"), cacheSubDir="tss400bp_bin", reload=TRUE, assignToVariable="mat")
	allCM = allCM+mat
}
mat =allCM

pdf(paste0(resDir, "tssSignals/tss400_bin_cm_combined.pdf"), width=16)
par(mfrow=c(2,4))
	for (type in names(tssGroupIds)) {
	plot(-199:200, apply(mat[tssGroupIds[[type]],], 2, sum), type="l", main = paste0("CM combined"), xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	plot(-50:50, apply(mat[tssGroupIds[[type]],150:250], 2, sum), type="l", xlab="TSS")
	legend('topleft', type, bty='n', cex=.8)
	abline(v=0, col="gray", lty="dotted")
	} # for
dev.off()
















# DEPRECATED:
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












