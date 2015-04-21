project.init2("chipmentation")
psa = loadPSA()

# Look for combined samples
psa[biologicalReplicate==0 & technicalReplicate==0,]


# Divide the CAGE peaks into TATA and Cpg Island groups:
tss = loadCageTSS();
tss
tssGroupIds = split(1:nrow(tss), tss$group)
names(tssGroupIds)
tssGroupIds$TATA.CpG

msa = psa[biologicalReplicate==0 & technicalReplicate==0,]
psa[1,]

pdf(paste0(resDir, "tss400.pdf"))
for (i in 1:nrow(msa)) {
	message(msa[i, sampleName])
	simpleCache(msa[i, sampleName], cacheSubDir="tss400bp")
	mat = get(msa[i, sampleName])[tssGroupIds$TATA.CpG,]
	plot(apply(mat, 2, sum), type="l", main = msa[i, sampleName])
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












