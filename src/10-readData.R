project.init2("chipmentation")
psa = loadPSA()

# Look for combined samples

psa[biologicalReplicate==0 & technicalReplicate==0,]

tss = loadCageTSS();
# Divide the CAGE peaks into TATA and Cpg Island groups:

tss[, group:=interaction(V10, V11)]
tss
bedGR = dtToGr(tss, "V1", "V2", "V3", strand="V6", splitFactor="group")
bedGR=GRangesList(bedGR)

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





#Now take a look at exact cuts:
psa[,xcfile:= paste0("media/nsheffield/red6/chipmentation/data/shiftedExactCuts/", sampleName, ".bigWig")]
psa[1,]


am = as.matrix(a[, -1, with=FALSE])
dim(am)

plot(apply(am, 2, sum), type="l")



