# Addressing sequencing bias:
project.init("chipmentation", "chipmentation")
utility("funcMotif.R")
library(seqLogo)
eload(loadPSA())
eload(loadCageTSS());
eload(loadTransposePWM())
SV$tss

eload(loadPeaks())


for i in 
simpleCache("ctcfSeq", { getSeq(BSgenome.Hsapiens.UCSC.hg19.masked, SV$ctcf) }, cacheSubDir="sequence")
tssConsensus = (consensusMatrix(ctcfSeq, as.prob=TRUE)[DNA_BASES,])
	seqLogoLabeled(makePWM(tssConsensus), ic.scale=FALSE, main=type)
}


#s = extractSequences(tssGR[1]) #broken
#tssSequences = getSeq(Hsapiens, tssGR[1:50000])

# getSeq correctly accounts for strand:
testRegionNegStrand = GRanges(seqnames="chr1", IRanges(start=100013159, end=100013658), strand="-")
testRegionPosStrand = GRanges(seqnames="chr1", IRanges(start=100013159, end=100013658), strand="+")
testRegionNegStrand; testRegionPosStrand;
getSeq(Hsapiens, testRegionNegStrand)
reverseComplement(getSeq(Hsapiens, testRegionPosStrand))


simpleCache("tssSequences_500s", { getSeq(Hsapiens, tssGR) }, assignToVariable="tssSequences_500")


tssSequences_500

# For all of them.
tssConsensus = (consensusMatrix(tssSequences_500, as.prob=TRUE)[DNA_BASES,])
a = seqLogo(makePWM(tssConsensus[,200:300]), ic.scale=FALSE)

pdfrd("tssSignals/tss_logo_divided.pdf", width=16)
for (type in names(SV$tssGroupIdsNoExp)) {
	message(type);
	subsetIds = SV$tssGroupIdsNoExp[[type]]
	tssConsensus = (consensusMatrix(tssSequences_500[subsetIds,], as.prob=TRUE)[DNA_BASES,])
	seqLogoLabeled(makePWM(tssConsensus[,200:300]), ic.scale=FALSE, main=type)
}
dev.off()



pdfrd("tssSignals/tss_tpswm_score_divided.pdf", width=16)
for (type in names(SV$tssGroupIdsNoExp)) {
	message(type);
	subsetIds = SV$tssGroupIdsNoExp[[type]]
	tssSeqSubset = tssSequences_500[subsetIds,]
	pwmScores = lapplyAlias(as.character(tssSeqSubset), matchPWMscoreI, SV$tpwm)
	pwmScoreMatrix = do.call(rbind, pwmScores)
plot(-240:239, scale(colSums(pwmScoreMatrix)), type="l", main=type)
plot(-90:110, smooth.spline(scale(colSums(pwmScoreMatrix[,150:350])))$y, type="l", lwd=2)
lines(-50:50, smooth.spline(tssCMModels[[paste0(type, ".Exp")]])$y, col="darkgreen", lwd=1)
lines(-50:50, smooth.spline(tssCMModels[[paste0(type, ".NoExp")]])$y, col="maroon", lwd=1)
}
dev.off()

sv()


# TSS ###############################
# Get sequences:
Hsapiens.masked = BSgenome.Hsapiens.UCSC.hg19.masked

simpleCache("tssSeq", { getSeq(Hsapiens , SV$tssGR) }, assignToVariable="seqs")
simpleCache("tssSeqMasked", { extractSequences(SV$tssGR, Hsapiens.masked); }, assignToVariable="seqMasked", recreate=TRUE)

seqs

# Calculate PWM matches: 
setLapplyAlias(8)
pwmScoreMatrix = pwmMatchSignal(seqs, SV$tpwm)

# Grab signal data:
simpleCache("combinedCMData", { mergeCachedMatrices(paste0(SV$msa[cmIds, sampleName], "_cap5"), "signal/tss_cap5/") } )

dim(combinedCMData)
m = summarizeMatrixModel(combinedCMData)

# Plot them on the same axis:

pdfrd("tssSignals/tss_tpswm_score2.pdf", width=16)

dev.off()


# Subset if desired:
for (type in names(SV$tssGroupIds)) {
	# add expression condition:
	#type = paste0(type, ".NoExp")
	subsetIds = SV$tssGroupIds[[type]]
}



subsetIds = SV$tssGroupIds[["TATA.CpG.Exp"]]

subsetIds = SV$tssGroupIds[["TATA-less.CpG-less.Exp"]]

seqBias("tss_TATA-less.CpG-less.Exp", SV$tssGR[subsetIds], combinedCMData[subsetIds,])
 
# and a function to do the same:
factorName="tss"
seqBias = function(factorName, GR, signalData) {
	message("Get sequences...")
	simpleCache(paste0(factorName, "Seq"), { getSeq(Hsapiens , GR) }, assignToVariable="seqs", buildEnvir=nlist(GR))
	simpleCache(paste0(factorName, "SeqMasked"), { extractSequences(GR, Hsapiens.masked) }, assignToVariable="seqsMasked", buildEnvir=nlist(GR))
	# TODO: extractSequences needs to control for strand.

	message("Calculate PWM matches...")
	setLapplyAlias(8)
	simpleCache(paste0(factorName, "PWM"), {
		pwmMatchSignal(seqs, SV$tpwm)
	}, buildEnvir=nlist(seqs), assignToVariable="pwmScoreMatrix")
	simpleCache(paste0(factorName, "PWMMasked"), {
		pwmMatchSignal(seqsMasked, SV$tpwm)
	}, buildEnvir=nlist(seqsMasked), assignToVariable="pwmScoreMatrixMasked")

	m = summarizeMatrixModel(signalData)
	plot(gscale(pwmScoreMatrix), pwmScoreMatrix, type="l", ylim=c(-3, 3))
	lines(gscale(pwmScoreMatrixMasked), pwmScoreMatrixMasked, type="l", col="red")
	lines(gscale(m), m, type="l", main = paste0("CM combined"), xlab="TSS", col="blue")
	legend('topleft', c("PWM match (rpts)", "PWM match (masked)", "CM signal"), col=c("black", "red", "blue"), lty=1)
	crosshair()
}
# TODO: Now run this on every category of interest...






# CTCF: 
simpleCache("mat_ctcf", { mergeCachedMatrices(paste0(SV$msa[cmIds, sampleName], "_cap5"), "signal/ctcf_cap5/") } )
mod_ctcf = summarizeMatrixModel(mat_ctcf)
lines(mod_ctcf, type='l')
 



sapply(c(100,101,102,105), function(x) { length(gscale(x)) } )





