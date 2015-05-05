# Addressing sequencing bias:

project.init("chipmentation", "projects/chipmentation")
utility("funcMotif.R")
library(seqLogo)
eload(loadPSA())
eload(loadCageTSS());
eload(loadTransposePWM())
SV$tss
tssGR = dtToGr(SV$tss, "V1", "V2", "V3", strand="V6")
tssGR = resize(tssGR, 500, fix='center')

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



tssSequences_500
setLapplyAlias(8)
pwmScores = lapplyAlias(as.character(tssSequences_500), matchPWMscoreI, SV$tpwm)
pwmScoreMatrix = do.call(rbind, pwmScores)
scale(colSums(pwmScoreMatrix))
pdfrd("tssSignals/tss_tpswm_score2.pdf", width=16)
plot(-240:239, scale(colSums(pwmScoreMatrix)), type="b")
plot(-90:110, scale(colSums(pwmScoreMatrix[,150:350])), type="b", lwd=2)
lines(-50:50, scalemod, col="blue", lwd=2)
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

