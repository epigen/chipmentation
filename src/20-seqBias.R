# Addressing sequencing bias:

project.init("chipmentation", "projects/chipmentation")
eload(loadPSA())
eload(loadCageTSS());
SV$tss
tssGR = dtToGr(SV$tss, "V1", "V2", "V3", strand="V6")
tssGR = resize(tssGR, 500, fix='center')

rp()
eload(loadTransposePWM())
tssGR[1]
SV$tss
#s = extractSequences(tssGR[1]) #broken
tssSequences = getSeq(Hsapiens, tssGR[1:50000])


# getSeq correctly accounts for strand:
testRegionNegStrand = GRanges(seqnames="chr1", IRanges(start=100013159, end=100013658), strand="-")
testRegionPosStrand = GRanges(seqnames="chr1", IRanges(start=100013159, end=100013658), strand="+")
testRegionNegStrand; testRegionPosStrand;
getSeq(Hsapiens, testRegionNegStrand)
reverseComplement(getSeq(Hsapiens, testRegionPosStrand))


simpleCache("tssSequences_500s", { getSeq(Hsapiens, tssGR) })

utility("funcMotif.R")
tssSequences_500

tssConsensus = (consensusMatrix(tssSequences_500, as.prob=TRUE)[DNA_BASES,])

#' A convenience function to stick the file into the results directory
#' (resDir).
pdfrd = function(file, ...) {
	file = paste0(resDir, file)
	pdf(file, ...)
}

#pdf(paste0(resDir, "tssSignals/tss_logo.pdf"), width=16)
pdfrd("tssSignals/tss_logo.pdf", width=16)
seqLogo(makePWM(tssConsensus[,200:300]), ic.scale=FALSE)
dev.off()

tssSequences_500
setLapplyAlias(8)
pwmScores = lapplyAlias(as.character(tssSequences_500), matchPWMscoreI, SV$tpwm)
SV$tpwm
pwmScoreMatrix = do.call(rbind, pwmScores)
scale(colSums(pwmScoreMatrix))
pdfrd("tssSignals/tss_tpswm_score.pdf", width=16)
plot(-240:239, scale(colSums(pwmScoreMatrix)), type="b")
plot(-90:110, scale(colSums(pwmScoreMatrix[,150:350])), type="b", lwd=2)
lines(-50:50, scalemod, col="blue", lwd=2)
dev.off()




