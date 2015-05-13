# Addressing sequencing bias:
project.init("chipmentation", "chipmentation")
utility("funcMotif.R")
library(seqLogo)
eload(loadPSA())
eload(loadCageTSS());
eload(loadTransposePWM())
utility("funcMotif.R")
SV$tss

eload(loadPeaks())

# TSS ###############################
# Get sequences:
Hsapiens.masked = BSgenome.Hsapiens.UCSC.hg19.masked

simpleCache("tssSeq", { getSeq(Hsapiens , SV$tssGR) }, assignToVariable="seqs")
simpleCache("tssSeqMasked", { extractSequences(SV$tssGR, Hsapiens.masked); }, assignToVariable="seqMasked", recreate=TRUE)

seqs

# Calculate PWM matches: 
setLapplyAlias(8)
pwmScoreMatrix = pwmMatchSignal(seqs, SV$tpwm)

subsetIds = SV$tssGroupIds[["TATA.CpG.Exp"]]

subsetIds = SV$tssGroupIds[["TATA-less.CpG-less.Exp"]]

preamble = substitute({
project.init2("chipmentation")
eload(loadPSA())
eload(loadTransposePWM())
eload(loadCageTSS());
eload(loadPeaks());
Hsapiens.masked = BSgenome.Hsapiens.UCSC.hg19.masked
utility("funcMotif.R")
setLapplyAlias(4)
})


# Grab signal data:


# tss subsets:

simpleCache("combinedCM_tss", { mergeCachedMatrices(paste0(SV$msa[SV$cmIds, sampleName], "_cap5"), "signal/tss_cap5/") } )
simpleCache("combinedCM_tss_igg", { mergeCachedMatrices(paste0(SV$msa[SV$cmIggIds, sampleName], "_cap5"), "signal/tss_cap5/") } )
simpleCache("K562_500K_CM_H3K4ME3_nan_nan_0_0_hg19_cap5", cacheSubDir=paste0("signal/tss_cap5/") , assignToVariable="cmSig") # ChIPmentation sample
simpleCache("K562_500K_ATAC_INPUT_nan_01ULTN5_PE_1_1_hg19_cap5", cacheSubDir=paste0("signal/tss_cap5/") , assignToVariable="nakedSig") # naked dna

pdfrd("tssSignals/transposase_bias_TSS.pdf", width=16)
iType = 2;
for (iType in 1:length(SV$tssGroupIds)) {
	type = names(SV$tssGroupIds)[iType]
	message(type)
	subsetIds = SV$tssGroupIds[[type]]
	sbo = seqBias(paste0("tss_", type), SV$tssGR[subsetIds], SV$tpwm, Hsapiens, Hsapiens.masked)
	m = summarizeMatrixModel(combinedCM_tss[subsetIds,])
	m_bg = summarizeMatrixModel(combinedCM_tss_igg[subsetIds,])
	m_nakedSig = summarizeMatrixModel(nakedSig[subsetIds,])
	m_cmSig = summarizeMatrixModel(cmSig[subsetIds,])
	with(sbo, seqBiasPlot(factorName, models, modelsMasked, list(CM_Combined=m, CM_IGG=m_bg, CM_H3K4me3=m_cmSig, ATAC_INPUT=m_nakedSig)))
	lines(gscale(m_nakedSig), m_nakedSig, col="magenta", lty="dashed")
}
dev.off()

cor(m_nakedSig, m_cmSig)


# FACTORS
SV$factors
pdfrd("tssSignals/transposase_bias_factor.pdf", width=16)
iFactor=1
for (iFactor in 1:length(SV$factors)) { 
	factorName = SV$factors[[iFactor]];
	message(factorName);
	simpleCache(paste0("combinedCM_", factorName), { mergeCachedMatrices(paste0(SV$msa[SV$cmIds, sampleName], "_cap5"), paste0("signal/", factorName, "_cap5/")) } , assignToVariable="combinedCM")

	simpleCache(paste0("combinedCM_", factorName, "_igg"),  { mergeCachedMatrices(paste0(SV$msa[SV$cmIggIds, sampleName], "_cap5"), paste0("signal/", factorName, "_cap5/")) } , assignToVariable="combinedCM_igg")

 	simpleCache("K562_500K_CM_H3K4ME3_nan_nan_0_0_hg19_cap5", cacheSubDir=paste0("signal/", factorName, "_cap5/") , assignToVariable="cmSig") # ChIPmentation sample
simpleCache("K562_500K_ATAC_INPUT_nan_01ULTN5_PE_1_1_hg19_cap5", cacheSubDir=paste0("signal/", factorName, "_cap5/") , assignToVariable="nakedSig") # naked dna

	SV[[factorName]] = resize(SV[[factorName]], 200, fix="center")
	sbo = seqBias(factorName, SV[[factorName]], SV$tpwm, Hsapiens, Hsapiens.masked)
	m = summarizeMatrixModel(combinedCM)
	m_bg = summarizeMatrixModel(combinedCM_igg)
	m_nakedSig = summarizeMatrixModel(nakedSig)
	m_cmSig = summarizeMatrixModel(cmSig)
	with(sbo, seqBiasPlot(factorName, models, modelsMasked, list(CM_Combined=m, CM_IGG=m_bg, CM_H3K4me3=m_cmSig, ATAC_INPUT=m_nakedSig)))
}
dev.off()

# Produce text files;
for (iFactor in 1:length(SV$factors)) { 
	factorName = SV$factors[[iFactor]];
	message(factorName)
	sbo = seqBias(factorName, SV[[factorName]], SV$tpwm, Hsapiens, Hsapiens.masked)
	sbo
write(sbo$models$modelPWM, rd(paste0("sequence_models/", factorName, "_tranpososePWM.txt")), ncol=1)
write(scale(sbo$models$modelDinuc), rd(paste0("sequence_models/", factorName, "_ATDinuc.txt")), ncol=1)
}


# TODO: 
# - add naked DNA control to TF plots
# - provide PWM and AT signal tracks to Andre
# - produce TSS plot (all TSS, just H3K4me3 signal), with background plots:
# (AT, PWM, and Naked DNA). To show it's not a bias.
# NEXT:
# - what is up with TSS signal? Cluster it?



################################################################################
# Some previous work, looking at consensus motifs, which I am 
# no longer using:

for i in 
simpleCache("ctcfSeq", { getSeq(BSgenome.Hsapiens.UCSC.hg19.masked, SV$ctcf) }, cacheSubDir="sequence")
tssConsensus = (consensusMatrix(ctcfSeq, as.prob=TRUE)[DNA_BASES,])
	seqLogoLabeled(makePWM(tssConsensus), ic.scale=FALSE, main=type)
}


#s = extractSequences(tssGR[1]) #broken
#tssSequences = getSeq(Hsapiens, tssGR[1:50000])

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


################################################################################
# Some slurm ideas (not really used):
dir.create("slurm")
buildSlurmScript(rcode, preamble, submit=TRUE)
eval(preamble)


