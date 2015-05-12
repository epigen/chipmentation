# R Project initialization


# This is how I rsynced the raw data to my local hard drive during the corruption
# rsync -av /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/ /media/nsheffield/red6/chipmentation/data/mapped

# all data:
# rsync -av rsync -av /fhgfs/groups/lab_bock/shared/projects/chipmentation/data /media/nsheffield/red6/chipmentation/

utility("funcGenomeSignals.R")
utility("funcGenomeLocations.R")
utility("funcLoadSharedData.R")
utility("funcDNASequence.R")

loadBSgenome("hg19")

# Set up some quick alias variables to refer to project directories easily.
dataDir = paste0(getOption("PROJECT.DATA.DIR"), "/data/")
annoDir = paste0(getOption("PROJECT.DATA.DIR"), "/annotation/")
resDir = paste0(getOption("PROJECT.DATA.DIR"), "/results/")
dir.create(resDir, showWarnings=FALSE);

# Paths to communal data files
dat = list()
dat$tss = paste0(dataDir, "hg19.cage_peak_coord_robust.400bp.bed")

loadPeaks = function() {
	ctcfBedFile = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_10M_CM_CTCF_nan_nan_0_0_hg19_peaks.motifCentered.bed"
	gata1BedFile = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_10M_CM_GATA1_nan_nan_0_0_hg19_peaks.motifCentered.bed"
	pu1BedFile = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_10M_CM_PU1_nan_nan_0_0_hg19_peaks.motifCentered.bed"
	restBedFile = "/fhgfs/groups/lab_bock/shared/projects/chipmentation/data/peaks/K562_10M_CM_REST_nan_nan_0_0_hg19_peaks.motifCentered.bed"

	ctcf = readBed(ctcfBedFile)
	gata1 = readBed(gata1BedFile)
	pu1 = readBed(pu1BedFile)
	rest = readBed(restBedFile)
	
	factors = c("ctcf", "pu1", "gata1", "rest")

	return(nlist(ctcf, ctcfBedFile, gata1, gata1BedFile, pu1, pu1BedFile, rest, restBedFile, factors))
}

loadPSA = function() {
	#psa = fread("metadata/chipmentation.sample_annotation.csv")
	#psa[, sampleName := paste(cellLine, numberCells, technique, ip, patient, treatment, biologicalReplicate, technicalReplicate, sep="_")]
	psa = fread("metadata/chipmentation.replicates.annotation_sheet.csv")
	psa[, filePath := sub("/fhgfs/groups/lab_bock/shared/projects/", dataDir, filePath)]
	psa[which(!file.exists(filePath)),]
	psa[,xcfile:= paste0(dataDir, "shiftedExactCuts/", sampleName, ".bigWig")]
	psa[which(!file.exists(xcfile)),]

	# Our new algorithm to select the one best merged sample for each
# experiment type
	msa = psa[order(biologicalReplicate,technicalReplicate), .SD[1,], by=c("cellLine", "numberCells", "technique", "ip", "treatment", "genome")]
	msa = msa[-51,]


	# COMBINE ALL CM DATA:
	msa[technique=="CM" & substr(ip, 1, 2) == "H3" & numberCells != "1K",]
	#msa[technique=="CM",]
	#cmIds = msa[technique=="CM", which=TRUE]
	cmIds = msa[technique=="CM" & substr(ip, 1, 2) == "H3" & numberCells != "1K", which=TRUE]
	cmIggIds = msa[technique=="CM" & ip =="IGG", which=TRUE]
	return(nlist(psa, msa, cmIds, cmIggIds))
}


loadTransposePWM = function() {
	downloadCache("transposasePWM", "https://raw.githubusercontent.com/GreenleafLab/NucleoATAC/master/pyatac/pwm/Human2.PWM.txt")
	transposasePWM = data.frame(transposasePWM);
	rownames(transposasePWM) = DNA_BASES
	return(list(tpwm = as.matrix(transposasePWM)))
}

loadCageTSS = function() {
	tssBedFile = paste0(annoDir, "/hg19.cage_peak_coord_robust.TATA_Annotated.bed")
	tss = fread(tssBedFile, sep="\t")
	# or switch to this set and look at robust ones?
	
	#tssall = fread(paste0(annoDir, "/hg19.cage_peak_all_annotation.tsv"), sep="\t")

	#sum(tssall[["Robust set"]]>0)

	setkey(tss, "V5")

	k562tss = fread(paste0(annoDir, "/hg19.cage_peak_K562_normexpression.tsv"), sep="\t")
	setnames(k562tss, c("tss_id", "description", "gene_id", "uniprot_id", "exp1", "exp2", "exp3"))

	setkey(k562tss, "tss_id")
	k562tss[,expMean:=pmean(exp1, exp2, exp3)]
	tss[k562tss, expMean := expMean]
	tss[, expGroup := ifelse(expMean > 0.5, "Exp", "NoExp") ]
	tss[, group:=interaction(V10, V11, expGroup)]
	tss[, groupNoExp:=interaction(V10, V11)]
	tss
	tss[,.N, by=group]
	tssGroupIds = split(1:nrow(tss), tss$group)
	tssGroupIdsNoExp = split(1:nrow(tss), tss$groupNoExp)
	tssGR = dtToGr(tss, "V1", "V2", "V3", strand="V6")
	tssGR = resize(tssGR, 400, fix='center')
	return(nlist(tss, tssGroupIds, tssGroupIdsNoExp, tssBedFile, tssGR))
	#distance distribution:
	summary(diff(tss$V3))
}


makeWindowAroundTSS = function() {
	# Produce a bed file with 400bp surrounding each cage peak,
	# to extract exact cuts.

	tss = fread(paste0(annoDir, "/hg19.cage_peak_coord_robust.TATA_Annotated.bed"), sep="\t")

	tss[, V2:=pmax(0, V2-200)]
	tss[, V3:=V3+200]
	tss[, uniqueName:=1:nrow(tss)]
	tss
	tss[, c(1,2,3, 12), with=FALSE]
	tssGR = dtToGr(tss, "V1", "V2", "V3")

	write.tsv(tss[, c(1,2,3, 12), with=FALSE], file=paste0(dataDir, "/hg19.cage_peak_coord_robust.400bp.bed"), col.names=FALSE)
}


#' Given a list of caches, load them up one-by-one and just merge (sum)
#' them. 
#' Might be worth moving this into simpleCache package.
#' @param cacheNames Vector of named caches
#' @param cacheSubDir C0mmand passed to simpleCache.
mergeCachedMatrices = function(cacheNames, cacheSubDir) {
	allCM=NULL
	for (i in cacheNames) {
		message(i)
		simpleCache(i, cacheSubDir=cacheSubDir, reload=TRUE, assignToVariable="mat")
		if (is.null(allCM)) {
			allCM = matrix(0, nrow=nrow(mat), ncol=ncol(mat))
		}
		allCM = allCM+mat
	}
	return(allCM)
}

#' Given some sequences (DNAString or character), and a pwm,
#' returns a scaled summarized match scoring matrix.
#' @param sequences DNAStringSet or character set of sequences to search
#' @param pwm Motif to scan.
pwmMatchSignal = function(sequences, pwm) { 
	if (is.character(sequences)) {
		sequences = as.character(sequences);
	}
	pwmScores = lapplyAlias(as.character(sequences), matchPWMscoreI, pwm)
	pwmScoreMatrix = do.call(rbind, pwmScores)
	scale(colSums(pwmScoreMatrix))
}

#' Given a matrix, (and optionally subset ids for nows and cols), just
#' smashes the matrix into a sum, scaled model.
summarizeMatrixModel = function(signalMatrix, subsetRows=NULL, subsetCols=NULL) {
	if (is.null(subsetRows)) {
		subsetRows = 1:nrow(signalMatrix);
	}
	if (is.null(subsetCols)) {
		subsetCols = 1:ncol(signalMatrix);
	}
	mod = apply(signalMatrix[subsetRows,subsetCols], 2, sum)
	scalemod = scale(mod)[,1]
}



# and a function to do the same:
#' Give a name for the analysis; then give a GRanges, a signal Matrix, and
#' genomes. It will extract the DNA sequence for those GRanges, and
#' calculate the tranpose enzyme PWM match; as well as model
#' the signal you give.
seqBias = function(factorName, GR, PWM, BSG=Hsapiens, BSGM=Hsapiens.masked) {
	message("Get sequences...")
	if (length(GR) < 2) { stop("Uhh, give me some regions, man"); }
	simpleCache(paste0(factorName, "_Models"), { 
		seqs = getSeq(BSG, GR)
		modelPWM = pwmMatchSignal(seqs, PWM)
		modelDinuc = getATDinucSignalInStringSet(seqs)
		nlist(modelPWM, modelDinuc);
	}, assignToVariable="models", buildEnvir=nlist(GR, BSG, PWM))

	simpleCache(paste0(factorName, "_ModelsMasked"), { 
		seqsMasked = extractSequences(GR, BSGM, reorder=FALSE) 
		modelPWM = pwmMatchSignal(seqsMasked, PWM)
		modelDinuc = getATDinucSignalInStringSet(seqsMasked)
		nlist(modelPWM, modelDinuc);
	}, assignToVariable="modelsMasked", buildEnvir=nlist(GR, BSGM, PWM))
	
	return(nlist(models, modelsMasked, factorName))
}

seqBiasPlot = function(factorName, models, modelsMasked, signals) {
	maxLim = max(abs(c(models$modelPWM, modelsMasked$modelPWM)))
	plot(gscale(models$modelPWM)-10, smooth.spline(models$modelPWM)$y, type="l", ylim=c(-maxLim, maxLim), lty="dashed",	xlab="Genome")
	lines(gscale(modelsMasked$modelPWM)-10, smooth.spline(modelsMasked$modelPWM)$y, type="l", col="gray")
	#lines(gscale(m[1:380])-10, m[1:380]-models$modelPWM[1:380], type="l", xlab="TSS", col="blue", lty="dashed")
	dn = scale(models$modelDinuc)
	lines(gscale(dn), smooth.spline(dn)$y, type="l", col="red", lty="dashed")

	for (i in 1:length(signals)) {
		lines(gscale(signals[[i]]), signals[[i]], type="l", col=i)
	}
	legend('topleft', c("PWM match (rpts)", "PWM match (masked)", "AT", names(signals)), col=c("blue", "gray", "red", 1:length(signals)), lty=c(2,2,2, rep(1, length(signals))), bty='n', cex=.9)
	crosshair()
	mtext(factorName)
}
#' with(sbo, seqBiasPlot(pwmScoreMatrix, pwmScoreMatrixMasked, m, factorName))









# Some notes on early command-line versions of exact cuts:


#/fhgfs/groups/lab_bock/arendeiro/share/chipmentationBamFiles

#bamToBed -i /fhgfs/groups/lab_bock/arendeiro/share/chipmentationBamFiles/H3K4me3_K562_500k_CM.bam > temp.bed

# bedToExactWig.pl temp.bed ~/fhgfs/share/data/ucsc/chromInfo.hg19.txt out


#/fhgfs/groups/lab_bock/arendeiro/share/chipmentationBedFiles/H3K4me3_K562_500k_ChIP.5prime.bw


#bigWigSummaryOverBed /fhgfs/groups/lab_bock/arendeiro/share/chipmentationBedFiles/H3K4me3_K562_500k_ChIP.5prime.bw promoters.bed out.tab 200





