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
dat$cage = paste0(dataDir, "hg19.cage_peak_coord_robust.400bp.bed")


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
	return(nlist(psa, msa))
}



loadCageTSS = function() {
	tss = fread(paste0(annoDir, "/hg19.cage_peak_coord_robust.TATA_Annotated.bed"), sep="\t")
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
	tss
	tss[,.N, by=group]
	return(nlist(tss))
	#distance distribution:
	summary(diff(tss$V3))
}



loadTransposePWM = function() {
	downloadCache("transposasePWM", "https://raw.githubusercontent.com/GreenleafLab/NucleoATAC/master/pyatac/pwm/Human2.PWM.txt")
	rownames(tpwm) = DNA_BASES
	return(list(tpwm = as.matrix(transposasePWM)))
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

write.tsv(tss[, c(1,2,3, 12), with=FALSE], file=paste0(data, "/hg19.cage_peak_coord_robust.400bp.bed"), col.names=FALSE)

}










# Some notes on early command-line versions of exact cuts:


#/fhgfs/groups/lab_bock/arendeiro/share/chipmentationBamFiles

#bamToBed -i /fhgfs/groups/lab_bock/arendeiro/share/chipmentationBamFiles/H3K4me3_K562_500k_CM.bam > temp.bed

# bedToExactWig.pl temp.bed ~/fhgfs/share/data/ucsc/chromInfo.hg19.txt out


#/fhgfs/groups/lab_bock/arendeiro/share/chipmentationBedFiles/H3K4me3_K562_500k_ChIP.5prime.bw


#bigWigSummaryOverBed /fhgfs/groups/lab_bock/arendeiro/share/chipmentationBedFiles/H3K4me3_K562_500k_ChIP.5prime.bw promoters.bed out.tab 200





