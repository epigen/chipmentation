# R Project initialization


# This is how I rsynced the raw data to my local hard drive during the corruption
# rsync -av /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/ /media/nsheffield/red6/chipmentation/data/mapped

# all data:
# rsync -av rsync -av /fhgfs/groups/lab_bock/shared/projects/chipmentation/data /media/nsheffield/red6/chipmentation/

utility("funcGenomeSignals.R")
utility("funcGenomeLocations.R")
utility("funcLoadSharedData.R")
utility("funcEnrichment.R") #for topn
loadBSgenome("hg19")

dataDir = paste0(getOption("PROJECT.DATA.DIR"), "/data/")
resDir = paste0(getOption("PROJECT.DATA.DIR"), "/results/")

dir.create(resDir, showWarnings=FALSE);

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
	tss = fread(paste0(dataDir, "/hg19.cage_peak_coord_robust.TATA_Annotated.bed"), sep="\t")
	tss[, group:=interaction(V10, V11)]
	tss
	#tss400 = 
}

# Paths to communal data files

dat = list()
dat$cage = paste0(dataDir, "hg19.cage_peak_coord_robust.400bp.bed")


