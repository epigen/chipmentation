# R Project initialization


# This is how I rsynced the raw data to my local hard drive during the corruption
# rsync -av /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/ /media/nsheffield/red6/chipmentation/data/mapped

# all data:
# rsync -av rsync -av /fhgfs/groups/lab_bock/shared/projects/chipmentation/data /media/nsheffield/red6/chipmentation/

utility("funcGenomeSignals.R")
utility("funcGenomeLocations.R")
loadBSgenome("hg19")



baseDir = getOption("PROJECT.DATA.BASE")

loadPSA = function() {
	#psa = fread("metadata/chipmentation.sample_annotation.csv")
	#psa[, sampleName := paste(cellLine, numberCells, technique, ip, patient, treatment, biologicalReplicate, technicalReplicate, sep="_")]
	psa = fread("metadata/chipmentation.replicates.annotation_sheet.csv")
	psa[, filePath := sub("/fhgfs/groups/lab_bock/shared/projects/", baseDir, filePath)]
	psa[which(!file.exists(filePath)),]
	psa
}


