project.init2("chipmentation")
getwd()

# This is how I rsynced the raw data to my local hard drive during the corruption
# rsync -av /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/ /media/nsheffield/red6/chipmentation/data/mapped

baseDir = getOption("PROJECT.DATA.BASE")

loadPSA = function() {
	#psa = fread("metadata/chipmentation.sample_annotation.csv")
	#psa[, sampleName := paste(cellLine, numberCells, technique, ip, patient, treatment, biologicalReplicate, technicalReplicate, sep="_")]
	psa = fread("metadata/chipmentation.replicates.annotation_sheet.csv")
	psa[, filePath := sub("/fhgfs/groups/lab_bock/shared/projects/", baseDir, filePath)]
	psa[which(!file.exists(filePath)),]
	psa
}


psa = loadPSA()

# Look for combined samples
psa[biologicalReplicate==0 & technicalReplicate==0,]


utility("funcGenomeSignals.R")
utility("funcGenomeLocations.R")


bed =fread(paste0(baseDir, "chipmentation/data/hg19.cage_peak_coord_robust.TATA_Annotated.bed"), sep="\t")

bed
bedGR = dtToGr(bed, "V1", "V2", "V3")
bedGR
seqLength = 50
bedGR= promoters(dtToGr(bed, "V1", "V2", "V3", strand="V6"), upstream=seqLength, downstream=0)


loadBSgenome("hg19")

bam = bamSlurp(psa[biologicalReplicate==0 & technicalReplicate==0 & ip=="",filePath], bedGR)



