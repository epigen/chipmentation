project.init2("chipmentation")

# Produce a bed file with 400bp surrounding each cage peak,
# to extract exact cuts.

tss =fread(paste0(baseDir, "chipmentation/data/hg19.cage_peak_coord_robust.TATA_Annotated.bed"), sep="\t")

tss[, V2:=pmax(0, V2-200)]
tss[, V3:=V3+200]
tss[, uniqueName:=1:nrow(bed)]
tss
tss[, c(1,2,3, 12), with=FALSE]

write.tsv(tss[, c(1,2,3, 12), with=FALSE], file=paste0(baseDir, "chipmentation/data/hg19.cage_peak_coord_robust.400bp.bed"), col.names=FALSE)

