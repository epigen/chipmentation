


ROBUST=/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.bed

sed '/^#/ d' /home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.osc \
| tail -n +2 > /home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.bed

TATANOT=reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.bed

bedtools intersect -wa \
-a $ROBUST \
-b $TATANOT \
| head

robust = read.table("/fhgfs/groups/lab_bock/shared/data/cage_tss/hg19.cage_peak_coord_robust.bed", sep = "\t", header = FALSE)
annotation = read.table("/home/arendeiro/reference/Homo_sapiens/DPIcluster_hg19_20120116.permissive_set.TATA_CpG_annotated.bed", sep = "\t", header = FALSE)
colnames(robust)[4] = "peak"
colnames(annotation)[4] = "peak"

annotation = annotation[, c()]

