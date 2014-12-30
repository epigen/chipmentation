####
#
# After this, samples should be renamed to $SAMPLE_NAME.bam
#
####

PROJECTDIR=/home/arendeiro/data/human/chipmentation

for SAMPLE_NAME in H3K27me3_K562_500k_ChIP_CM11-10_R1 H3K27me3_K562_500k_ChIP_CM12-10_R2 H3K27me3_K562_500k_CM_CM11-2_R1 \
H3K27me3_K562_10k_CM_CM11-6_R1 H3K27me3_K562_10k_CM_CM11-6_R2 \
H3K27me3_K562_500k_CM_CM12-2_R1 H3K4me3_K562_500k_ChIP_CM11-9_R1 H3K4me3_K562_500k_ChIP_CM12-9_R1 \
H3K4me3_K562_500k_CM_CM11-1_R1 H3K4me3_K562_500k_CM_CM12-1_R2 IgG_K562_10k_CM_CM11-7_R1 IgG_K562_10k_CM_CM12-7_R2 \
 IgG_K562_500k_ChIP_CM12-11_R2 IgG_K562_500k_CM_CM11-3_R1 IgG_K562_500k_CM_CM12-3_R2 ATAC-seq_ASP14_50k_untreated_ATAC15-4_R2
do
samtools index $SAMPLE_NAME.bam
done

cp $PROJECTDIR/mapped/IgG_K562_500k_ChIP_CM11-11_R1.trimmed.bowtie2.sorted.dup.bam $PROJECTDIR/mapped/IgG_K562_500k_ChIP_CM11-11_R1.bam



for SAMPLE_NAME in H3K27me3_K562_10k_CM_CM11-6_R1 H3K27me3_K562_10k_CM_CM11-6_R2
do
cp $PROJECTDIR/mapped/${SAMPLE_NAME}.bam /fhgfs/groups/lab_bock/shared/data/mapped/
done

