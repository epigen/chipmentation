CONTROL=K562_500K_CM_IGG_nan_nan_0_0_hg19

for SAMPLE in K562_500K_CM_CTCF_nan_nan_0_0_hg19 K562_500K_CM_PU1_nan_nan_0_0_hg19
do
  mkdir -p /projects/chipmentation/data/peaks/$SAMPLE

  macs2 callpeak \
  -t /projects/chipmentation/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam \
  -c /projects/chipmentation/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam \
  --bw 200 \
  -g hs \
  -n $SAMPLE \
  --outdir /projects/chipmentation/data/peaks/${SAMPLE}
done


CONTROL=K562_100K_CM_IGG_nan_nan_1_1_hg19

for SAMPLE in K562_100K_CM_GATA1_nan_nan_0_0_hg19 K562_100K_CM_REST_nan_nan_0_0_hg19
do
  mkdir -p /projects/chipmentation/data/peaks/$SAMPLE

  macs2 callpeak \
  -t /projects/chipmentation/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam \
  -c /projects/chipmentation/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam  \
  --bw 200 \
  -g hs \
  -n $SAMPLE\
  --outdir /projects/chipmentation/data/peaks/${SAMPLE}
done
