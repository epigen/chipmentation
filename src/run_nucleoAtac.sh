# Get nucleossomes
# 1K4
ROOT=/home/arendeiro/chipmentation
SAMPLE=K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19
PEAKS=${ROOT}/data/peaks/${SAMPLE}/${SAMPLE}_peaks.narrowPeak
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/home/arendeiro/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/$SAMPLE
mkdir -p $OUTPUT

source .bashrc
source activate env

samtools index $BAM
nucleoatac run --bed $PEAKS --bam $BAM --fasta $GENOME --out $OUTPUT --cores 16

# 3K4
ROOT=/home/arendeiro/chipmentation
SAMPLE=K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19
PEAKS2=${ROOT}/data/peaks/${SAMPLE}/${SAMPLE}_peaks.narrowPeak
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/home/arendeiro/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/$SAMPLE
mkdir -p $OUTPUT

source .bashrc
source activate env

samtools index $BAM
nucleoatac run --bed $PEAKS --bam $BAM --fasta $GENOME --out $OUTPUT --cores 32

# 1K4 + 3K4
ROOT=/home/arendeiro/chipmentation
SAMPLE=K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19-K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19
PEAKS3=${ROOT}/data/peaks/${SAMPLE}/${SAMPLE}_peaks.narrowPeak
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/home/arendeiro/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/$SAMPLE
mkdir -p $OUTPUT

source .bashrc
source activate env

samtools index $BAM
nucleoatac run --bed $PEAKS --bam $BAM --fasta $GENOME --out $OUTPUT --cores 16


# VPLOTS
# vplot 1k4 on predicted nucleosomes
ROOT=/media/afr/cemm-backup/chipmentation
SAMPLE=K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/${SAMPLE}.nucleosomes.vplot
mkdir -p $OUTPUT

cd $ROOT/data/nucleoATAC/10M_CM_H3K4ME1_PE

sort -k6nr K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.nucpos.bed | cut -f 1,2,3 | bedtools slop -b 250 -i - -g ../../../../reference/hg19/hg19.chrom-sizes.tsv | sortBed > $OUTPUT/$SAMPLE.nucleosomes.bed

pyatac vplot --cores 4 --lower 0 --upper 250 --flank 200 --plot_extra --bed $OUTPUT/$SAMPLE.nucleosomes.bed --bam $BAM --out $OUTPUT/$SAMPLE

# pu1 motif centered
ROOT=/media/afr/cemm-backup/chipmentation
SAMPLE=K562_10M_CM_PU1_nan_PE_1_1_hg19
PEAKS=${ROOT}/data/peaks/K562_10M_CM_PU1_nan_nan_0_0_hg19/K562_10M_CM_PU1_nan_nan_0_0_hg19_peaks.motifCentered.bed
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/$SAMPLE.motifCentered
mkdir -p $OUTPUT

pyatac vplot --cores 4 --lower 0 --upper 250 --flank 200 --plot_extra --bed $PEAKS --bam $BAM --out $OUTPUT



##### ON THE CLUSTER #######
# vplot 1k4 on predicted nucleosomes
ROOT=~/chipmentation
BAM=K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.trimmed.bowtie2.shifted.dups.bam
OUTPUT=10M_CM_H3K4ME1_PE.vplot

pyatac vplot --cores 32 --lower 0 --upper 250 --flank 200 --plot_extra --bed nucleosomes.bed --bam $BAM --out $OUTPUT


ROOT=/media/afr/cemm-backup/chipmentation
SAMPLE=K562_50K_ATAC_nan_nan_nan_0_0_hg19
PEAKS=${ROOT}/data/peaks/${SAMPLE}/${SAMPLE}_peaks.narrowPeak
BAM=${ROOT}/data/mapped/${SAMPLE}.trimmed.bowtie2.shifted.dups.bam
GENOME=/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa
OUTPUT=${ROOT}/data/nucleoATAC/$SAMPLE
mkdir -p $OUTPUT

nucleoatac run --bed $PEAKS --bam $BAM --fasta $GENOME --out $OUTPUT --cores 16





sort -k6nr K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.nucpos.bed | head -n 5000 | cut -f 1,2,3 | sortBed > K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.nucpos.top.bed







# Calculate distance between features (dyads and nucleosomes)
# dyads
awk '{ if ($4 >= 5) print }' K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.nucpos.bed > nucs.filtered.bed


bedtools closest -d -N -t all \
-a K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.nucpos.bed \
-b K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.nucpos.bed \
> dyad_distances.bed


cut -f 25 dyad_distances.bed > dyad_distances.counts.txt

# Calculate distance between features (dyads and nucleosomes)
bedtools slop -b 74 -i K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19.nucpos.bed -g ../../../../reference/hg19/hg19.chrom-sizes.tsv | sortBed > nucleosomes.bed

bedtools closest -d -N -t all \
-a nucleosomes.bed \
-b nucleosomes.bed \
> nucleosomes_distances.bed

cut -f 25 nucleosomes_distances.bed > nucleosomes_distances.counts.txt




import pandas as pd 
from collections import Counter

df = pd.read_csv("dyad_distances.counts.txt", header=None)
df2 = pd.read_csv("nucleosomes_distances.counts.txt", header=None)

c = Counter(df[0])
c2 = Counter(df2[0])

plt.plot(c.keys(), c.values())

plt.savefig("")


plt.plot(c2.keys(), c2.values())
