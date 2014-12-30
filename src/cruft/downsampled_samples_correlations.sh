### Activate virtual environment
source /home/arendeiro/venv/bin/activate

### Specify paths
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19

# count reads
for QUANT in 100 10 2
do
    sambamba flagstat -t 16 $PROJECTDIR/mapped/H3K4me3_500k_PBMC_CM_next_${QUANT}pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.bam \
    > H3K4me3_500k_PBMC_CM_next_${QUANT}pg_ChIP8-8_R1.flagstat

    sambamba flagstat -t 16 $PROJECTDIR/mapped/Input_500k_PBMC_CM_next_${QUANT}pg_ChIP8-Input_R1.trimmed.bowtie2.sorted.shifted.dup.bam \
    > Input_500k_PBMC_CM_next_${QUANT}pg_ChIP8-Input_R1.flagstat
done

sambamba flagstat -t 16 $PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.bam \
> H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.flagstat

sambamba flagstat -t 16 $PROJECTDIR/mapped/Input_500k_PBMC_ChIP_ChIP8-Input_R1.trimmed.bowtie2.sorted.dup.bam \
> Input_500k_PBMC_ChIP_ChIP8-Input_R1.trimmed.bowtie2.sorted.dup.flagstat
#3923340
#11143912
#11099266
#18743527
#8336833
#11247474
#44963030
#46392828

# subsample
cp $PROJECTDIR/mapped/H3K4me3_500k_PBMC_CM_next_10pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.bam $PROJECTDIR/mapped/H3K4me3_500k_PBMC_CM_next_10pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled.bam

sambamba view -f bam -t 16 -s 0.352 --subsampling-seed=2014 \
$PROJECTDIR/mapped/H3K4me3_500k_PBMC_CM_next_10pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.bam \
-o $PROJECTDIR/mapped/H3K4me3_500k_PBMC_CM_next_10pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled.bam

sambamba view -f bam -t 16 -s 0.353 --subsampling-seed=2014 \
$PROJECTDIR/mapped/H3K4me3_500k_PBMC_CM_next_2pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.bam \
-o $PROJECTDIR/mapped/H3K4me3_500k_PBMC_CM_next_2pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled.bam

sambamba view -f bam -t 16 -s 0.209 --subsampling-seed=2014 \
$PROJECTDIR/mapped/Input_500k_PBMC_CM_next_100pg_ChIP8-Input_R1.trimmed.bowtie2.sorted.shifted.dup.bam \
-o $PROJECTDIR/mapped/Input_500k_PBMC_CM_next_100pg_ChIP8-Input_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled.bam

sambamba view -f bam -t 16 -s 0.471 --subsampling-seed=2014 \
$PROJECTDIR/mapped/Input_500k_PBMC_CM_next_10pg_ChIP8-Input_R1.trimmed.bowtie2.sorted.shifted.dup.bam \
-o $PROJECTDIR/mapped/Input_500k_PBMC_CM_next_10pg_ChIP8-Input_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled.bam

sambamba view -f bam -t 16 -s 0.349 --subsampling-seed=2014 \
$PROJECTDIR/mapped/Input_500k_PBMC_CM_next_2pg_ChIP8-Input_R1.trimmed.bowtie2.sorted.shifted.dup.bam \
-o $PROJECTDIR/mapped/Input_500k_PBMC_CM_next_2pg_ChIP8-Input_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled.bam

sambamba view -f bam -t 16 -s 0.087 --subsampling-seed=2014 \
$PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.bam \
-o $PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.subsampled.bam

sambamba view -f bam -t 16 -s 0.085 --subsampling-seed=2014 \
$PROJECTDIR/mapped/Input_500k_PBMC_ChIP_ChIP8-Input_R1.trimmed.bowtie2.sorted.dup.bam \
-o $PROJECTDIR/mapped/Input_500k_PBMC_ChIP_ChIP8-Input_R1.trimmed.bowtie2.sorted.dup.subsampled.bam


# call peaks
source /home/arendeiro/venv/bin/activate
PROJECTDIR=/home/arendeiro/data/human/chipmentation

for QUANT in 100 10 2
do
    macs2 callpeak -t $PROJECTDIR/mapped/H3K4me3_500k_PBMC_CM_next_${QUANT}pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled.bam \
    -c $PROJECTDIR/mapped/Input_500k_PBMC_CM_next_${QUANT}pg_ChIP8-Input_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled.bam \
    --bw 200 \
    -g hs -n H3K4me3_500k_PBMC_CM_next_${QUANT}pg --outdir $PROJECTDIR/peaks/H3K4me3_500k_PBMC_CM_next_${QUANT}pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled_peaks
done

macs2 callpeak -t $PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.subsampled.bam \
-c $PROJECTDIR/mapped/Input_500k_PBMC_ChIP_ChIP8-Input_R1.trimmed.bowtie2.sorted.dup.subsampled.bam \
--bw 200 \
-g hs -n H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1 --outdir $PROJECTDIR/peaks/$PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.subsampled_peaks

deactivate

# count reads on windows genome-wide
for QUANT in 100 10 2
do
    bedtools bamtobed -i $PROJECTDIR/mapped/H3K4me3_500k_PBMC_CM_next_${QUANT}pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled.bam | \
    bedtools intersect -c -a $WINDOWS -b stdin | \
    bedtools sort \
    > $PROJECTDIR/bed/correlations/H3K4me3_500k_PBMC_CM_next_${QUANT}pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled_1kb_windows.bed

    cut -f 4 $PROJECTDIR/bed/correlations/H3K4me3_500k_PBMC_CM_next_${QUANT}pg_ChIP8-8_R1.trimmed.bowtie2.sorted.shifted.dup.subsampled_1kb_windows.bed | \
    sed "1s/^/$SAMPLE_NAME\n/" \
    > $PROJECTDIR/bed/correlations/H3K4me3_500k_PBMC_CM_next_${QUANT}pg_ChIP8-8_R1_.trimmed.bowtie2.sorted.shifted.dup.subsampled.1kb_windows_1col.bed
done

bedtools bamtobed -i $PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.subsampled.bam | \
bedtools intersect -c -a $WINDOWS -b stdin | \
bedtools sort \
> $PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.subsampled_1kb_windows.bed

cut -f 4 $PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.subsampled_1kb_windows.bed | \
sed "1s/^/$SAMPLE_NAME\n/" \
> $PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.subsampled.1kb_windows_1col.bed

# concatenate counts

COUNTS=$PROJECTDIR/downsampled_counts_1kb_windows.tsv
cp $PROJECTDIR/mapped/H3K4me3_500k_PBMC_ChIP_ChIP8-8_R1.trimmed.bowtie2.sorted.dup.subsampled.1kb_windows_1col.bed $COUNTS
paste $COUNTS $PROJECTDIR/bed/correlations/H3K4me3_500k_PBMC_CM_next_100pg_ChIP8-8_R1_.trimmed.bowtie2.sorted.shifted.dup.subsampled.1kb_windows_1col.bed > tmp
mv tmp $COUNTS
paste $COUNTS $PROJECTDIR/bed/correlations/H3K4me3_500k_PBMC_CM_next_10pg_ChIP8-8_R1_.trimmed.bowtie2.sorted.shifted.dup.subsampled.1kb_windows_1col.bed > tmp
mv tmp $COUNTS
paste $COUNTS $PROJECTDIR/bed/correlations/H3K4me3_500k_PBMC_CM_next_2pg_ChIP8-8_R1_.trimmed.bowtie2.sorted.shifted.dup.subsampled.1kb_windows_1col.bed > tmp
mv tmp $COUNTS



R

require(lattice)
library(ggplot2)
library(reshape2)
library(LSD)

projectDir <- "/home/arendeiro/projects/chipmentation/"
dataDir <- "/home/arendeiro/data/human/chipmentation/"

counts <- read.table(paste(dataDir, "downsampled_counts_1kb_windows.tsv", sep = ""), header = TRUE)
names(counts) <- c("ChIP","CM_100pg","CM_10pg","CM_2pg")

#### Scatterplots with all data between replicates

### Raw correlations
rawCounts <- counts

i = 1
for (sample in seq(1, length(counts))) {
    rawCounts[i] <- counts[ , sample] / sum(counts[ , sample])
    i = i + 1
}

#rawCor <- cor(rawCounts[seq(1, length(counts), 2)])
rawCor <- cor(rawCounts)
write.table(rawCor, paste(projectDir, "results/downsampled_counts_1kb_windows.tsv", sep = ""))

#### Plot
rawCor <- melt(rawCor)
rawCor <- rawCor[order(rawCor$Var1),]
p <- qplot(x=Var1, y=Var2, data = rawCor, fill=value, geom="tile")
ggsave(filename = paste(projectDir, "results/plots/downsampled_correlations_10kb_windows_raw.pdf", sep = ""), plot=p, width=15, height=15, units="in")
