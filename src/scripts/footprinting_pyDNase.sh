

SAMPLES_FILE=/home/arendeiro/projects/chipmentation/samples_peaks.txt
PROJECTDIR=/home/arendeiro/data/human/chipmentation
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt
GENOMEREF=/fhgfs/prod/ngs_resources/genomes/hg19/forBowtie2/hg19.fa
BEDTOOLSDIR=/home/arendeiro/.local/software/bedtools2/bin
HOMERDIR=/home/arendeiro/.local/software/homer-4.6/bin
CONSERVATION=/home/arendeiro/reference/Homo_sapiens/phyloP/placentalMammals
CAGE=/home/arendeiro/reference/Homo_sapiens/cage_K562_cell_tss_f10.bed

cut -f 1,2,3 $PROJECTDIR/peaks/PU1_K562_10mio_CM_CM15-1.narrowPeak $PROJECTDIR/peaks/PU1_K562_10mio_CM_CM15-1.bed

wellington_footprints.py $PROJECTDIR/peaks/PU1_K562_10mio_CM_CM15-1.bed $PROJECTDIR/mapped/PU1.merged.bam $PROJECTDIR/pyDNASE_PU1

dnase_average_profile.py PU1_peaks.bed PU1_K562_500k_CM_CM11-4.bam PU1_profile


