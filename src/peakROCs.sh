# ROC curves of peaks recovered from Encode

PEAKS=/home/afr/cemm-backup/chipmentation/data/peaks
OUTPUT=${PEAKS}/roc.txt

# H3K4ME3
IP=H3K4ME3
CM=${PEAKS}/K562_500K_CM_${IP}_nan_nan_0_0_hg19/K562_500K_CM_${IP}_nan_nan_0_0_hg19_peaks.narrowPeak
CHIP=${PEAKS}/K562_500K_CHIP_${IP}_nan_nan_0_0_hg19/K562_500K_CHIP_${IP}_nan_nan_0_0_hg19_peaks.narrowPeak
ENCODE=${PEAKS}/wgEncodeBroadHistoneK562H3k4me3StdAln/wgEncodeBroadHistoneK562H3k4me3StdAln_peaks.narrowPeak

TRUE=`bedtools intersect -u -a $CM -b $ENCODE | wc -l`
TOTAL=`wc -l $CM`
TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
echo $IP CM $TRUE $TOTAL >> $OUTPUT

TRUE=`bedtools intersect -u -a $CHIP -b $ENCODE | wc -l`
TOTAL=`wc -l $CHIP`
TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
echo $IP CHIP $TRUE $TOTAL >> $OUTPUT


# H3K4ME1
IP=H3K4ME1
CM=${PEAKS}/K562_10M_CM_${IP}_nan_nan_0_0_hg19/K562_10M_CM_${IP}_nan_nan_0_0_hg19_peaks.narrowPeak
CHIP=${PEAKS}/K562_10M_CHIP_${IP}_nan_nan_0_0_hg19/K562_10M_CHIP_${IP}_nan_nan_0_0_hg19_peaks.narrowPeak
ENCODE=${PEAKS}/wgEncodeBroadHistoneK562H3k4me1StdAln/wgEncodeBroadHistoneK562H3k4me1StdAln_peaks.narrowPeak

TRUE=`bedtools intersect -u -a $CM -b $ENCODE | wc -l`
TOTAL=`wc -l $CM`
TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
echo $IP CM $TRUE $TOTAL >> $OUTPUT

TRUE=`bedtools intersect -u -a $CHIP -b $ENCODE | wc -l`
TOTAL=`wc -l $CHIP`
TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
echo $IP CHIP $TRUE $TOTAL >> $OUTPUT


# # H3K27AC
# IP=H3K27AC
# CM=${PEAKS}/K562_10M_CM_${IP}_nan_nan_0_0_hg19/K562_10M_CM_${IP}_nan_nan_0_0_hg19_peaks.broadPeak
# CHIP=${PEAKS}/K562_10M_CHIP_${IP}_nan_nan_1_0_hg19/K562_10M_CHIP_${IP}_nan_nan_1_0_hg19_peaks.broadPeak
# ENCODE=${PEAKS}/wgEncodeBroadHistoneK562H3k36me3StdAln/wgEncodeBroadHistoneK562H3k36me3StdAln_peaks.broadPeak

# TRUE=`bedtools intersect -u -a $CM -b $ENCODE | wc -l`
# TOTAL=`wc -l $CM`
# TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
# echo $IP CM $TRUE $TOTAL >> $OUTPUT

# TRUE=`bedtools intersect -u -a $CHIP -b $ENCODE | wc -l`
# TOTAL=`wc -l $CHIP`
# TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
# echo $IP CHIP $TRUE $TOTAL >> $OUTPUT


# H3K27ME3
IP=H3K27ME3
CM=${PEAKS}/K562_500K_CM_${IP}_nan_nan_0_0_hg19/K562_500K_CM_${IP}_nan_nan_0_0_hg19_peaks.broadPeak
CHIP=${PEAKS}/K562_500K_CHIP_${IP}_nan_nan_0_0_hg19/K562_500K_CHIP_${IP}_nan_nan_0_0_hg19_peaks.broadPeak
ENCODE=${PEAKS}/wgEncodeBroadHistoneK562H3k27me3StdAln/wgEncodeBroadHistoneK562H3k27me3StdAln_peaks.broadPeak

TRUE=`bedtools intersect -u -a $CM -b $ENCODE | wc -l`
TOTAL=`wc -l $CM`
TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
echo $IP CM $TRUE $TOTAL >> $OUTPUT

TRUE=`bedtools intersect -u -a $CHIP -b $ENCODE | wc -l`
TOTAL=`wc -l $CHIP`
TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
echo $IP CHIP $TRUE $TOTAL >> $OUTPUT



# H3K36ME3
IP=H3K36ME3
CM=${PEAKS}/K562_10M_CM_${IP}_nan_nan_0_0_hg19/K562_10M_CM_${IP}_nan_nan_0_0_hg19_peaks.broadPeak
CHIP=${PEAKS}/K562_10M_CHIP_${IP}_nan_nan_1_1_hg19/K562_10M_CHIP_${IP}_nan_nan_1_1_hg19_peaks.broadPeak
ENCODE=${PEAKS}/wgEncodeBroadHistoneK562H3k36me3StdAln/wgEncodeBroadHistoneK562H3k36me3StdAln_peaks.broadPeak

TRUE=`bedtools intersect -u -a $CM -b $ENCODE | wc -l`
TOTAL=`wc -l $CM`
TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
echo $IP CM $TRUE $TOTAL >> $OUTPUT

TRUE=`bedtools intersect -u -a $CHIP -b $ENCODE | wc -l`
TOTAL=`wc -l $CHIP`
TOTAL=`echo $TOTAL | grep -o "^[0-9]* "`
echo $IP CHIP $TRUE $TOTAL >> $OUTPUT

# scikitlearn
# make roc curves
