while read SAMPLE_NAME CONTROL_NAME; do

SAMPLE_NAME=H3K27me3_K562_10k_CM
CONTROL_NAME=IgG_K562_10k_CM

ss $SAMPLE_NAME

#salloc --partition=mediumq --time=12:00:00 --job-name=$SAMPLE_NAME --nodes=1 --cpus-per-task=8 --mem=120000 sh /home/arendeiro/projects/chipmentation/src/scripts/call_peaks_spp.R $SAMPLE_NAME $CONTROL_NAME

SAMPLE_NAME=H3K27me3_K562_10k_CM
CONTROL_NAME=IgG_K562_10k_CM

module load R
module load openmpi/gcc/64/1.8.2-mlnx-ofed2
module load openmpi/gcc/64/1.8.1

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

PROJECTDIR=/home/arendeiro/data/human/chipmentation

### Call Peaks
mkdir -p $PROJECTDIR/spp_peaks/

Rscript /home/arendeiro/projects/chipmentation/src/scripts/call_peaks_spp.R $SAMPLE_NAME $CONTROL_NAME

date

screen -d 

done < $PEAKS_FILE_MERGED

