#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Optional parameters
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8000
#SBATCH --nodes=1
#SBATCH --job-name=call_peaks
#SBATCH --output=/home/arendeiro/logs/call_peaks_%j.out

# *** setup environment ***

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments

### Activate virtual environment
source /home/arendeiro/venv/bin/activate

PROJECTDIR=/home/arendeiro/data/human/chipmentation

SAMPLE_NAME=H3K4me3_K562_500k_ChIP_Encode
CONTROL_NAME=Input_K562_500k_ChIP_Encode


macs2 callpeak -t $PROJECTDIR/mapped/merged/$SAMPLE_NAME.bam \
-c $PROJECTDIR/mapped/merged/$CONTROL_NAME.bam \
--bw 200 \
-g hs -n ${SAMPLE_NAME} --outdir $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks_bw200

deactivate
date