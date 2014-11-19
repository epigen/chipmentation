#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --ntasks=8
#SBATCH --time=12:00:00

# Optional parameters
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10000
#SBATCH --nodes=1
#SBATCH --job-name=spp_peaks
#SBATCH --output=/home/arendeiro/logs/spp_peaks_%j.out

# *** setup environment ***
# load the required environmental modules that your script depends upon
module load R
module load openmpi/gcc/64/1.8.2-mlnx-ofed2
module load openmpi/gcc/64/1.8.1

source /home/arendeiro/.bash_profile
source /home/arendeiro/.bashrc

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

PROJECTDIR=/home/arendeiro/data/human/chipmentation

### Get sample info from arguments
SAMPLE_NAME=$1
CONTROL_NAME=$2

### Call Peaks
mkdir -p $PROJECTDIR/spp_peaks/

Rscript /home/arendeiro/projects/chipmentation/src/scripts/call_peaks_spp.R $SAMPLE_NAME $CONTROL_NAME

date
