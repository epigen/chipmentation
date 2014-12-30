#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --time=48:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32000
#SBATCH --nodes=1

#SBATCH --job-name=cageParsePeaksCoverage
#SBATCH --output=/home/arendeiro/logs/cageParsePeaksCoverage_%j.out

# *** setup environment ***
# load the required environmental modules that your script depends upon
module load python/2.7.6

# set the temporary folder to your scratch location (avoid using node local /tmp)
export SCRATCH_PATH=/fhgfs/scratch/users/arendeiro
export TMPDIR=$SCRATCH_PATH/tmp

# *** run the job ***
hostname
date

### Get sample info from arguments
CAGE=$1
TECH=$2
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss/

### Specify paths

### Start work
source /home/arendeiro/venv/bin/activate

echo $TECH
echo "Running python parseBedCoverage.py $TECH"
python /home/arendeiro/projects/chipmentation/src/lib/parseBedCoverage-peaks-strand.py $CAGEDIR/$CAGE.120bpSlop.${TECH}coverage.bed

deactivate

date

