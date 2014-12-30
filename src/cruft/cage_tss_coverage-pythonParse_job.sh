#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
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
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss/

### Specify paths

### Start work
source /home/arendeiro/venv/bin/activate

echo $CAGE
echo "Running python /home/arendeiro/projects/chipmentation/src/lib/parseBedCoverage.py $CAGE"
python /home/arendeiro/projects/chipmentation/src/lib/parseBedCoverage-peaks-strand.py $CAGE

deactivate

date

