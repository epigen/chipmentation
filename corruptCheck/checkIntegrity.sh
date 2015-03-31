#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=5:00:00

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#SBATCH --nodes=1

#SBATCH --job-name=checkIntegrity
#SBATCH --output=/home/arendeiro/checkIntegrity.log

# Start running the job
hostname
date

FILE=$1
LOG=`basename "${FILE%.*}.log"`

sambamba index -t 16 $FILE >> ~/$LOG 2>&1

date
