#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=6:00:00

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000
#SBATCH --nodes=1

#SBATCH --job-name=wget
#SBATCH --output=/home/arendeiro/logs/wget_%j.out

hostname
date

URL=$1
OUTPUT=$2

echo $URL

wget -O $OUTPUT $URL

date
