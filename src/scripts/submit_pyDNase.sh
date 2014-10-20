PROJECTDIR=/home/arendeiro/data/human/chipmentation

for SAMPLE_NAME in PU1_K562_10mio_CM
do
    # get simpler bed file
    cut -f 1,2,3 $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks/${SAMPLE_NAME}_peaks.narrowPeak > $PROJECTDIR/peaks/${SAMPLE_NAME}_peaks/${SAMPLE_NAME}_peaks.bed

    # Submit jobs with different parameters
    FS=(8 11 14 17 20 8 11 14 17 20 8 11 14 17 20 8 11 14 17 20 8 11 14 17 20)
    FE=(38 41 44 47 50 38 41 44 47 50 38 41 44 47 50 38 41 44 47 50 38 41 44 47 50)
    SS=(10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10)
    SE=(11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11)

    i=0
    while [ $i -lt ${#FE[*]} ]; do
        echo "Doing sample" $SAMPLE_NAME "with options: ${SS[$i]}_${SE[$i]}_${FE[$i]}"
        sbatch /home/arendeiro/projects/chipmentation/src/scripts/run_pyDNase.sh $SAMPLE_NAME ${SS[$i]} ${SE[$i]} ${FS[$i]} ${FE[$i]}
        i=$(( $i + 1));
    done

done

# first run
#FE=(21 26 31 36 41 21 26 31 36 41 21 26 31 36 41 21 26 31 36 41 21 26 31 36 41)
#SS=(10 10 10 10 10 20 20 20 20 20 35 35 35 35 35 55 55 55 55 55 70 70 70 70 70)
#SE=(11 11 11 11 11 21 21 21 21 21 36 36 36 36 36 56 56 56 56 56 71 71 71 71 71)

