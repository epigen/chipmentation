#!/bin/bash

DIRS=(
    /fhgfs/groups/lab_bsf/samples/BSF_0145_C5PL4ACXX/BSF_0145_C5PL4ACXX_1_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0081_C3N4DACXX/BSF_0081_C3N4DACXX_2_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0083_C3MNWACXX/BSF_0083_C3MNWACXX_2_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0086_C3NY0ACXX/BSF_0086_C3NY0ACXX_2_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0113_C558FACXX/BSF_0113_C558FACXX_3_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0111_C568EACXX/BSF_0111_C568EACXX_4_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0127_C5CJ1ACXX/BSF_0127_C5CJ1ACXX_6_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0139_C64UYACXX/BSF_0139_C64UYACXX_1_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0120_C5CVMACXX/BSF_0120_C5CVMACXX_3_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0146_C5NUPACXX/BSF_0146_C5NUPACXX_7_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0127_C5CJ1ACXX/BSF_0127_C5CJ1ACXX_5_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0142_C55ALACXX/BSF_0142_C55ALACXX_1_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0120_C5CVMACXX/BSF_0120_C5CVMACXX_1_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0107_C4GLTACXX/BSF_0107_C4GLTACXX_5_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0120_C5CVMACXX/BSF_0120_C5CVMACXX_2_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0109_C4W36ACXX/BSF_0109_C4W36ACXX_3_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0129_C5HM7ACXX/BSF_0129_C5HM7ACXX_8_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0146_C5NUPACXX/BSF_0146_C5NUPACXX_8_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0088_C3U90ACXX/BSF_0088_C3U90ACXX_2_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0110_C561TACXX/BSF_0110_C561TACXX_7_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0090_C3TYGACXX/BSF_0090_C3TYGACXX_3_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0102_C49ATACXX/BSF_0102_C49ATACXX_2_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0082_C3JPBACXX/BSF_0082_C3JPBACXX_2_samples/
    /fhgfs/groups/lab_bsf/samples/BSF_0051_D2GW1ACXX/BSF_0051_D2GW1ACXX_4_samples/
)

# get corrupted bams from bsf folder
for DIR in ${DIRS[*]}
do
    for FILE in `ls $DIR | grep .bam$`
    do
        cat ${DIR}/${FILE}.md5 >> ~/md5.txt
        echo " " ${DIR}/$FILE >> ~/md5.txt
    done
    md5sum -c ~/md5.txt >> ~/md5.log
    rm ~/md5.txt 
done
grep FAILED md5.log


# get corrupted bams from own folder
while read LINE
do
    sbatch ~/checkIntegrity.sh $LINE
done < ~/allBams.txt

ll | grep K562.*\.log | grep 82
