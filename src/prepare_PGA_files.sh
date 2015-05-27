mkdir ~/PGA1Nextera
cd ~/PGA1Nextera

# Check files are not corrupt

for FILE in /fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_3_PGA_0001_GD_Nextera_recalibrated.bam /fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_3_PGA_0001_TD_Nextera_recalibrated.bam /fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_4_PGA_0001_GD_Nextera_recalibrated.bam /fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_4_PGA_0001_TD_Nextera_recalibrated.bam
do
sambamba flagstat -t 12 $FILE >> flagstat.log
done

# This was corrupt, removed from list
# /fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_4_PGA_0001_GD_Nextera_recalibrated.bam


# Merge
# sambamba merge -p -t 16 PGA_0001_Nextera.bam /fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_3_PGA_0001_GD_Nextera_recalibrated.bam /fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_3_PGA_0001_TD_Nextera_recalibrated.bam /fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_4_PGA_0001_TD_Nextera_recalibrated.bam

java -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/MergeSamFiles.jar \
USE_THREADING=TRUE \
INPUT=/fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_3_PGA_0001_GD_Nextera_recalibrated.bam \
INPUT=/fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_3_PGA_0001_TD_Nextera_recalibrated.bam \
INPUT=/fhgfs/groups/lab_bsf/projects/BSA_0010_PGA/b37/variant_calling_process_lane_SET_53_C248CACXX_4_PGA_0001_TD_Nextera_recalibrated.bam \
OUTPUT=PGA_0001_Nextera-2.bam

# Index
sambamba index -t 16 PGA_0001_Nextera-2.bam
