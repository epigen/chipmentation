### CAGE
# get data
cd /home/arendeiro/reference/Homo_sapiens/
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeRikenCageK562CellPapTssHmm.txt.gz
gzip -d wgEncodeRikenCageK562CellPapTssHmm.txt.gz
# filter on ~score, make bed
awk -v OFS='\t' '{(if $8 > 10) print $2, $3, $4, $5, $6, $7}' /home/arendeiro/reference/Homo_sapiens/wgEncodeRikenCageK562CellPapTssHmm.txt > /home/arendeiro/reference/Homo_sapiens/cage_K562_cell_tss_f10.bed
