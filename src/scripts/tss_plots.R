##### Calculate coverage for each annotation set
CAGEDIR=/fhgfs/groups/lab_bock/shared/data/cage_tss
GENOMESIZE=/fhgfs/prod/ngs_resources/genomes/hg19/hg19_chromLengths_sorted.txt

declare -a arr=("hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed")

# Get window around CAGE peak
for CAGE in "${arr[@]}"
do
    bedtools slop -b 1000 \
    -i $CAGEDIR/$CAGE.bed \
    -g $GENOMESIZE \
    > $CAGEDIR/$CAGE.1000bpSlop.bed
done


# Calculate coverage
## produce bed files from the 3 techniques
# bedtools bamtobed -i /home/arendeiro/data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam \
# | python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
# > /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_CM.bed

# bedtools bamtobed -i /home/arendeiro/data/human/chipmentation/mapped/merged/IgG_K562_500k_CM.bam \
# | python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
# > /home/arendeiro/data/human/chipmentation/bed/IgG_K562_500k_CM.bed

# bedtools bamtobed -i /home/arendeiro/data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_ChIP.bam \
# | python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
# > /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_ChIP.bed

# bedtools bamtobed -i /home/arendeiro/data/human/chipmentation/mapped/merged/DNase_UWashington_K562_mergedReplicates.bam \
# | python /home/arendeiro/projects/chipmentation/src/lib/get5primePosition.py \
# > /home/arendeiro/data/human/chipmentation/bed/DNase_UWashington_K562_mergedReplicates.bed

# calculate coverage
for CAGE in "${arr[@]}"
do
    for TECH in CM IgG ChIP DNase 
    do
        sbatch /home/arendeiro/projects/chipmentation/src/scripts/cage_tss_1000bpcoverage_job.sh $CAGE $TECH
        # bedtools coverage -d -a \
        # /home/arendeiro/data/human/chipmentation/bed/H3K4me3_K562_500k_CM.bed \
        # -b $CAGEDIR/$CAGE.120bpSlop.bed \
        # > $CAGEDIR/$CAGE.120bpSlop.CMcoverage.bed
    done
done

# parse coverage to a csv (also orients TSSs!)
for CAGE in "${arr[@]}"
do
    for TECH in IgG #CM ChIP DNase 
    do
    	echo $CAGEDIR/$CAGE.1000bpSlop.${TECH}coverage.bed
        sbatch /home/arendeiro/projects/chipmentation/src/scripts/cage_tss_coverage-pythonParse_job.sh \
        $CAGEDIR/$CAGE.1000bpSlop.${TECH}coverage.bed
        # python /home/arendeiro/projects/chipmentation/src/lib/parseBedCoverage.py $CAGEDIR/$CAGE.120bpSlop.CMcoverage.bed
    done
done





R
cageDir = "/fhgfs/groups/lab_bock/shared/data/cage_tss/"
cage = "hg19.cage_peak_coord_robust.TATA_Annotated.K562_expressed"

cm = read.csv(paste(cageDir, cage, ".1000bpSlop.CMcoverage.csv", sep = ''))
rownames(cm) = cm[,1]
cm = cm[,-1]
colnames(cm) = seq(-1000, 1000)

igg = read.csv(paste(cageDir, cage, ".1000bpSlop.IgGcoverage.csv", sep = ''))
rownames(igg) = igg[,1]
igg = igg[,-1]
colnames(igg) = seq(-1000, 1000)

chip = read.csv(paste(cageDir, cage, ".1000bpSlop.ChIPcoverage.csv", sep = ''))
rownames(chip) = chip[,1]
chip = chip[,-1]
colnames(chip) = seq(-1000, 1000)

dnase = read.csv(paste(cageDir, cage, ".1000bpSlop.DNasecoverage.csv", sep = ''))
rownames(dnase) = dnase[,1]
dnase = dnase[,-1]
colnames(dnase) = seq(-1000, 1000)


ChIPmentation = colMeans(cm) / (64054226 / 64054226.)
ChIP = colMeans(chip) / (76894657 / 64054226.)
IgG = colMeans(igg) / (66863728 / 64054226.)
DNase = colMeans(dnase) / (66663835 / 64054226.)

colors = c("#3FAB35", "#4169E1", "#696969", "#8B0000")

### Plot average profiles
library(ggplot2)
library("reshape2")
df = cbind(ChIPmentation, ChIP, IgG, DNase)
df = melt(df)

p = ggplot(df, aes(Var1, value, colour = Var2)) +
	#geom_line() +
	stat_smooth(method = "gam", formula = y ~ s(x, k = 240), se = FALSE) + 
	facet_wrap( ~ Var2) + # , ncol = 2, scales = "free"
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank()) #+
	#scale_size(range=c(0.0, 1.2), guide=FALSE)

ggsave(filename = "tss_signal.pdf", plot = p, height = 2, width = 7)
ggsave(filename = "tss_signal_wide.pdf", plot = p, height = 2, width = 4)

p = ggplot(df, aes(Var1, value, colour = Var2)) +
	#geom_line() +
	stat_smooth(method = "gam", formula = y ~ s(x, k = 240), se = FALSE) + 
	#facet_wrap( ~ Var2) + # , ncol = 2, scales = "free"
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank()) #+
	#scale_size(range=c(0.0, 1.2), guide=FALSE)

ggsave(filename = "tss_signal_allinone_wide.pdf", plot = p, height = 2, width = 7)
ggsave(filename = "tss_signal_allinone.pdf", plot = p, height = 2, width = 4)

window = -1000:1000
a = df[df$Var1 %in% window, ]
p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 240), se = FALSE) + 
	facet_grid(Var2 ~ ., scales = "free") + 
	#coord_cartesian(xlim = c(-1000, 1000)) +
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank()) 

ggsave(filename = "tss_signal_1kb.pdf", plot = p, height = 3, width = 6)


window = -60:60
a = df[df$Var1 %in% window, ]

library(grid)
p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5) + 
	#facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	geom_segment(x = -28, xend = -28, y = 0.4, yend = 0.61, colour = "black", size = 0.8, arrow = arrow(length=unit(0.3,"cm"))) + 
	geom_segment(x = 0, xend = 0, y = 0.4, yend = 0.61, colour = "black", size = 0.8, arrow = arrow(length=unit(0.3,"cm"))) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_highlight_arrows.pdf", plot = p, height = 4, width = 4)
ggsave(filename = "tss_signal_120bp_wide_highlight_arrows.pdf", plot = p, height = 2.5, width = 6)


p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5) + 
	facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	#geom_segment(x = -35, xend = -35, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -20, xend = -20, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_facet_highlight.pdf", plot = p, height = 4, width = 4)
ggsave(filename = "tss_signal_120bp_facet_wide_highlight.pdf", plot = p, height = 2.5, width = 6)




# INDEPENDENT


df = cbind(ChIPmentation)
df = melt(df)
window = -60:60
a = df[df$Var1 %in% window, ]

p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5, colour = colors[1]) + 
	#facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	geom_segment(x = -35, xend = -35, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	geom_segment(x = -20, xend = -20, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_wide_CM.pdf", plot = p, height = 2, width = 6)


df = cbind(ChIP)
df = melt(df)
window = -60:60
a = df[df$Var1 %in% window, ]

p = ggplot(a, aes(Var1, value, colour = Var2)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5, colour = colors[2]) + 
	#facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	#geom_segment(x = -35, xend = -35, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -20, xend = -20, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_wide_ChIP.pdf", plot = p, height = 2, width = 6)



df = cbind(DNase)
df = melt(df)
window = -60:60
a = df[df$Var1 %in% window, ]

p = ggplot(a, aes(Var1, value)) +
	#geom_line() + 
	stat_smooth(method = "gam", formula = y ~ s(x, k = 40), se = FALSE, size = 1.5, colour = colors[4]) + 
	#facet_grid(Var2 ~ ., scales = "free") + 
	scale_x_continuous(limits=c(-60, 60)) +
	#coord_cartesian(xlim = c(-200, 200)) +
	#geom_vline(xintercept = -40) + 
	#geom_vline(xintercept = -20) + 
	#geom_segment(x = -35, xend = -35, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -20, xend = -20, y = 0, yend = 1.25, colour = "grey", size = 0.8) + 
	#geom_segment(x = -5, xend = -5, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	#geom_segment(x = 10, xend = 10, y = 0, yend = 1.25, colour = "grey", size = 0.2) + 
	xlab("Distance to peak") +
	ylab("Tags") +
	scale_color_manual(values = colors) +
	theme_bw() +
	theme(legend.title=element_blank())

ggsave(filename = "tss_signal_120bp_wide_DNase.pdf", plot = p, height = 2, width = 6)


### Plot heatmap of CM
library(gplots)
require(made4)

pdf("tss_heatmap_CM.pdf")
heatplot(cm, dend = 'row', labRow = NA, , labCol = NA)
dev.off()

pdf("tss_heatmap_ChIP.pdf")
heatplot(chip, dend = 'row', labRow = NA, , labCol = NA)
dev.off()


window = seq(1000,3000)

pdf("tss_heatmap_1kb_CM.pdf")
heatplot(cm[,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()

pdf("tss_heatmap_1kb_ChIP.pdf")
heatplot(chip[,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()


window = seq(1600,2400)

pdf("tss_heatmap_400bp_CM.pdf")
heatplot(cm[,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()

pdf("tss_heatmap_400bp_ChIP.pdf")
heatplot(chip[,window], dend = 'row', labRow = NA, , labCol = NA)
dev.off()

