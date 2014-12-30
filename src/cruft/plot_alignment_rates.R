align = read.table("~/projects/chipmentation/results/alignment rates.csv", sep = ",", header = TRUE)
#align = read.table("~/workspace/ChIPmentation/alignment rates.csv", sep = ",", header = TRUE)
align$name <- paste(align$technique, align$antibody, sep = '_')

library(ggplot2)

p <- ggplot(align, aes(aligned, name)) +
  geom_point(aes(colour = factor(technique))) +
  coord_cartesian(xlim = c(85, 100)) +
  facet_grid(algorithm ~ variation) +
  theme_bw()

ggsave(filename = "~/projects/chipmentation/results/plots/alignment_test_percentAligned.png", plot = p, height = 3, width = 12)

p <- ggplot(align, aes(unique_aligned, name)) +
  geom_point(aes(colour = factor(technique))) +
  coord_cartesian(xlim = c(30, 100)) +
  facet_grid(algorithm ~ variation) +
  theme_bw()
ggsave(filename = "~/projects/chipmentation/results/plots/alignment_test_percentUniquelyAligned.png", plot = p, height = 3, width = 12)
