# 12.01.2024
# This script generates general statistics and figures on distribution of DE genes (robust only) across the genome of two species. 

if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("circlize", "stringr", "RColorBrewer", "gridGraphics", "gridExtra", "cowplot", "ggplot2")

setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure/Supp_DE_density/")
setwd('/Users/clara/Documents/GitHub.nosync/Tillandsia-compgenomics/Supplementary_figures/II. DE gene distribution/')

DE_Tlei <- read.table("DE_genes.mapped_to_Tlei.PER-WINDOW.density.txt", sep = "\t")
DE_Tfas <- read.table("DE_genes.mapped_to_Tfas.PER-WINDOW.density.txt", sep = "\t")

mean(DE_Tfas$DE_counts) # 1.47 DE genes per kb window in Tfas on average
mean(DE_Tfas$proportion) # 3 % of genes per kb window is DE in Tfas on average

mean(DE_Tlei$DE_counts) # 0.93 DE genes per kb window in Tlei on average
mean(DE_Tlei$proportion) # 3.3 % of genes per kb window are DE in Tlei on average
chrom_sizes <- read.table("chromosomes_Tfas_Tlei.coordinates-circle.txt")
chrom_sizes_Tfas <- chrom_sizes[c(1:25),c(1:2)]
chrom_sizes_Tlei <- chrom_sizes[c(26:51),c(1:2)]
names(chrom_sizes_Tfas)[1] <- 'chrom'
names(chrom_sizes_Tlei)[1] <- 'chrom'
total_gene_counts_Tfas <- aggregate(DE_counts ~ chrom, data = DE_Tfas, sum)
total_gene_counts_Tlei <- aggregate(DE_counts ~ chrom, data = DE_Tlei, sum)

merged_Tfas <- merge(total_gene_counts_Tfas, chrom_sizes_Tfas, by = 'chrom')
merged_Tlei <- merge(total_gene_counts_Tlei, chrom_sizes_Tlei, by = 'chrom')

shapiro.test(merged_Tfas$DE_counts) # => p-value = 0.8917
shapiro.test(merged_Tlei$DE_counts) # => p-value = 0.01748
res <- cor.test(merged_Tfas$DE_counts, merged_Tfas$size) # parametric test 0.007065
res <- cor.test(merged_Tfas$DE_counts, merged_Tfas$size,  method="kendall") # non-parametric 0.01084
res
res <- cor.test(merged_Tlei$DE_counts, merged_Tlei$size,  method="kendall") # non-parametric 0.01084
res # 0.001349

# Create a new column for colors based on 'chrom' values
merged_Tfas$color_group <- 'default' # Default color
merged_Tfas$color_group[merged_Tfas$chrom %in% c('Tfas_chr3', 'Tfas_chr13')] <- 'goldenrod'
merged_Tfas$color_group[merged_Tfas$chrom %in% c('Tfas_chr10', 'Tfas_chr24')] <- 'darkgreen'

merged_Tlei$color_group <- 'default' # Default color
merged_Tlei$color_group[merged_Tlei$chrom %in% c('Tlei_chr3', 'Tlei_chr19')] <- 'goldenrod'
merged_Tlei$color_group[merged_Tlei$chrom %in% c('Tlei_chr13', 'Tlei_chr23')] <- 'darkgreen'
merged_Tlei$color_group[merged_Tlei$chrom %in% c('Tlei_chr14')] <- 'blue'

merged_Tfas$label_group <- NA # Default to NA (which will be ignored in geom_text)
rows_to_label <- c(19, 5, 2, 17)
label_chroms <- c('scaffold 3', 'scaffold 13', 'scaffold 10', 'scaffold 24')
merged_Tfas$label_group[rows_to_label] <- label_chroms

p1 <- ggplot(merged_Tfas, aes(x = (size/1000000), y = DE_counts)) +
  geom_point(aes(color = color_group), size = 2) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.5) +
  labs(x = "Scaffold size (Mb)", y = "DE gene count", title = "T. fasciculata") +
  scale_color_manual(values = c(goldenrod = 'goldenrod', darkgreen = 'forestgreen', default = 'black')) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "italic") # Italicize the title
  )

p2 <- ggplot(merged_Tlei, aes(x = (size/1000000), y = DE_counts)) +
  geom_point(aes(color = color_group), size = 2) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", size = 0.5) +
  labs(x = "Scaffold size (Mb)", y = "DE gene count", title = "T. leiboldiana") +
  scale_color_manual(values = c(goldenrod = 'goldenrod', darkgreen = 'forestgreen', blue = 'royalblue', default = 'black')) +
  scale_alpha_manual(values = ifelse(merged_Tlei$color_group == 'default', 1, 0.5)) + # Define alpha for each group
  theme(
    legend.position = "none",
    plot.title = element_text(face = "italic") # Italicize the title
  )

pdf("FigureS16_DEgenecount_ScaffoldSize.pdf", width = 8, height = 4)
grid.arrange(p1,p2, ncol = 2)
dev.off()
