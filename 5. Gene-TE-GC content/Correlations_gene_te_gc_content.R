if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
setwd("/Users/clara/Documents/GitHub.nosync/Tillandsia-compgenomics/I. Genome overview/A-Circular-Figure/")

genic_content <- read.table("Gene_counts_per_1MB_windows.Tfas-Tlei.mainScaffolds.curatedOGs.txt", header = T)
genic_content_Tfas <- genic_content[startsWith(genic_content$chrom, "Tfas"),]
genic_content_Tlei <- genic_content[startsWith(genic_content$chrom, "Tlei"),]

TE_content <- read.table("TE_content_Tfas-Tlei_per1MB-window_python.txt", header = T)
TE_content_Tfas <- TE_content[startsWith(TE_content$chrom, "Tfas"),]
TE_content_Tlei <- TE_content[startsWith(TE_content$chrom, "Tlei"),]

GC_content_Tfas <- read.table("GC_content_per1000000window_Tfas.txt", header = T)
GC_content_Tlei <- read.table("GC_content_per1000000window_Tlei.txt", header = T)

### IN T. FASCICULATA ###

all_Tfas <- cbind(genic_content_Tfas, TE_content_Tfas$perc_te, GC_content_Tfas$perc_GC)
colnames(all_Tfas) <- c("chrom", "start_window", "end_window",
                        "gene_counts", "TE_content", "GC_content")

# Shapiro-Wilk normality test for gene count
shapiro.test(all_Tfas$gene_counts) # => p-value < 2.2e-16
# Shapiro-Wilk normality test for te content
shapiro.test(all_Tfas$TE_content) # => p-value < 2.2e-16
# Shapiro-Wilk normality test for te content
shapiro.test(all_Tfas$GC_content) # => p-value < 2.2e-16

# Our data is not normally distributed, so we use a non-parametric 
# test for correlation

# Kendall correlation test between genic and repetitive content
res <- cor.test(all_Tfas$gene_counts, all_Tfas$TE_content,  method="kendall")
res # correlation coefficient = -0.7866569, p-value < 2.2e-16

# Kendall correlation test between genic and GC content
res <- cor.test(all_Tfas$gene_counts, all_Tfas$GC_content,  method="kendall")
res # correlation coefficient = -0.6815129, p-value < 2.2e-16

# Kendall correlation test between repetitive and GC content
res <- cor.test(all_Tfas$TE_content, all_Tfas$GC_content,  method="kendall")
res # correlation coefficient = 0.7930836, p-value < 2.2e-16

### IN T. LEIBOLDIANA ###

all_Tlei <- cbind(genic_content_Tlei, TE_content_Tlei$perc_te, GC_content_Tlei$perc_GC)
colnames(all_Tlei) <- c("chrom", "start_window", "end_window",
                        "gene_counts", "TE_content", "GC_content")

# Shapiro-Wilk normality test for gene count
shapiro.test(all_Tlei$gene_counts) # => p-value < 2.2e-16
# Shapiro-Wilk normality test for te content
shapiro.test(all_Tlei$TE_content) # => p-value < 2.2e-16
# Shapiro-Wilk normality test for te content
shapiro.test(all_Tlei$GC_content) # => p-value < 2.2e-16

# Our data is not normally distributed, so we use a non-parametric 
# test for correlation

# Kendall correlation test between genic and repetitive content
res <- cor.test(all_Tlei$gene_counts, all_Tlei$TE_content,  method="kendall")
res # correlation coefficient = -0.8220984, p-value < 2.2e-16

# Kendall correlation test between genic and GC content
res <- cor.test(all_Tlei$gene_counts, all_Tlei$GC_content,  method="kendall")
res # correlation coefficient = -0.7079335, p-value < 2.2e-16

# Kendall correlation test between repetitive and GC content
res <- cor.test(all_Tlei$TE_content, all_Tlei$GC_content,  method="kendall")
res # correlation coefficient = 0.7776291, p-value < 2.2e-16
