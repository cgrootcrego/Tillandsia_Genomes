#!/usr/bin/Rscript --vanilla

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("tidyverse", "rstatix", "stringr", "ggplot2")

args <- commandArgs(trailingOnly = TRUE)

#setwd('/Users/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/te_intersect/')
#setwd('/home/clara/Documents/GitHub/Tillandsia-compgenomics/7. Rna-seq experiment 6 timepoints/Co-expression_networks_MaSigPro/mapped_to_Tfas/te_intersect/')


#counts <- read.table('Tillandsia_fasciculata_GENE-TE-intersection.ALLgenes.counts.edited.txt', sep = '\t', header = TRUE)
counts <- read.table(args[[1]], sep = "\t", header = TRUE)
species <- paste0(str_split(args[[1]], "_")[[1]][1], " ", str_split(args[[1]], "_")[[1]][2])
geneset <- dim(counts)[1]
#colnames(counts) <- c('chr', 'start', 'end','name', 'te_counts')
#genes <- scan('List_DE_Genes_Tfas.txt', character())
genes <- scan(args[[2]], character())
DE_counts <- subset(counts, counts$name %in% genes)
DEset <- dim(DE_counts)[1]
counts$type <- ifelse(counts$name %in% genes, "subset", "total")

cat(paste0("\nTE INSERTION STATS\nSpecies: ", species, "\nTotal gene set: ", geneset, "\nDE genes: ", DEset, "\n\n"))

# Normalize by gene length, or intronic length when available
if (grepl("intron", as.character(args[[1]]), fixed = TRUE) == TRUE){
	cat("...working with intronic lengths...\n")
	counts$norm_te_counts <- counts$te_counts / counts$intronic_length * 100
	na <- dim(counts[is.na(counts$norm_te_counts),])[1]
	inf <- dim(counts[is.infinite(counts$norm_te_counts),])[1]
	remove <- na + inf
	cat(paste0("Removing ", remove, " genes with no introns\n\n"))
	counts <- counts[!is.na(counts$norm_te_counts),]
	counts <- counts[!is.infinite(counts$norm_te_counts),]
	med_DE <- median(counts[counts$type == "subset",10])
	med_non_DE <- median(counts[counts$type == "total",10])
	# Test whether there are length differences between DE and non-DE genes
	res_length <- wilcox.test(intronic_length ~ type, data = counts,
	                          exact = FALSE)
	med_length_DE <- median(counts[counts$type == "subset","intronic_length"])
	med_length_nonDE <- median(counts[counts$type == "total","intronic_length"])
	cat(paste0("\nMann-Whitney U test of length differences between DE and non-DE genes:\np-value: ",
	           res_length$p.value,"\nmedian intronic length DE genes:",med_length_DE,"\nmedian intronic length non-DE genes:", med_length_nonDE, "\n\n"))
} else {
	cat("...working with full transcript lengths...\n")
	counts$length <- counts$end - counts$start
	counts$norm_te_counts <- counts$te_counts / counts$length * 100
	med_DE <- median(counts[counts$type == "subset",8])
	med_non_DE <- median(counts[counts$type == "total",8])
	# Test whether there are length differences between DE and non-DE genes
	res_length <- wilcox.test(length ~ type, data = counts,
	                          exact = FALSE)
	med_length_DE <- median(counts[counts$type == "subset","length"])
	med_length_nonDE <- median(counts[counts$type == "total","length"])
	cat(paste0("\nMann-Whitney U test of length differences between DE and non-DE genes:\np-value: ",
	           res_length$p.value,"\nmedian length DE genes:",med_length_DE,"\nmedian length non-DE genes:", med_length_nonDE, "\n\n"))
}

# Calculate whether differentially expressed genes have on average higher rates of intronic TEs than
# non-DE genes using the Mann-Whitney U test
res <- wilcox.test(norm_te_counts ~ type, data = counts,
                   exact = FALSE)

res2 <- wilcox.test(te_counts ~ type, data = counts,
                   exact = FALSE)
med_DE2 <- median(counts[counts$type == "subset",5])
med_non_DE2 <- median(counts[counts$type == "total",5])
mea_DE2 <- mean(counts[counts$type == "subset",5])
mea_non_DE2 <- mean(counts[counts$type == "total",5])

cat(paste0("Mann-Whitney U test of intronic TE insertions between DE and non-DE genes:\np-value: ",
             res$p.value,"\nInsertion median DE genes:",med_DE,"\nInsertion median non-DE genes:", med_non_DE, "\n"))

cat(paste0("\nWITHOUT LENGTH CORRECTION:\nMann-Whitney U test of intronic TE insertions between DE and non-DE genes:\np-value: ",
           res2$p.value,"\nInsertion median DE genes:",med_DE2,"\tmean:", mea_DE2,"\nInsertion median non-DE genes:",
           med_non_DE2, "\tmean:", mea_non_DE2,"\n"))

# Calculate whether the percentage of genes with a TE insertion is significantly higher in DE genes
contingency <- data.frame(rbind(c(sum(counts$te_counts == 0),
                                    sum(counts[counts$type == "subset",]$te_counts == 0)),
                                  c(sum(counts$te_counts != 0),
                                    sum(counts[counts$type == "subset",]$te_counts != 0))))
colnames(contingency) <- c("total", "subset")
rownames(contingency) <- c("no_te", "te")

chisq <- chisq.test(contingency)

cat(paste0("\nChi-square test of proportion of genes with a TE insertion, DE subset versus all genes:","\n","Contingency table:\n"))
print(contingency)
cat(paste0("\n", "p-value: ", chisq$p.value, "\n"))

#ggplot(counts, aes(x=te_counts, color=type)) +
#  geom_histogram(fill="white", binwidth = .5) + xlim(c(0,20))

#ggplot(counts, aes(x=te_counts, color=type)) +
#  geom_density() + xlim(c(0,15))
