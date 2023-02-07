#!/usr/bin/Rscript --vanilla

setwd("/Users/clara/Documents/GitHub.nosync/Tillandsia-compgenomics/Supplementary_figures/V.Expression_curves/")

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("ggplot2", "reshape2", "stringr","grid", "gridExtra",
               "RColorBrewer", "forcats")

# Load arguments
# 1 is counts, 2 is the subset of genes of each module, 3 is the GOterms
args <- commandArgs(trailingOnly = TRUE)

# Load data
#counts <- read.table("counts.Tfas_Tlei_6_timepoints.exons.sum.normalized-cpm.EdgeR.logtransformed.mediancentered.txt", header = T, row.names = 1)
counts <- read.table(args[[1]], header = T, row.names = 1)
genes <- scan(args[[2]], character(), quote = "")
#genes <- scan("Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic-cluster1.txt", character(), quote = "")
specgen <- read.table(args[[3]], header = T, sep = "\t")
#specgen <- read.table("genes_of_interest.txt", header = T, sep = "\t")

module <- str_split(args[[2]], "\\-|\\_|\\.")[[1]][12]
#module <- str_split("Genes_Significant_Tfas-vs-Tlei_0.7-trimmed_TLEI-REF.exonic-cluster1.txt", "\\-|\\_|\\.")[[1]][12]
mod <- gsub('cluster','',module)

# Make colors for functional categories
#nb.cols <- 11
#mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
#mycolors <- as.data.frame(cbind(unique(specgen$Category), mycolors))
#colnames(mycolors) <- c("Category", "Color")
mycolors <- read.table(args[[4]], header = T, sep = "\t")

# Obtain genes belonging to modules, split by species
module_counts <- subset(counts, rownames(counts) %in% genes)
module_counts_Tfas <- module_counts[, c(1:36)]
module_counts_Tlei <- module_counts[, c(37:72)]
mod_size <- as.character(nrow(module_counts))

mean_count_Tfas <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tfas[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tfas) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")
mean_count_Tlei <- as.data.frame(sapply(seq(1, 6, 1), function(j) rowMeans(module_counts_Tlei[, c(j,j+6,j+12,j+18,j+24,j+30)])))
colnames(mean_count_Tlei) <- c("N+9", "D+1", "D+5","D+9", "N+1", "N+5")

mean_count_Tfas$gene_id <- row.names(mean_count_Tfas)
mean_count_Tlei$gene_id <- row.names(mean_count_Tlei)

mean_count_Tfas_m <- melt(mean_count_Tfas, id.vars = "gene_id")
colnames(mean_count_Tfas_m) <- c("gene_id", "time", "count")
mean_count_Tlei_m <- melt(mean_count_Tlei, id.vars = "gene_id")
colnames(mean_count_Tlei_m) <- c("gene_id", "time", "count")
write.table(mean_count_Tfas_m, file = paste0("MedianCentered_Mean_Expression_counts_T.fasciculata-", module, ".txt"))
write.table(mean_count_Tlei_m, file = paste0("MedianCentered_Mean_Expression_counts_T.leiboldiana-", module, ".txt"))

specgen_in_mod <- specgen[specgen$Cluster == mod,]
cats <- unique(specgen[specgen$Cluster == mod,"Category"])

specgen_list_Tfas <- list()
specgen_list_Tlei <- list()
for (cat in cats){
  g <- specgen_in_mod[specgen_in_mod$Category == cat, "gene_id"]
  expr_Tfas <- subset(mean_count_Tfas_m, gene_id %in% g)
  expr_Tlei <- subset(mean_count_Tlei_m, gene_id %in% g)
  expr_Tfas$Category <- cat
  expr_Tlei$Category <- cat
  specgen_list_Tfas <- append(specgen_list_Tfas, list(expr_Tfas))
  specgen_list_Tlei <- append(specgen_list_Tlei, list(expr_Tlei))
}

# Calculate scale
max_Tfas <- max(mean_count_Tfas_m$count)
min_Tfas <- min(mean_count_Tfas_m$count)
max_Tlei <- max(mean_count_Tlei_m$count)
min_Tlei <- min(mean_count_Tlei_m$count)

if (max_Tfas > max_Tlei){
  ymax = max_Tfas + 1
} else {
  ymax = max_Tlei + 1
}

if (min_Tfas < min_Tlei){
  ymin = min_Tfas - 1
} else {
  ymin = min_Tlei - 1
}

# Calculate total mean curve
means_Tfas <- data.frame(time=c("N+9", "D+1", "D+5","D+9", "N+1", "N+5"),
                         count = c(mean(mean_count_Tfas$`N+9`), mean(mean_count_Tfas$`D+1`),
                                   mean(mean_count_Tfas$`D+5`), mean(mean_count_Tfas$`D+9`),
                                   mean(mean_count_Tfas$`N+1`), mean(mean_count_Tfas$`N+5`)))

means_Tlei <- data.frame(time=c("N+9", "D+1", "D+5","D+9", "N+1", "N+5"),
                         count = c(mean(mean_count_Tlei$`N+9`), mean(mean_count_Tlei$`D+1`),
                                   mean(mean_count_Tlei$`D+5`), mean(mean_count_Tlei$`D+9`),
                                   mean(mean_count_Tlei$`N+1`), mean(mean_count_Tlei$`N+5`)))

mean_count_Tfas_m$order <- c(rep(3,dim(mean_count_Tfas)[1]),
                             rep(4,dim(mean_count_Tfas)[1]),
                             rep(5,dim(mean_count_Tfas)[1]),
                             rep(6,dim(mean_count_Tfas)[1]),
                             rep(1,dim(mean_count_Tfas)[1]),
                             rep(2,dim(mean_count_Tfas)[1]))
mean_count_Tfas_m$time <- fct_reorder(mean_count_Tfas_m$time, mean_count_Tfas_m$order)

mean_count_Tlei_m$order <- c(rep(3,dim(mean_count_Tlei)[1]),
                             rep(4,dim(mean_count_Tlei)[1]),
                             rep(5,dim(mean_count_Tlei)[1]),
                             rep(6,dim(mean_count_Tlei)[1]),
                             rep(1,dim(mean_count_Tlei)[1]),
                             rep(2,dim(mean_count_Tlei)[1]))
mean_count_Tlei_m$time <- fct_reorder(mean_count_Tlei_m$time, mean_count_Tlei_m$order)

# Highlight genes of interest
pdf(paste("Expression_curve_", module, "_logtransformed.mediancentered.pdf", sep = ""), width = 12, height = 8)
p1 <- ggplot(mean_count_Tfas_m, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means_Tfas, aes(group = 1), size = 1, color = "black",
            linetype = "dashed") +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  ylim(c(ymin, ymax)) +
  ylab("Median-centered log(CPM)") +
  xlab("Time") +
  theme( axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = 'black', size = 1),
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

for (i in 1:length(cats)){
  cat <- cats[i]
  p1 <- p1 + geom_line(data = specgen_list_Tfas[[i]], aes(group = gene_id),
                       size = 1, color = mycolors[mycolors$Category == cat, "Color"])
}

grob <- grobTree(textGrob("T. fasciculata", x=0.02,  y=0.97, hjust=0,
  gp=gpar(col="black", fontsize=15, fontface="italic")))
p1 <- p1 + annotation_custom(grob)

p2 <- ggplot(mean_count_Tlei_m, aes(x=time, y=count, group = gene_id)) +
  geom_point(color = "grey") +
  geom_line(color = "grey") +
  geom_line(data = means_Tlei, aes(group = 1), size = 1, color = "black",
            linetype = "dashed") +
  geom_vline(xintercept = 3.5, linetype = "dashed") +
  ylim(ymin, ymax) +
  ylab("Median-centered log(CPM)") +
  xlab("Time") +
  theme( axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1), axis.ticks = element_line(colour = 'black', size = 1),
        axis.text.x = element_text(colour="black",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=13,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="black",size=15,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=15,angle=90,hjust=.5,vjust=.5,face="italic"))

for (i in 1:length(cats)){
  cat <- cats[i]
  p2 <- p2 + geom_line(data = specgen_list_Tlei[[i]], aes(group = gene_id),
                       size = 1, color = mycolors[mycolors$Category == cat, "Color"])
}

grob <- grobTree(textGrob("T. leiboldiana", x=0.02,  y=0.97, hjust=0, gp=gpar(col="black", fontsize=15, fontface="italic", fill = "white")))
p2 <- p2 + annotation_custom(grob)

grid.arrange(p1, p2, nrow = 1, top = textGrob(paste0("CLUSTER ",mod, " - ", mod_size, " genes"), gp=gpar(fontsize=22,font=8)))
dev.off()
