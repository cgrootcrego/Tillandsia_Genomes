# 24.11.2023
### Per-gene family expression heatmaps

pacman::p_load("ggplot2", "reshape2", "RColorBrewer", "tidyr", "ComplexHeatmap", "circlize")

#setwd('/Users/clara/Documents/GitHub.nosync/Tillandsia-CAM-Drought/2. DGE/VI. Gene-specific Expression/')
#setwd('/home/botanik/Documents/GitHub/Tillandsia-CAM-Drought/2. DGE/VI. Gene-specific Expression/')

# 1. Expression data 2. sample info 3. Genes of interest
args <- commandArgs(trailingOnly = TRUE)

#cpm <- read.delim("../Tfas_gene22.counts.NormalizedTPM.All.csv", sep = "\t", header = T)
cpm <- read.delim(args[[1]], sep = "\t", header = T)

names(cpm) <- sub("^X", "", names(cpm))

#info <- read.table("Sampling_info.txt", header = T)
info <- read.table(args[[2]], header = T)
colnames(info) <- c("ID", "species", "Timepoint")

#genes <- read.table("test_genes.txt", sep = "\t", header = F)
genes <- read.table(args[[3]], sep = "\t", header = F)
colnames(genes) <- c("gene_id", "og_id", "newname", "group", "DE_in_")

# e.g. "Malate dehydrogenase"
descr <- args[[4]]

intersected_genes <- intersect(rownames(cpm), genes$gene_id)
genes_cpm <- as.data.frame((cpm[intersected_genes, ]))

# Create the matrix for the heatmap
 genes_cpm_log10 <- log10(genes_cpm+.1)
# heatmap_matrix <- as.matrix(genes_cpm_log10)
z_scores <- t(scale(t(genes_cpm_log10)))
heatmap_matrix <- as.matrix(z_scores)

# Create divergent color gradient for z-scores
midpoint <- 0
# Define the minimum and maximum of the breaks, making sure the range is symmetric around the midpoint
max_value <- max(abs(heatmap_matrix))
breaks <- seq(-max_value, max_value, length.out = 11)
# Get a divergent color palette from RColorBrewer
colors <- brewer.pal(11, "RdBu")
# Reverse the colors to match the order of the breaks (red for positive, blue for negative)
colors <- rev(colors)
# Define color mapping function with a color ramp
color_mapping <- colorRamp2(breaks, colors)

# Define color mapping for annotations
timepoint_colors <- c("D+1" = "#fdcd70", "D+5" = "#fcbf49", "D+9" = "#df9404", "N+1" = "#76415E", "N+5" = "#a75f86", "N+9" = "#d1abbf")

# Create annotation
# Reorder matrix by timepoint
info <- info[order(info$Timepoint),]
# Then, use the ordered sample names to reorder the columns of the heatmap matrix
heatmap_matrix <- heatmap_matrix[, info$ID]

info$species <- factor(info$species, levels = c("Tfas", "Tlei"))
info$Timepoint <- factor(info$Timepoint, levels = c("D+1", "D+5", "D+9", "N+1", "N+5", "N+9"))

# Create a factor vector specifying the grouping of columns
# Map original species labels to new ones
info$new_species <- factor(info$species,
                           levels = c("Tfas", "Tlei"),
                           labels = c("T. fasciculata", "T. leiboldiana"))

# Replace the row names of the heatmap_matrix with the new gene names
rownames(heatmap_matrix) <- genes$newname[match(rownames(heatmap_matrix), genes$gene_id)]
heatmap_matrix <- heatmap_matrix[order(rownames(heatmap_matrix)), ]

genes <- genes[match(rownames(heatmap_matrix), genes$newname), ]
genes_group_factor <- as.factor(genes$group)

# Create a color vector for the gene names
gene_name_colors <- ifelse(genes$DE_in_ == "DE", "darkgreen","grey40")
names(gene_name_colors) <- genes$newname

ha <- HeatmapAnnotation(df = data.frame(Timepoint = info$Timepoint),
                        col = list(Timepoint = timepoint_colors),
                        show_legend = TRUE,
                        show_annotation_name = FALSE)

rownames(heatmap_matrix) <- paste0("         ", rownames(heatmap_matrix))

ht <- Heatmap(heatmap_matrix,
              name = "z-score log(TPM)",
              col = color_mapping,
              top_annotation = ha,
              cluster_rows = FALSE,     # Disable row clustering
              cluster_columns = FALSE, # Disable column clustering
              show_column_names = FALSE,
              show_row_names = TRUE,  # Show row names (gene names)
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 6, col = gene_name_colors),
              column_split = info$new_species,
              row_split = genes_group_factor,
              row_title_gp = gpar(fontsize = 8),
              row_title_rot = 45,
              heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")))

pdf(file = paste0(descr,".pdf"), width = 4, height = 4)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

jpeg(file = paste0(descr,".jpg"), res = 300, width = 1800, height = 1800)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
