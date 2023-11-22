# 17.10.2023 
# This script plots malate abundances and performs statistical analyses
pacman::p_load("ggplot2", "tidyverse", "dplyr", "grid", "gridExtra", "cowplot", "patchwork", "ggpubr", "factoextra")

setwd('/home/botanik/Documents/GitHub/Tillandsia-compgenomics/8. Metabolic_measurements/')
setwd('/Users/clara/Documents/GitHub.nosync/Tillandsia-compgenomics/8. Metabolic_measurements/')

dat <- read.delim('Tillandsia_MTIC_normalized.csv', header = T, sep = ';')
dat_r <- read.delim('Tillandsia_metabolomics_full.csv', header = T, sep = ';')

dat_MTIC <- dat[dat$Speciescode=="F" | dat$Speciescode=="L",]
dat_raw <- dat_r[dat_r$Speciescode=="F" | dat_r$Speciescode=="L",]

dat_full <- merge(dat_raw, dat_MTIC, by = "Code")
# Identify duplicate columns
duplicates <- duplicated(as.list(dat_full))
# Remove duplicate columns
dat_full <- dat_full[!duplicates]

dat <- dat_full %>% select(Code, Speceis_indiv.x, Speciescode.x, Species.x,	Individual.x,	Timepoint.x,	MTIC.x, Malic.acid..3TMS._MTIC, Malic.acid..3TMS., Malic.acid..3TMS._PE)
dat$Timepoint.x <- as.factor(dat$Timepoint.x)
dat$Timepoint.x <- factor(dat$Timepoint.x, levels = c("N1", "N5", "N9", "D1", "D5", "D9"))

colnames(dat) <- c("code", "sp_ind",	"sp_code",	"species",	"ind",	"timepoint",	"MTIC", "malate_MTIC",	"malate_raw",	"malate_PE")

p_all <- ggplot(dat, aes(
  x=timepoint,
  y=malate_MTIC, 
  fill=species)) +
  geom_boxplot(width = 0.6, color = "black", alpha = .8, outlier.shape = NA)+
  geom_point(aes(fill = species, alpha = .9),color = "black", shape = 21, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  theme_bw() +
  labs(title="",
       x ="Time", y = "Peak area of Malate, MTIC-normalized")+
  theme(axis.title = element_text(colour = "black", size = 10), 
        axis.text.x = element_text(colour="black",size=10,face="plain"),
        axis.text.y = element_text(colour="black",size=10,face="plain")) +
  scale_fill_manual(values = c("T. fasciculata" = "goldenrod", "T. leiboldiana"="darkgreen")) +
  scale_color_manual(values = c("T. fasciculata" = "goldenrod", "T. leiboldiana"="darkgreen")) +
  guides(alpha=FALSE) +
  #ylim(c(-0.0005,0.011)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")+
  scale_x_discrete(labels = c("D1" = "D+1", "D5" = "D+5", "D9" = "D+9", "N1" = "N+1", "N5" = "N+5", "N9" = "N+9"))

### Take extreme differences per sample, not per timepoint
# Calculate the most extreme differences for each sample
extreme_differences <- dat %>%
  group_by(ind, species) %>%
  mutate(min_time = timepoint[which.min(malate_MTIC)],
         max_time = timepoint[which.max(malate_MTIC)]) %>%
  filter(timepoint %in% c(min_time, max_time)) %>%
  summarise(Malate_diff = abs(diff(malate_MTIC)),
            min_timepoint = min_time,
            max_timepoint = max_time) %>%
  distinct()

# Calculate the average difference for each species
extreme_diff_avg <- extreme_differences %>%
  group_by(species) %>%
  summarise(mean_diff = mean(Malate_diff),
            sd_diff = sd(Malate_diff))

# Perform Wilcoxon rank-sum test
wilcoxon_test <- wilcox.test(Malate_diff ~ species, data = extreme_differences)

p_delta <- ggplot(extreme_diff_avg, aes(x=species, y=mean_diff, fill=species)) +
  theme_bw() +
  geom_bar(stat="identity", color="black", alpha = .8, 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_diff-sd_diff, ymax=mean_diff+sd_diff), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("T. fasciculata" = "goldenrod", "T. leiboldiana" = "darkgreen")) +
  labs(title="", 
       y = "Mean difference in Malate Abundance",
       fill="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.x = element_blank(),
        legend.position = "bottom")

p_delta <- ggplot(extreme_differences, aes(x=species, y=Malate_diff, fill=species)) +
  theme_bw() +
  geom_boxplot(width = 0.6, color = "black", alpha = .8, outlier.shape = NA)+
  geom_point(aes(fill = species, alpha = .9),color = "black", shape = 21, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  theme_bw()+
  scale_fill_manual(values = c("T. fasciculata" = "goldenrod", "T. leiboldiana" = "darkgreen")) +
  labs(title="", 
       y = "Diel accumulation of Malate",
       fill="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=10),
        axis.text.x = element_blank(),
        legend.position = "bottom")+
  guides(alpha=FALSE) 

p_delta
# Combining the timewise boxplots + delta values
legend <- get_legend(p_all)
p_all_n <- p_all + theme(legend.position = "none")
p_delta_n <- p_delta + theme(legend.position = "none")

# Create the plot
plot <- grid.arrange(p_all_n, grid.arrange(p_delta_n, legend, nrow = 2, heights = c(6, 1)), ncol = 2, widths = c(3.25, 1.75))

### PCA on all metabolic measurements ###
dat_pca_full <- read.csv("Tillandsia_metabolomics_full.SUMMED.csv", header = T, sep = ";")
summary(dat_pca_full)
dat_pca <- dat_pca_full[,c(17:65)]

row.names(dat_pca) <- dat_pca_full$Code
dat_pca_full$species_timepoint <- interaction(dat_pca_full$Species, dat_pca_full$Timepoint)
dat_pca_full <- dat_pca_full[, c("species_timepoint", names(dat_pca_full)[-ncol(dat_pca_full)])]

dat_pca_full$species_timepoint <- factor(dat_pca_full$species_timepoint, levels = c("T. fasciculata.D1", "T. fasciculata.D5", "T. fasciculata.D9", "T. fasciculata.N1", "T. fasciculata.N5", "T. fasciculata.N9", "T. leiboldiana.D1", "T. leiboldiana.D5", "T. leiboldiana.D9", "T. leiboldiana.N1", "T. leiboldiana.N5", "T. leiboldiana.N9"))

res.pca <- prcomp(dat_pca, scale = TRUE)
# Extract individual level results
res.ind <- get_pca_ind(res.pca)

explained_var <- summary(res.pca)$importance["Proportion of Variance", ]

# Access and scale the loadings
loadings <- as.data.frame(res.pca$rotation)

pca_plot <- ggplot(dat_pca_full, aes(x = -res.ind$coord[, 1], y = -res.ind$coord[, 2])) +
  theme_bw() +
  geom_point(aes(fill=species_timepoint), shape = 21, color = "black", size = 3) +
  scale_fill_manual(labels = c("T. fasciculata, D+1", "T. fasciculata, D+5", 
                               "T. fasciculata, D+9", "T. fasciculata, N+1", 
                               "T. fasciculata, N+5", "T. fasciculata, N+9", 
                               "T. leiboldiana, D+1", "T. leiboldiana, D+5",
                               "T. leiboldiana, D+9", "T. leiboldiana, N+1", 
                               "T. leiboldiana, N+5", "T. leiboldiana, N+9"),
                    values = c(
                      "T. fasciculata.D1" = "lightgoldenrod1", 
                      "T. fasciculata.D5" = "lightgoldenrod2", 
                      "T. fasciculata.D9" = "lightgoldenrod3", 
                      "T. fasciculata.N1" = "darkgoldenrod1", 
                      "T. fasciculata.N5" = "darkgoldenrod2", 
                      "T. fasciculata.N9" = "darkgoldenrod3", 
                      "T. leiboldiana.D1" = "darkolivegreen1", 
                      "T. leiboldiana.D5" = "darkolivegreen2",
                      "T. leiboldiana.D9" = "darkolivegreen3", 
                      "T. leiboldiana.N1" = "darkgreen", 
                      "T. leiboldiana.N5" = "#004B00", 
                      "T. leiboldiana.N9" = "#003100")) +
  xlab(paste0("PC1: ", round(explained_var[1] * 100, 2), "% variance")) +
  ylab(paste0("PC2: ", round(explained_var[2] * 100, 2), "% variance")) +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=10))+
  labs(fill = "Species, timepoint")
dat_pca$Glu
pca_plot
# Select the CAM metabolites
metabolites <- c("Malic.acid", "Citric.acid", "Fumaric.acid", "Glucose","Fructose", "Fructose.1.6.diphosphate", "Sucrose", "Oxaloacetate", "Succinic.acid", "Maltose", "Glucose.6.phosphate", "Phenylalanine", "Threonine", "Gluconic.acid", "Pyruvic.acid", "Proline","Raffinose", "Trehalose")
selected_rows <- sapply(rownames(loadings), function(x) any(sapply(metabolites, grepl, x = x)))
selected_loadings <- loadings[selected_rows, ]

new_names <- c("Citrate", "Fructose", "F-1,6-DP", "Fumarate", "Gluconic acid", "Glucose", "G-6-P", "Malate", "Maltose", "Oxaloacetate", "Phenylalalnine", "Proline", "Pyruvic acid","Raffinose", "Succinate", "Sucrose", "Threonine", "Trehalose")

rownames(selected_loadings) <- new_names

# Define specific coordinates for each label
scaling_factor <- max(res.ind$coord) / max(loadings) * 1.8
label_coordinates <- data.frame(
  Label = rownames(selected_loadings),
  x = c(-2.3, -1.0188139,	0.5, -1.5, -1.7, 0.2, -1.8, 1.2, -2.3, 	
        -1.6, -1.8),
  y = c(.8, 2.1,0.75, 0.3, 1.85, -0.5, 1.5, 0.05, -1.3, -2.49812001, 	
        -1.97161838)
)

# Define the offsets
label_offset_x <- .5
label_offset_y <- .5

# Apply the offsets conditionally
selected_loadings$label_offset_x <- ifelse(selected_loadings$PC1 > 0, -label_offset_x, label_offset_x)
selected_loadings$label_offset_y <- ifelse(selected_loadings$PC2 > 0, -label_offset_y, label_offset_y)

selected_loadings["F-1,6-DP", "label_offset_x"] <- selected_loadings["F-1,6-DP", "label_offset_x"] - .5
selected_loadings["Oxaloacetate", "label_offset_x"] <- selected_loadings["Oxaloacetate", "label_offset_x"] + 1

pca_plot_m <- pca_plot +
  geom_segment(data = as.data.frame(selected_loadings), 
               aes(x = 0, y = 0, xend = scaling_factor * -PC1, yend = scaling_factor * -PC2),
               arrow = arrow(length = unit(0.2, "cm")), 
               alpha = 0.5, color = "blue") +
  geom_text(data = as.data.frame(selected_loadings), 
            aes(x = scaling_factor * -PC1 + label_offset_x, y = scaling_factor * -PC2 + label_offset_y, label = rownames(selected_loadings)), 
            size = 3, vjust = 1, hjust = 0.5)

pca_plot_m

combined_plots <- grid.arrange(arrangeGrob(pca_plot_m), arrangeGrob(plot), ncol = 1, heights = c(3, 2))
ggsave("Figure2_Metabolomics_boxplot_dots.09112023.png", plot = combined_plots, height = 8, width = 9)
ggsave("Figure2_Metabolomics_boxplot_dots.09112023.pdf", plot = combined_plots, height = 8, width = 9)

loadings <- read.table("pca_loadings_summed_controlle.csv", header = T, sep = ";")

# Order Variable by descending order of PC1
loadings <- loadings %>%
  arrange(desc(PC1)) %>%
  mutate(Variable = factor(Variable, levels = Variable))

# Create the plot
loadings_pc1 <- ggplot(loadings, aes(x = Variable, y = PC1, color = Category)) +
  geom_point(size = 3) +  # Draw the circles
  geom_segment(aes(x = Variable, xend = Variable, y = 0, yend = PC1), size = 1) +  # Draw the lines
  scale_color_manual(labels= c("Aminoacid", "Marker", "Organic Acid", "Sugar"),values = c("#EDC79B", "#713E5A", "#58A4B0", "#CC5803"))+  # Use a color palette
  ylab("Loadings of PC1") +
  xlab("") +
  coord_flip() +  # Make the plot horizontal
  theme_minimal()  # Use a minimal theme

# Order Variable by descending order of PC2
loadings <- loadings %>%
  arrange(desc(PC2)) %>%
  mutate(Variable = factor(Variable, levels = Variable))

# Create the plot
loadings_pc2 <- ggplot(loadings, aes(x = Variable, y = PC2, color = Category)) +
  geom_point(size = 3) +  # Draw the circles
  geom_segment(aes(x = Variable, xend = Variable, y = 0, yend = PC2), size = 1) +  # Draw the lines
  scale_color_manual(labels= c("Aminoacid", "Marker", "Organic Acid", "Sugar"),values = c("#EDC79B", "#713E5A", "#58A4B0", "#CC5803"))+  # Use a color palette
  ylab("Loadings of PC2") +
  xlab("") +
  coord_flip() +  # Make the plot horizontal
  theme_minimal()  # Use a minimal theme

legend <- get_legend(loadings_pc1)
loadings_pc1_nl <- loadings_pc1 + theme(legend.position = "none")
loadings_pc2_nl <- loadings_pc2 + theme(legend.position = "none")
loadings_plot <- grid.arrange(loadings_pc1_nl, loadings_pc2_nl, legend, ncol = 3, widths = c(3, 3, 1))

ggsave("FigureS1_Loadings_PC1-2.pdf", plot = loadings_plot, height = 9, width = 9)
ggsave("FigureS1_Loadings_PC1-2.png", plot = loadings_plot, height = 9, width = 9)
