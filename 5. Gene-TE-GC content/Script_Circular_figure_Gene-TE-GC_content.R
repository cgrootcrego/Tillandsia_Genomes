#!/usr/bin/Rscript --vanilla

# Code for large circular figure - Genome Paper T.fasciculata / T. leiboldiana
# This Rcode creates a circular figure with the 25 / 24 main scaffolds of both assemblies as halves of the circle. The circle has 3 tracks: Gene count, TE and GC content

# Installation and loading
# Pacman will only install missing packages
if (!require("pacman")) install.packages("pacman", repos = "
http://cran.us.r-project.org")
pacman::p_load("circlize", "stringr", "RColorBrewer")

#setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/I. Circular figure")

# Load arguments
# 1 is chromosome list, 2 is the gene density, 3 is TE content, 4 is GC content, 5 is output name
args <- commandArgs(trailingOnly = TRUE)
output_name <- args[[5]]
# Read in the complete list of chromosomes with all needed info:
chrom <- read.table(args[[1]], header = T, sep = "\t")
# Make matrix of start and end position to initialize the circular plot
pos <- cbind(rep(0, dim(chrom)[1]), chrom$size)
# Make plot
pdf(paste0(output_name, ".pdf"), width = 10, height = 8)
print("Initializing the plot...")
circos.clear()
# Plot initialization
circos.par("track.height"=0.05, cell.padding=c(0.02, 0, 0.02, 0),
           start.degree=90,gap.degree=c(rep(1.5, 24), 6, rep(1.5, 23), 6),points.overflow.warning=FALSE,
           track.margin=c(0,0), circle.margin = c(.1,.1,.1,.1))
circos.initialize(sectors = chrom$name, xlim = pos)

# Add blocks representing the chromosomes
nb.cols <- 25
mycolors <- c()
mycolors_tfas <- c("#E8D430","#7F5477","#E97421","#F781BF","#00afb9","#B47229", "#4C71A4","#50AA4D", "#FFAF13", "#419485", "#B13749", "#ccd5ae", "#FFD422", "#8ecae6", "#E41A1C", "#C76766", "#cdb4db", "#FF8904", "#91569B", "#FFF930", "#DF7492", "#CEA32D", "#ce4257", "#a1c181", "#3A86A5")
mycolors_tlei <- c(mycolors_tfas[4], mycolors_tfas[1], mycolors_tfas[13],
                   mycolors_tfas[15], mycolors_tfas[12],mycolors_tfas[3],
                   mycolors_tfas[9], mycolors_tfas[8], mycolors_tfas[6],
                   mycolors_tfas[5], mycolors_tfas[14], mycolors_tfas[16],
                   mycolors_tfas[10], mycolors_tfas[17], mycolors_tfas[11],
                   mycolors_tfas[18], mycolors_tfas[7], mycolors_tfas[19],
                   mycolors_tfas[2], mycolors_tfas[21], mycolors_tfas[23],
                   mycolors_tfas[22], mycolors_tfas[24], mycolors_tfas[20])
mycolors <- c(mycolors_tfas, mycolors_tlei)
circos.track(chrom$name, ylim = c(0, 1), bg.border = mycolors, bg.col = mycolors, track.height = .05)

# Add chromosome names
is.even <- function(x) x %% 2 == 0
n = 0
m = 0
for (i in 1:nrow(chrom)){
  name=chrom[i,1]
  if (is.even(i) == TRUE){
    if (chrom[i,3] == "Tfas") {
      circos.text(chrom$size/2 - n, 2.3, str_split(name, "_sca")[[1]][2], sector.index=name,
                  col="grey40",cex=0.5, facing = "inside", niceFacing = T)
      n = n + 220000
      if (i > 21 & i < 26){
        n = n + 600000
      }
    } else {
      circos.text(chrom$size - n, 2.3, str_split(name, "_sca")[[1]][2], sector.index=name,
                  col="darkgreen",cex=0.5, facing = "inside", niceFacing = T)
      n = n + 20000
      if (i > 38){
        n = n + 2400000
      }
    }
  } else {
    if ((i > 21 & i < 26) | i > 45){
      if (chrom[i,3] == "Tfas") {
        circos.text(chrom$size/2 - n, 3.4, str_split(name, "_sca")[[1]][2], sector.index=name,
                    col="grey40",cex=0.5, facing = "inside", niceFacing = T)
        n = n + 600000
      } else {
        circos.text(chrom$size - n, 3.4, str_split(name, "_sca")[[1]][2], sector.index=name,
                    col="darkgreen",cex=0.5, facing = "inside", niceFacing = T)
        n = n + 600000
      }
    } else {
      if (chrom[i,3] == "Tfas") {
        circos.text(chrom$size/2 - n , 2.3, str_split(name, "_sca")[[1]][2], sector.index=name,
                    col="grey40",cex=0.5, facing = "inside", niceFacing = T)
        n = n + 220000
      } else {
        circos.text(chrom$size - n, 2.3, str_split(name, "_sca")[[1]][2], sector.index=name,
                    col="darkgreen",cex=0.5, facing = "inside", niceFacing = T)
        n = n + 20000
        if (i > 38){
          n = n + 2400000
        }
      }
    }
  }
}

# Add species lines
draw.sector(90, 301, rou1 = 1.09, rou2 = 1.10, col = "grey40", border="NA")
draw.sector(297, 95, rou1 = 1.09, rou2 = 1.10, col = "darkgreen", border="NA")
circos.text(chrom$size/2, 9, "T. fasciculata", sector.index="Tfas_sca11",col="grey40",cex=1.5, facing = "bending.inside")
circos.text(chrom$size/2, 9, "T. leiboldiana", sector.index="Tlei_sca11",col="darkgreen",cex=1.5, facing = "bending.inside")

#-------------------TRACK 1: GENE DENSITY-------------------#

## Read in gene content files
# All gene counts
#gene_counts_per_mb_windows <- read.table("Gene_counts_per_1MB_windows.Tfas-Tlei.mainScaffolds.curatedOGs.txt.no-Tlei_sca2526", header = T)
gene_counts_per_mb_windows <- read.table(args[[2]], header = T)
print("Drawing first track: Gene density...")
color_gene_density <- c(Tfas_sca1 = "palegreen3", Tfas_sca2 = "palegreen3", Tfas_sca3 = "palegreen3",
                        Tfas_sca4 = "palegreen3", Tfas_sca5 = "palegreen3", Tfas_sca6 = "palegreen3",
                        Tfas_sca7 = "palegreen3", Tfas_sca8 = "palegreen3", Tfas_sca9 = "palegreen3",
                        Tfas_sca10 = "palegreen3", Tfas_sca11 = "palegreen3", Tfas_sca12 = "palegreen3",
                        Tfas_sca13 = "palegreen3", Tfas_sca14 = "palegreen3", Tfas_sca15 = "palegreen3",
                        Tfas_sca16 = "palegreen3", Tfas_sca17 = "palegreen3", Tfas_sca18 = "palegreen3",
                        Tfas_sca19 = "palegreen3", Tfas_sca20 = "palegreen3", Tfas_sca21 = "palegreen3",
                        Tfas_sca22 = "palegreen3", Tfas_sca23 = "palegreen3", Tfas_sca24 = "palegreen3",
                        Tfas_sca25 = "palegreen3", Tlei_sca1 = "seagreen", Tlei_sca2 = "seagreen",
                        Tlei_sca3 = "seagreen", Tlei_sca4 = "seagreen", Tlei_sca5 = "seagreen",
                        Tlei_sca6 = "seagreen", Tlei_sca7 = "seagreen", Tlei_sca8 = "seagreen",
                        Tlei_sca9 = "seagreen", Tlei_sca10 = "seagreen", Tlei_sca11 = "seagreen",
                        Tlei_sca12 = "seagreen", Tlei_sca13 = "seagreen", Tlei_sca14 = "seagreen",
                        Tlei_sca15 = "seagreen",Tlei_sca16 = "seagreen", Tlei_sca17 = "seagreen",
                        Tlei_sca18 = "seagreen",Tlei_sca19 = "seagreen", Tlei_sca20 = "seagreen",
                        Tlei_sca21 = "seagreen",Tlei_sca22 = "seagreen", Tlei_sca23 = "seagreen",
                        Tlei_sca24 = "seagreen")
circos.track(gene_counts_per_mb_windows$chrom, y = gene_counts_per_mb_windows$gene_counts,
             x = gene_counts_per_mb_windows$start_window,
             bg.col = "grey92", panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = color_gene_density[CELL_META$sector.index])
               circos.yaxis(side = "left", sector.index = "Tfas_sca1", tick = T, tick.length = 3, labels = T,
                            at = c(0,50, 100), labels.cex = 0.40, labels.col="gray28")
             }, track.height = 0.15, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 50)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}
#-------------------TRACK 3: TE DENSITY-------------------#

#TE_content_per_mb_windows <- read.table("TE_content_Tfas-Tlei_per1MB-window_python.txt.no-Tlei_sca2526", header = T, sep = "\t")
TE_content_per_mb_windows <- read.table(args[[3]], header = T,sep = "\t")
print("Drawing second track: TE density...")

color_te_density <- c(Tfas_sca1 = "#f7d486", Tfas_sca2 = "#f7d486", Tfas_sca3 = "#f7d486",
                        Tfas_sca4 = "#f7d486", Tfas_sca5 = "#f7d486", Tfas_sca6 = "#f7d486",
                        Tfas_sca7 = "#f7d486", Tfas_sca8 = "#f7d486", Tfas_sca9 = "#f7d486",
                        Tfas_sca10 = "#f7d486", Tfas_sca11 = "#f7d486", Tfas_sca12 = "#f7d486",
                        Tfas_sca13 = "#f7d486", Tfas_sca14 = "#f7d486", Tfas_sca15 = "#f7d486",
                        Tfas_sca16 = "#f7d486", Tfas_sca17 = "#f7d486", Tfas_sca18 = "#f7d486",
                        Tfas_sca19 = "#f7d486", Tfas_sca20 = "#f7d486", Tfas_sca21 = "#f7d486",
                        Tfas_sca22 = "#f7d486", Tfas_sca23 = "#f7d486", Tfas_sca24 = "#f7d486",
                        Tfas_sca25 = "#f7d486", Tlei_sca1 = "#edae49", Tlei_sca2 = "#edae49",
                        Tlei_sca3 = "#edae49", Tlei_sca4 = "#edae49", Tlei_sca5 = "#edae49",
                        Tlei_sca6 = "#edae49", Tlei_sca7 = "#edae49", Tlei_sca8 = "#edae49",
                        Tlei_sca9 = "#edae49", Tlei_sca10 = "#edae49", Tlei_sca11 = "#edae49",
                        Tlei_sca12 = "#edae49", Tlei_sca13 = "#edae49", Tlei_sca14 = "#edae49",
                        Tlei_sca15 = "#edae49",Tlei_sca16 = "#edae49", Tlei_sca17 = "#edae49",
                        Tlei_sca18 = "#edae49",Tlei_sca19 = "#edae49", Tlei_sca20 = "#edae49",
                        Tlei_sca21 = "#edae49",Tlei_sca22 = "#edae49", Tlei_sca23 = "#edae49",
                        Tlei_sca24 = "#edae49")

circos.track(TE_content_per_mb_windows$chrom, x = TE_content_per_mb_windows$start_window,
             y = TE_content_per_mb_windows$perc_te, bg.col = "grey92", ylim = c(0,100),
             panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = color_te_density[CELL_META$sector.index])
               circos.yaxis(c("left"), sector.index = "Tfas_sca1", tick = T, tick.length = 3,
                            labels = T, at = seq(0, CELL_META$ylim[2], by = 25),
                            labels.cex = 0.40, labels.col="gray28")
             }, track.height = 0.17, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(0, CELL_META$ylim[2], by = 25)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

#-------------------TRACK 3: GC Content-------------------#

GC_content_per_mb_windows <- read.table(args[[4]], header = T,sep = "\t")
#GC_content_per_mb_windows <- read.table("GC_content_per1000000window-noTlei-chr2526.txt", header = T,sep = "\t")
print("Drawing third track: GC content...")

color_te_density <- c(Tfas_sca1 = "#f27a7d", Tfas_sca2 = "#f27a7d", Tfas_sca3 = "#f27a7d",
                        Tfas_sca4 = "#f27a7d", Tfas_sca5 = "#f27a7d", Tfas_sca6 = "#f27a7d",
                        Tfas_sca7 = "#f27a7d", Tfas_sca8 = "#f27a7d", Tfas_sca9 = "#f27a7d",
                        Tfas_sca10 = "#f27a7d", Tfas_sca11 = "#f27a7d", Tfas_sca12 = "#f27a7d",
                        Tfas_sca13 = "#f27a7d", Tfas_sca14 = "#f27a7d", Tfas_sca15 = "#f27a7d",
                        Tfas_sca16 = "#f27a7d", Tfas_sca17 = "#f27a7d", Tfas_sca18 = "#f27a7d",
                        Tfas_sca19 = "#f27a7d", Tfas_sca20 = "#f27a7d", Tfas_sca21 = "#f27a7d",
                        Tfas_sca22 = "#f27a7d", Tfas_sca23 = "#f27a7d", Tfas_sca24 = "#f27a7d",
                        Tfas_sca25 = "#f27a7d", Tlei_sca1 = "#d1495b", Tlei_sca2 = "#d1495b",
                        Tlei_sca3 = "#d1495b", Tlei_sca4 = "#d1495b", Tlei_sca5 = "#d1495b",
                        Tlei_sca6 = "#d1495b", Tlei_sca7 = "#d1495b", Tlei_sca8 = "#d1495b",
                        Tlei_sca9 = "#d1495b", Tlei_sca10 = "#d1495b", Tlei_sca11 = "#d1495b",
                        Tlei_sca12 = "#d1495b", Tlei_sca13 = "#d1495b", Tlei_sca14 = "#d1495b",
                        Tlei_sca15 = "#d1495b",Tlei_sca16 = "#d1495b", Tlei_sca17 = "#d1495b",
                        Tlei_sca18 = "#d1495b",Tlei_sca19 = "#d1495b", Tlei_sca20 = "#d1495b",
                        Tlei_sca21 = "#d1495b",Tlei_sca22 = "#d1495b", Tlei_sca23 = "#d1495b",
                        Tlei_sca24 = "#d1495b")

circos.track(GC_content_per_mb_windows$chrom, x = GC_content_per_mb_windows$start_window,
             y = GC_content_per_mb_windows$perc_GC, bg.col = "grey92", ylim = c(30,50),
             panel.fun = function(x, y) {
               circos.lines(x, y, area = T, col = color_te_density[CELL_META$sector.index])
               circos.yaxis(c("left"), sector.index = "Tfas_sca1", labels = T, tick = T, tick.length = 3,
                            at = seq(30, 50, by = 10),
                            labels.cex = 0.40, labels.col="gray28")
             }, track.height = 0.17, bg.border = "black")
for(sn in get.all.sector.index()) {
  set.current.cell(sector.index = sn, track.index = get.current.track.index())
  breaks = seq(30, 50, by = 10)
  for(b in breaks) {
    circos.lines(CELL_META$cell.xlim, rep(b, 2), lty = 3, col = "#00000040")
  }
}

dev.off()
print("Done!")
