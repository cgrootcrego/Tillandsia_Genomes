# TL - 111022
setwd("/media/tleroy/Data/dell_mtp/Collab/Clara/genespace")
library("GENESPACE")
library(RColorBrewer)
library("viridis")

####
setwd("/media/tleroy/Data/dell_mtp/Collab/Clara/genespace/TillandsiaGenespace/")
runwd <- file.path("/media/tleroy/Data/dell_mtp/Collab/Clara/genespace/TillandsiaGenespace")
#make_exampleDataDir(writeDir = runwd)
list.files(runwd, recursive = T, full.names = F)

gpar <- init_genespace(
  genomeIDs = c("Acomosus","Tfasciculata","Tleiboldiana"),
  speciesIDs = c("Acomosus","Tfasciculata","Tleiboldiana"),
  versionIDs = c("Acomosus","Tfasciculata","Tleiboldiana"),
  ploidy = rep(1,3),
  diamondMode = "fast",
  orthofinderMethod = "fast",
  wd = runwd,
  nCores = 4,
  minPepLen = 50,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "/home/tleroy/Softwares/OrthoFinder/orthofinder",
  path2mcscanx = "/home/tleroy/Softwares/MCScanX-master",
  rawGenomeDir = file.path(runwd, "rawGenomes"))


parse_annotations(
  gsParam = gpar,
  gffEntryType = "gene",
  gffIdColumn = "locus",
  gffStripText = "locus=",
  headerEntryIndex = 1,
  headerSep = " ",
  headerStripText = "locus=",
  troubleshoot = TRUE )


gpar <- run_orthofinder(
  gsParam = gpar)


gpar <- synteny(gsParam = gpar)

# indicate the names of the chromosomes to invert for the plot
chrinvert=read.table("chromosomes_to_invert_riparian_plot.txt",h=F)

ripdat <- plot_riparian(
  gpar,
  colByChrs = c(RColorBrewer::brewer.pal(13, "Set1"),RColorBrewer::brewer.pal(13, "Set2")),
  #chrFill="darkgrey",
  #blackBg=FALSE,
  #invertTheseChrs=data.frame(genome = "Tleiboldiana", chr = "Tlei_chr16"))
  invertTheseChrs=data.frame(genome = chrinvert$V1, chr = chrinvert$V2))
  #colByChrs = c("#BC4F43", "#F67243"))

ripdat <- plot_riparian(
  gpar,
  colByChrs = c(RColorBrewer::brewer.pal(8, "Set1"),RColorBrewer::brewer.pal(8, "Set2"),RColorBrewer::brewer.pal(9, "Set1")),
                #braidAlpha=0.6,
                #braidBorderLwd=1.6,
                refGenome="Acomosus",
                genomeIDs=c("Tleiboldiana","Tfasciculata","Acomosus"),
                #highlightRef="coral",
                chrLabCex=0.6,
                excludeChrOutOfRegion=TRUE,
                #chrBorder = "darkorange",
                #gapProp = 0,
                minGenes2plot=200,
                returnSourceData=TRUE,
                #chrFill="darkgrey",
                blackBg=FALSE,
                #invertTheseChrs=data.frame(genome = "Tleiboldiana", chr = "Tlei_chr16"))
                invertTheseChrs=data.frame(genome = chrinvert$V1, chr = chrinvert$V2))
  #colByChrs = c("#BC4F43", "#F67243"))
  



