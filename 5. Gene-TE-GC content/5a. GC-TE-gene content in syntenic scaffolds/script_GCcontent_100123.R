# TL - 201022
# GC content Clara
library(ggplot2)
library(grid)
library(gridExtra)
setwd("~/Collab/Clara/GCcontent/")

#### TE & GC%

Acos=read.table("GCF_001540865.1_ASM154086v1_genomic.fna.sed_TE-GC.content_slwindows100kb_newversion_pourR.sed",sep="\t",header = TRUE)
Tfas=read.table("Tillandsia_fasciculata_25_scaffolds.fasta.masked.soft_TE-GC.content_slwindows100kb_newversion_pourR.clean",sep="\t",header = TRUE)
Tlei=read.table("Tillandsia_leiboldiana_26_scaffolds.fasta.masked.soft_TE-GC.content_slwindows100kb_newversion_pourR.clean",sep="\t",header = TRUE)

pdf("TE-GCcontent_summary_3firstchrs_of_eachsp.pdf",width = 14,height = 9)

Acoschr1=subset(Acos,Acos$sequence=="1")
plotAcos1 <- ggplot(Acoschr1)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr1 - A. comosus")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


Acoschr2=subset(Acos,Acos$sequence=="2")
plotAcos2 <- ggplot(Acoschr2)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr2 - A. comosus")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

Acoschr3=subset(Acos,Acos$sequence=="3")
plotAcos3 <- ggplot(Acoschr3)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr3 - A. comosus")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotAcos4 <- ggplot(Acoschr1)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr1 - A. comosus")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotAcos5 <- ggplot(Acoschr2)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr2 - A. comosus")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotAcos6 <- ggplot(Acoschr3)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr3 - A. comosus")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plotAcos1, vp = define_region(row = 1, col = 1))   # Span over two columns
print(plotAcos2, vp = define_region(row = 1, col = 2))
print(plotAcos3, vp = define_region(row = 1, col = 3))
print(plotAcos4, vp = define_region(row = 2, col = 1))   # Span over two columns
print(plotAcos5, vp = define_region(row = 2, col = 2))
print(plotAcos6, vp = define_region(row = 2, col = 3))


# Tfas
Tfaschr1=subset(Tfas,Tfas$sequence=="Tfas_chr1")
plotTfas1 <- ggplot(Tfaschr1)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr1 - T. fasciculata")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


Tfaschr2=subset(Tfas,Tfas$sequence=="Tfas_chr2")
plotTfas2 <- ggplot(Tfaschr2)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr2 - T. fasciculata")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

Tfaschr3=subset(Tfas,Tfas$sequence=="Tfas_chr3")
plotTfas3 <- ggplot(Tfaschr3)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr3 - T. fasciculata")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotTfas4 <- ggplot(Tfaschr1)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr1 - T. fasciculata")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotTfas5 <- ggplot(Tfaschr2)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr2 - T. fasciculata")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotTfas6 <- ggplot(Tfaschr3)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr3 - T. fasciculata")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plotTfas1, vp = define_region(row = 1, col = 1))   # Span over two columns
print(plotTfas2, vp = define_region(row = 1, col = 2))
print(plotTfas3, vp = define_region(row = 1, col = 3))
print(plotTfas4, vp = define_region(row = 2, col = 1))   # Span over two columns
print(plotTfas5, vp = define_region(row = 2, col = 2))
print(plotTfas6, vp = define_region(row = 2, col = 3))


# Tlei
Tleichr1=subset(Tlei,Tlei$sequence=="Tlei_chr1")
plotTlei1 <- ggplot(Tleichr1)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr1 - T. leiboldiana")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))


Tleichr2=subset(Tlei,Tlei$sequence=="Tlei_chr2")
plotTlei2 <- ggplot(Tleichr2)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr2 - T. leiboldiana")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

Tleichr3=subset(Tlei,Tlei$sequence=="Tlei_chr3")
plotTlei3 <- ggplot(Tleichr3)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr3 - T. leiboldiana")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotTlei4 <- ggplot(Tleichr1)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr1 - T. leiboldiana")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotTlei5 <- ggplot(Tleichr2)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr2 - T. leiboldiana")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plotTlei6 <- ggplot(Tleichr3)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr3 - T. leiboldiana")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plotTlei1, vp = define_region(row = 1, col = 1))   # Span over two columns
print(plotTlei2, vp = define_region(row = 1, col = 2))
print(plotTlei3, vp = define_region(row = 1, col = 3))
print(plotTlei4, vp = define_region(row = 2, col = 1))   # Span over two columns
print(plotTlei5, vp = define_region(row = 2, col = 2))
print(plotTlei6, vp = define_region(row = 2, col = 3))


dev.off()


####   Syntenic chromosomes ########

pdf("TE-GCcontent_summary_3synthenicchromosomes.pdf",width = 14,height = 9)

# Acos6 / Tfas11 /Tlei15
Acoschr6=subset(Acos,Acos$sequence=="6")
Tfaschr11=subset(Tfas,Tfas$sequence=="Tfas_chr11")
Tleichr15=subset(Tlei,Tlei$sequence=="Tlei_chr15")

plot1 <- ggplot(Acoschr6)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr6 - A. comosus")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot2 <- ggplot(Tfaschr11)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr11 - T. fasciculata")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot3 <- ggplot(Tleichr15)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr15 - T. leiboldiana")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))



plot4 <- ggplot(Acoschr6)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr6 - A. comosus")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot5 <- ggplot(Tfaschr11)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr11 - T. fasciculata")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot6 <- ggplot(Tleichr15)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr15 - T. leiboldiana")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plot1, vp = define_region(row = 1, col = 1))   # Span over two columns
print(plot2, vp = define_region(row = 1, col = 2))
print(plot3, vp = define_region(row = 1, col = 3))
print(plot4, vp = define_region(row = 2, col = 1))   # Span over two columns
print(plot5, vp = define_region(row = 2, col = 2))
print(plot6, vp = define_region(row = 2, col = 3))




# Acos3 / Tfas4 /Tlei1
Acoschr3=subset(Acos,Acos$sequence=="3")
Tfaschr4=subset(Tfas,Tfas$sequence=="Tfas_chr4")
Tleichr1=subset(Tlei,Tlei$sequence=="Tlei_chr1")

plot1 <- ggplot(Acoschr3)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr3 - A. comosus")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot2 <- ggplot(Tfaschr4)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr4 - T. fasciculata")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot3 <- ggplot(Tleichr1)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr1 - T. leiboldiana")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))



plot4 <- ggplot(Acoschr3)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr3 - A. comosus")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot5 <- ggplot(Tfaschr4)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr4 - T. fasciculata")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot6 <- ggplot(Tleichr1)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr1 - T. leiboldiana")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plot1, vp = define_region(row = 1, col = 1))   # Span over two columns
print(plot2, vp = define_region(row = 1, col = 2))
print(plot3, vp = define_region(row = 1, col = 3))
print(plot4, vp = define_region(row = 2, col = 1))   # Span over two columns
print(plot5, vp = define_region(row = 2, col = 2))
print(plot6, vp = define_region(row = 2, col = 3))




# Acos11 / Tfas12 /Tlei5
Acoschr11=subset(Acos,Acos$sequence=="11")
Tfaschr12=subset(Tfas,Tfas$sequence=="Tfas_chr12")
Tleichr5=subset(Tlei,Tlei$sequence=="Tlei_chr5")

plot1 <- ggplot(Acoschr11)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr11 - A. comosus")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot2 <- ggplot(Tfaschr12)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr12 - T. fasciculata")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot3 <- ggplot(Tleichr5)+
  geom_point(aes(x=start_pos_windows,y=GCrate),size=1.5,color="red")+
  geom_smooth(aes(x=start_pos_windows,y=GCrate),span=0.2,color="red")+
  geom_point(aes(x=start_pos_windows,y=TErate),size=1.5,color="navyblue")+
  geom_smooth(aes(x=start_pos_windows,y=TErate),span=0.2,color="navyblue")+
  ylim(0,1)+
  xlab("Position chr5 - T. leiboldiana")+
  ylab("TE% (blue),GC%(red)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))



plot4 <- ggplot(Acoschr11)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr11 - A. comosus")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot5 <- ggplot(Tfaschr12)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr12 - T. fasciculata")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot6 <- ggplot(Tleichr5)+
  geom_point(aes(x=start_pos_windows,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(aes(x=start_pos_windows,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(aes(x=start_pos_windows,y=GCrateinTE),span=0.2,color="forestgreen")+
  ylim(0.2,0.6)+
  xlab("Position chr5 - T. leiboldiana")+
  ylab("GC% in non-TE (black), in TE (green)")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plot1, vp = define_region(row = 1, col = 1))   # Span over two columns
print(plot2, vp = define_region(row = 1, col = 2))
print(plot3, vp = define_region(row = 1, col = 3))
print(plot4, vp = define_region(row = 2, col = 1))   # Span over two columns
print(plot5, vp = define_region(row = 2, col = 2))
print(plot6, vp = define_region(row = 2, col = 3))

dev.off()


### merged plot 
plot1tripletA <- ggplot()+
  geom_point(data=Acoschr6,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),size=1.5,color="black")+
  geom_smooth(data=Acoschr6,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),span=0.2,color="black")+
  geom_point(data=Tfaschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),span=0.2,color="grey")+
  geom_point(data=Tleichr15,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr15,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),span=0.2,color="forestgreen")+
  ylim(0,1)+
  xlab("")+
  ylab("TE%")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot2tripletA <- ggplot()+
  geom_point(data=Acoschr6,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),size=1.5,color="black")+
  geom_smooth(data=Acoschr6,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),span=0.2,color="black")+
  geom_point(data=Tfaschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),span=0.2,color="grey")+
  geom_point(data=Tleichr15,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr15,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),span=0.2,color="forestgreen")+
  #ylim(0,1)+
  xlab("")+
  ylab("GC%")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot3tripletA <- ggplot()+
  geom_point(data=Acoschr6,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),size=1.5,color="black")+
  geom_smooth(data=Acoschr6,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),span=0.2,color="black")+
  geom_point(data=Tfaschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),span=0.2,color="grey")+
  geom_point(data=Tleichr15,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr15,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),span=0.2,color="forestgreen")+
  #ylim(0,1)+
  xlab("% chromosome length (Acos6/Tfas11/Tlei15)")+
  ylab("GC% in TEs")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot4tripletA <- ggplot()+
  geom_point(data=Acoschr6,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(data=Acoschr6,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(data=Tfaschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),span=0.2,color="grey")+
  geom_point(data=Tleichr15,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr15,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),span=0.2,color="forestgreen")+
  #ylim(0,1)+
  xlab("% chromosome length (Acos6/Tfas11/Tlei15)")+
  ylab("GC% in non-TEs")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot1tripletB <- ggplot()+
  geom_point(data=Acoschr3,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),size=1.5,color="black")+
  geom_smooth(data=Acoschr3,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),span=0.2,color="black")+
  geom_point(data=Tfaschr4,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr4,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),span=0.2,color="grey")+
  geom_point(data=Tleichr1,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr1,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),span=0.2,color="forestgreen")+
  ylim(0,1)+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot2tripletB <- ggplot()+
  geom_point(data=Acoschr3,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),size=1.5,color="black")+
  geom_smooth(data=Acoschr3,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),span=0.2,color="black")+
  geom_point(data=Tfaschr4,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr4,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),span=0.2,color="grey")+
  geom_point(data=Tleichr1,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr1,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),span=0.2,color="forestgreen")+
  #ylim(0,1)+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot3tripletB <- ggplot()+
  geom_point(data=Acoschr3,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),size=1.5,color="black")+
  geom_smooth(data=Acoschr3,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),span=0.2,color="black")+
  geom_point(data=Tfaschr4,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr4,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),span=0.2,color="grey")+
  geom_point(data=Tleichr1,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr1,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),span=0.2,color="forestgreen")+
  #ylim(0,1)+
  xlab("% chromosome length (Acos3/Tfas4/Tlei1)")+
  ylab("")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot4tripletB <- ggplot()+
  geom_point(data=Acoschr3,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(data=Acoschr3,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(data=Tfaschr4,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr4,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),span=0.2,color="grey")+
  geom_point(data=Tleichr1,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr1,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),span=0.2,color="forestgreen")+
  #ylim(0,1)+
  xlab("% chromosome length (Acos3/Tfas4/Tlei1)")+
  ylab("GC% in non-TEs")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot1tripletC <- ggplot()+
  geom_point(data=Acoschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),size=1.5,color="black")+
  geom_smooth(data=Acoschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),span=0.2,color="black")+
  geom_point(data=Tfaschr12,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr12,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),span=0.2,color="grey")+
  geom_point(data=Tleichr5,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr5,aes(x=start_pos_windows/max(start_pos_windows)*100,y=TErate),span=0.2,color="forestgreen")+
  ylim(0,1)+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot2tripletC <- ggplot()+
  geom_point(data=Acoschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),size=1.5,color="black")+
  geom_smooth(data=Acoschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),span=0.2,color="black")+
  geom_point(data=Tfaschr12,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr12,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),span=0.2,color="grey")+
  geom_point(data=Tleichr5,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr5,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrate),span=0.2,color="forestgreen")+
  #ylim(0,1)+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot3tripletC <- ggplot()+
  geom_point(data=Acoschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),size=1.5,color="black")+
  geom_smooth(data=Acoschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),span=0.2,color="black")+
  geom_point(data=Tfaschr12,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr12,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),span=0.2,color="grey")+
  geom_point(data=Tleichr5,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr5,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinTE),span=0.2,color="forestgreen")+
  #ylim(0,1)+
  xlab("% chromosome length (Acos11/Tfas12/Tlei5)")+
  ylab("")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

plot4tripletC <- ggplot()+
  geom_point(data=Acoschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),size=1.5,color="black")+
  geom_smooth(data=Acoschr11,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),span=0.2,color="black")+
  geom_point(data=Tfaschr12,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),size=1.5,color="grey")+
  geom_smooth(data=Tfaschr12,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),span=0.2,color="grey")+
  geom_point(data=Tleichr5,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),size=1.5,color="forestgreen")+
  geom_smooth(data=Tleichr5,aes(x=start_pos_windows/max(start_pos_windows)*100,y=GCrateinNonTE),span=0.2,color="forestgreen")+
  #ylim(0,1)+
  xlab("% chromosome length (Acos11/Tfas12/Tlei5)")+
  ylab("")+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=14,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="italic"))

grid.newpage()
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 3)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(plot1tripletA, vp = define_region(row = 1, col = 1))   
print(plot2tripletA, vp = define_region(row = 2, col = 1))
print(plot3tripletA, vp = define_region(row = 3, col = 1))
#print(plot4tripletA, vp = define_region(row = 4, col = 1))   
print(plot1tripletB, vp = define_region(row = 1, col = 2))   
print(plot2tripletB, vp = define_region(row = 2, col = 2))
print(plot3tripletB, vp = define_region(row = 3, col = 2))
#print(plot4tripletB, vp = define_region(row = 4, col = 2))
print(plot1tripletC, vp = define_region(row = 1, col = 3))   
print(plot2tripletC, vp = define_region(row = 2, col = 3))
print(plot3tripletC, vp = define_region(row = 3, col = 3))
#print(plot4tripletC, vp = define_region(row = 4, col = 3))  



### more ancient investigations: see script:  script_GCcontent_201022.R
