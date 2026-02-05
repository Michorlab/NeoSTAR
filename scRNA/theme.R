# Author:       Stephen-John Sammut
# Description:  ggplot theme for manuscript
#=========================================================================

library(pacman)

suppressPackageStartupMessages(
  p_load(data.table,
         ggplot2,
         ggpubr,
         grid,
         ggthemes,
         pheatmap,
         RColorBrewer,
         reshape2,
         tidyverse,
         viridis)
)

theme_manuscript <- function(base_size=12, base_family="arial") {
  (ggthemes::theme_foundation(base_size=base_size)
   + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = base_size),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title.y = element_text(angle = 90, vjust = 2, size = base_size),
           axis.title.x = element_text(vjust = -0.2, size = base_size),
           axis.text = element_text(size = base_size), 
           axis.line = element_line(colour = "black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour = "#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size = unit(0.5, "cm"),
           legend.spacing = unit(0, "cm"),
           strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
           strip.text = element_text(size = base_size)
   ))
}


#==========================================
# COLOUR SCHEMES
#==========================================


colours.tnbc <- list()

colours.tnbc$IntClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#8B0000', '#FFFF40', '#0000CD', '#FFAA00',
                           '#EE82EE', '#593C8F',"#AF5B5B")
names(colours.tnbc$IntClust) <- c(1:3,"4-",5:10,"Other")

colours.tnbc$TnbcType4        <- c("#DE1A1A","firebrick4","#F6AE2D","#006494","#006494","#ACADBC")
names(colours.tnbc$TnbcType4) <-c("BL1","BL2","LAR","MES","Mesenchymal","NC")

colours.tnbc$PAM50        <- c("#DE1A1A", "#FB9A99", "#1F78B4", "#66A61E", "#CBA328")
names(colours.tnbc$PAM50) <- c("Basal", "Her2", "Luminal", "Normal", "Claudin low")

colours.tnbc$express  <- c("#375E97","#FB6542")
names(colours.tnbc$express) <-c("E4","E10")

colours.tnbc$cell_colors <-c("#386641","#6A994E","#A7C957","#984EA3","#ff90b3","#FFC2E2","#BA2D0B", "#54B0E4","#377EB8", "#222F75","#FFC971", "#FB8500")
names(colours.tnbc$cell_colors) <- c('B cells','CD4 T cells','CD8 T cells','NK','Myeloid','Mast cells','Endothelial','iCAFs','myCAFs','Other CAFs','Other epithelial','Tumour cells')

colours.tnbc$All <- c(colours.tnbc$IntClust,colours.tnbc$TnbcType4,colours.tnbc$PAM50,colours.tnbc$express)



