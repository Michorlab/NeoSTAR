# ── Install required packages ────────────────────────────────────────────────

## CRAN packages
install.packages(c(
  # Tidyverse stack
  "tidyverse", "patchwork",
  # Visualization
  "circlize", "pheatmap", "RColorBrewer", "viridis",
  "ggrepel", "ggsignif", "ggpmisc", "ggpubr", "ggExtra",
  "ggnewscale", "ggsci", "ggalluvial", "gridExtra",
  # Enrichment / gene set
  "msigdbr",
  # Bulk / signatures / misc
  "gtools", "doMC",
  # NMF / misc
  "reshape2", "scales", "cowplot", "ggcorrplot",
  # Seurat dependencies
  "Seurat"
))

install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

## Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "SingleCellExperiment",
  "celldex",
  "sccomp",
  "MAST",
  "ComplexHeatmap",
  "clusterProfiler",
  "org.Hs.eg.db",
  "GSEABase",
  "GSVA",
  "edgeR",
  "preprocessCore",
  "glmGamPoi",
  "infercnv",
  "BiocNeighbors"
))

## GitHub packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("immunogenomics/harmony")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
remotes::install_github("cole-trapnell-lab/monocle3")
remotes::install_github("samuel-marsh/scCustomize")
remotes::install_github("sqjin/CellChat")
remotes::install_github("bhklab/genefu")


install.packages("s2")
install.packages("sf")
install.packages("rjags", force = TRUE)
remotes::install_github("cole-trapnell-lab/monocle3")

# using homebrew
Sys.setenv(HDF5_DIR = "/opt/homebrew/opt/hdf5")
Sys.setenv(PKG_CONFIG_PATH = "/opt/homebrew/opt/hdf5/lib/pkgconfig")
remotes::install_github("cole-trapnell-lab/monocle3")


# Retry CellChat (BiocNeighbors is now present)
remotes::install_github("sqjin/CellChat")

# Retry monocle3 after hdf5 is installed
remotes::install_github("cole-trapnell-lab/monocle3")

## NMF (may need special handling)
install.packages("NMF")
install.packages("rstatix")
## genefu (Bioconductor)
BiocManager::install("genefu")

R.version.string

# set and confirm path
setwd("/Users/isabellapabon/Documents/Ting Data")

# Library all packages needed for scRNAseq analysis
## Core
library(Seurat)
library(SingleCellExperiment)

## Tidyverse stack
library(tidyverse)     # includes: dplyr, ggplot2, tidyr, stringr, forcats, etc.
library(patchwork)

## Single-cell / integration / annotation
library(harmony)
library(celldex)
library(scCustomize)
library(sccomp)
library(DoubletFinder)
library(monocle3)
library(infercnv)
library(MAST)
library(rstatix)


## Visualization
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(ggsignif)
library(ggpmisc)
library(ggpubr)
library(ggExtra)
library(ggnewscale)
library(ggsci)
library(ggalluvial)
library(gridExtra)

## Enrichment / gene set
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(GSVA)
library(msigdbr)

## Bulk / signatures / misc
library(edgeR)
library(preprocessCore)
library(genefu)
library(gtools)
library(doMC)

# Addition
library(dplyr)
library(ggplot2)
library(ggsci)

library(RColorBrewer)
library(viridis)


# library sccomp
library(dplyr)
library(sccomp)
library(ggplot2)
library(forcats)
library(tidyr)

#NMF
library(reshape2)
library(NMF)
library(ggplot2)
library(scales)
library(glmGamPoi)
library(cowplot)

library(scales)
library(ggcorrplot)

setwd("/Users/isabellapabon/Documents/Ting Data")

# ======== Fig.1C sample info heatmap ========
sample<-read.csv("sample_info.csv",header = T)
head(sample)
sample$SG_pCR<-factor(sample$SG_pCR,levels=c("Y","N"))
df<-sample[order(sample$SG_pCR,sample$New_ID),]
df$WES_blood[df$WES_blood=="Y"]<-"pre"
df$WES_pre[df$WES_pre=="Y"]<-"pre"
df$WES_post[df$WES_post=="Y"]<-"post"
df$scRNA_pre[df$scRNA_pre=="Y"]<-"pre"
df$scRNA_post[df$scRNA_post=="Y"]<-"post"
df$Orion_pre[df$Orion_pre=="Y"]<-"pre"
df$Orion_post[df$Orion_post=="Y"]<-"post"

#sample matrix
exp_mat<-df[,8:14]
rownames(exp_mat)<-df$New_ID
#patient annotation
column_ha<- HeatmapAnnotation(
  df = data.frame("SG_pCR"=df$SG_pCR),
  col = list(SG_pCR = c("Y"="#CC0C00FF","N"="#5C88DAFF")#,
             #SG_RCB_value = colorRamp2(c(0, 4.4), c("white", "purple"))
  )
)
#set cell parameter
cell_fun = function(j, i, x, y, width, height, fill) {
  # draw dot
  grid.circle(
    x = x, y = y,
    r = min(width, height) * 1.8,
    gp = gpar(fill = fill, col = NA)
  )
  # draw horizontal line below rows 1, 3, 5, 7
  if (i %in% c(1, 3, 5, 7)) {
    grid.lines(x = unit(c(0, 1), "npc"),
               y = unit(y - height / 2, "npc"),
               gp = gpar(col = "black", lwd = 1))
  }
}
#draw
pdf("heatmap_sample_info.pdf",width=12,height=3)
hp<-ComplexHeatmap::Heatmap(t(exp_mat),
                            col = c("pre"="#F19465","post"="#94C465","N"="white"),
                            rect_gp = gpar(type = "none"),#remove default box
                            cell_fun = cell_fun,
                            #row_names_gp = gpar(fontsize = 10),
                            cluster_rows = FALSE,cluster_columns = FALSE,
                            #heatmap_width = unit(15, "npc"),
                            #width = 5,
                            #heatmap_height = unit(15, "npc")
                            #height = 15,
                            #row_km = 4,
                            border = NA,#remove border line
                            top_annotation = column_ha
)
draw(hp)
dev.off()


# ======== Fig.1G/S2A/S2B all cell UMAP and quantification ========
scRNA<-readRDS("All_cell_seurat.rds")
head(scRNA@meta.data)

scRNA$Celltype_Major <- factor(scRNA$Celltype_Major, 
                               levels=c("B cell",
                                        "Plasma cell",
                                        "T cell",
                                        "Endothelial cell",
                                        "Epithelial cell",
                                        "Mesenchymal cell",
                                        "Myeloid cell",
                                        "pDC")
                               )

#Set colors for UMAP plot
umapColor <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175', #36 colors
               '#FF7F0E', '#1F77B4', '#2CA02C', '#D62728', '#9467BD'#add 5 colors
)#36 colors

groupColor <- c("pre-pCR"="#CC0C00FF",
                "pre-RD"="#5C88DAFF",
                "post-pCR"="#84BD00FF",
                "post-RD"="#FFCD00FF")

#UMAP plot
p1<-DimPlot(scRNA, group.by = "Celltype_Major", label = T,pt.size=0.01,raster=FALSE)+
  scale_color_manual(values = umapColor)
p2<- DimPlot(scRNA, group.by = "new_id_4", label = F,pt.size=0.01,raster=FALSE)+
  scale_color_manual(values = umapColor)
p3<- DimPlot(scRNA, group.by = "new_group", label = F,pt.size=0.01,raster=FALSE)+
  scale_color_manual(values = groupColor)
plots = list(p1,p2,p3)
UMAP.plot <- wrap_plots(plots = plots, nrow=3)    
ggsave("output/all_UMAP_annotations.pdf", plot = UMAP.plot, width = 7, height = 13.5)




# Violin plot for QC
scRNA$new_id_4<-factor(scRNA$new_id_4,levels=mixedsort(levels(scRNA$new_id_4)))
# set theme for plot
theme.set2 = theme(axis.title.x=element_blank())
theme.set3 = theme(axis.text.x=element_text(vjust=1,size=8))
# set key parameters
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB")
group = "new_id_4"
# violin plot
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 +theme.set3+ NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)
ggsave("output/vlnplot_sampleID_qc.pdf", plot = violin, width = 13, height = 5)

# Dotplot for major cell type markers

## CHANGED SLOT TO LAYER
allGenes = row.names(GetAssayData(scRNA, layer="data"))


# Read markers
genegroup=read.table("all_anno_markers.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
annot.markerGenes<-rownames(genegroup)
# check these genes are in arrays
annot.markerGenes %in% allGenes
# set parameters
scRNA$Celltype_Major <- factor(scRNA$Celltype_Major, 
                               levels=c("T cell",
                                        "B cell",
                                        "Plasma cell",
                                        "Myeloid cell",
                                        "pDC",
                                        "Epithelial cell",
                                        "Mesenchymal cell",
                                        "Endothelial cell")
)
# draw
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="Celltype_Major"
)+coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
ggsave("output/Dotplot_Major_annot_markers.pdf", p, width = 7, height = 9)


# ======== Fig.1H-I Pseudobulk + TNBCtype calling + Major celltype comparison========

# pseudobulk count matrix: rows = genes, cols = samples
# Aggregate raw counts per group
pseudobulk_counts <- AggregateExpression(
  scRNA,
  group.by = "new_id_4",
  slot = "counts",
  return.seurat = FALSE
)$SCT #same result for SCT and RNA assays

write.table(pseudobulk_counts,"output/pseudobulk_counts.txt")

# get TNBC4type 80-gene centroids (quantile based)
centroid <- read.csv("Table.SX_2015CCR_TNBC_subtype.csv", header=TRUE, row.names=1)
head(centroid)
colnames(centroid)<-c("LAR","MES","BLIS","BLIA")

# Normalize pseudobulk counts
dge <- DGEList(counts=pseudobulk_counts)
dge <- calcNormFactors(dge, method="TMM")
logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=1)  # log2 CPM normalization, or Log2 Q3 normalization for GeoMx or array data

# Gene intersection
common.genes <- intersect(rownames(logCPM), rownames(centroid))
logCPM.sub <- logCPM[common.genes, ]
centroid <- centroid[common.genes, ]
cat("Common genes found:", length(common.genes), "of 80 genes\n")

# Quantile normalize to centroid space
# per-sample normalization rank scaling
rank_norm2 <- function(x) {
  z <- qnorm(rank(x) / (length(x) + 1))  # avoid Inf
  pnorm(z)  # scale to 0–1
}
logCPM.sub.qn <- apply(logCPM.sub, 2, rank_norm2)
rownames(logCPM.sub.qn) <- rownames(logCPM.sub)
colnames(logCPM.sub.qn) <- colnames(logCPM.sub)

summary(as.matrix(logCPM.sub.qn))
summary(as.matrix(centroid))

# Compute Euclidean distances
dist.mat <- sapply(colnames(centroid), function(subtype) {
  apply(logCPM.sub.qn, 2, function(sample) {
    sqrt(mean((sample - centroid[, subtype])^2))
  })
})

# Assign subtype with minimal distance
assigned <- apply(dist.mat, 1, function(x) names(which.min(x)))
min.dist <- apply(dist.mat, 1, min)

# ratio of distance to best vs 2nd-best centroid
second.min <- apply(dist.mat, 1, function(x) sort(x)[2])
confidence <- second.min / min.dist

# Permutation-based p-values
set.seed(123)  # for reproducibility
nperm <- 10000

perm.pvals <- sapply(1:ncol(logCPM.sub.qn), function(i) {
  sample.vec <- logCPM.sub.qn[, i]
  obs.dist <- min(sapply(colnames(centroid), function(subtype) {
    sqrt(mean((sample.vec - centroid[, subtype])^2))
  }))
  perm.dist <- replicate(nperm, {
    perm.sample <- sample(sample.vec)
    min(sapply(colnames(centroid), function(subtype) {
      sqrt(mean((perm.sample - centroid[, subtype])^2))
    }))
  })
  mean(perm.dist <= obs.dist)
})
names(perm.pvals) <- colnames(logCPM.sub.qn)

# Apply classification thresholds
final.subtype <- assigned
final.subtype[perm.pvals > 0.05 | min.dist > 0.35] <- "UNCLASSIFIED" #lenient distance cutoff (original 0.25)

quantile(perm.pvals)
quantile(min.dist)
quantile(confidence)

# Save output
result <- data.frame(
  Sample = colnames(logCPM.sub.qn),
  Subtype = final.subtype,
  MinDist = min.dist,
  Pvalue = perm.pvals,
  stringsAsFactors = FALSE
)
table(result$Subtype)
write.table(result, "output/TNBCtype4_subtype_assignment.txt", sep="\t", quote=FALSE, row.names=FALSE)

#project back to seurat
#update Allcell_TNBCtype
sample.group<-read.csv("sample_group.csv",header=T)
group_tmp <-sample.group$Allcell_TNBCtype
names(group_tmp) <- sample.group$new_id_4
Idents(scRNA) <- "new_id_4"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["Allcell_TNBCtype"]] <- Idents(object = scRNA)
head(scRNA@meta.data)#确认metadata信息

# stackbarplot showing major cell type proportions
meta <- scRNA@meta.data[,c("new_id_4", "Celltype_Major", "new_group", "Allcell_TNBCtype")]
meta$Celltype_Major <- factor(meta$Celltype_Major, 
                               levels=c("B cell",
                                        "Plasma cell",
                                        "T cell",
                                        "Endothelial cell",
                                        "Epithelial cell",
                                        "Mesenchymal cell",
                                        "Myeloid cell",
                                        "pDC")
)

# Order groups
meta$new_group <- factor(meta$new_group,
                         levels = c("pre-pCR","pre-RD","post-pCR","post-RD"))

# Compute per-sample relative composition of MP.type
df <- meta %>%
  dplyr::count(new_id_4, new_group, Allcell_TNBCtype, Celltype_Major, name = "n") %>%
  dplyr::group_by(new_id_4) %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::ungroup()

# Order samples by new_group
sample_order <- df %>%
  distinct(new_id_4, new_group) %>%
  arrange(new_group) %>%
  pull(new_id_4)
df$new_id_4 <- factor(df$new_id_4, levels = sample_order)

# Define colors
group_colors <- c("pre-pCR"="#CC0C00FF","pre-RD"="#5C88DAFF",
                  "post-pCR"="#84BD00FF","post-RD"="#FFCD00FF")
tnbc_colors <- c("BLIA"="skyblue","BLIS"="tomato",
                 "MES"="seagreen","LAR"="orange")
df$Allcell_TNBCtype<-factor(df$Allcell_TNBCtype,levels=c("BLIA","BLIS","MES","LAR"))

# Prepare top annotation data
annotation_df <- df %>%
  distinct(new_id_4, new_group, Allcell_TNBCtype)

# To create separate legends, we need to use `ggnewscale`
p <- ggplot() +
  # main stacked bar (Celltype_Major)
  geom_bar(data = df,
           aes(x = new_id_4, y = freq, fill = Celltype_Major),
           stat = "identity", width = 0.8) +
  scale_fill_manual(name = "Cell.type", values = umapColor) +
  # new scale for new_group
  ggnewscale::new_scale_fill() +
  geom_tile(data = annotation_df,
            aes(x = new_id_4, y = 1.05, fill = new_group),
            width = 0.8, height = 0.05, inherit.aes = FALSE) +
  scale_fill_manual(name = "new_group", values = group_colors) +
  # new scale for TNBCtype
  ggnewscale::new_scale_fill() +
  geom_tile(data = annotation_df,
            aes(x = new_id_4, y = 1.11, fill = Allcell_TNBCtype),
            width = 0.8, height = 0.05, inherit.aes = FALSE) +
  scale_fill_manual(name = "Allcell_TNBCtype", values = tnbc_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right") +
  labs(x = "Sample ID", y = "Fraction")
ggsave("output/StackBar_sample_proportion.pdf", plot = p, width = 10, height = 5)

# Compare TNBCtype proportions across groups
scRNA<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD","post-RD"))
meta <- scRNA@meta.data[,c("new_id_4", "Celltype_Major", "new_group", "Allcell_TNBCtype")]
meta$new_group <- factor(meta$new_group,
                         levels = c("pre-pCR","pre-RD","post-RD"))

# Summarize sample counts per TNBC subtype per new_group
tnbc_summary <- meta %>%
  dplyr::distinct(new_id_4, new_group, Allcell_TNBCtype) %>%  # one row per sample
  dplyr::group_by(new_group, Allcell_TNBCtype) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(new_group) %>%
  dplyr::mutate(freq = n / sum(n)) %>%  # fraction of samples per new_group
  ungroup()

# Define colors for TNBC subtypes
tnbc_colors <- c("BLIA"="skyblue","BLIS"="tomato",
                 "MES"="seagreen","LAR"="orange")

# Make stacked barplot
p2<-ggplot(tnbc_summary, aes(x = new_group, y = freq, fill = Allcell_TNBCtype)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = tnbc_colors, name = "TNBC subtype") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  labs(x = "Treatment Group", y = "Fraction of Samples") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
p2

# Chi-square test between pre-pCR and pre-RD
tab <- meta %>%
  dplyr::filter(new_group %in% c("pre-pCR","pre-RD")) %>%
  dplyr::distinct(new_id_4, new_group, Allcell_TNBCtype) %>% 
  dplyr::count(new_group, Allcell_TNBCtype) %>%
  tidyr::pivot_wider(names_from = Allcell_TNBCtype, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("new_group")

chi_res <- chisq.test(tab)
p_value <- chi_res$p.value

fisher_res <- fisher.test(tab)
p_value <- fisher_res$p.value #since not all cells >5, use fisher test

# Add p-value on figure
p2 <- p2 +
  annotate("text",
           x = 1.5, y = 1.05,
           label = paste0("p = ", signif(p_value, 3)),
           size = 5) +
  coord_cartesian(clip = "off")

ggsave("output/StackBar_sample_proportion_TNBCsubtype.pdf", plot = p2, width = 5, height = 5)

# Boxplot and sccomp comparing major cell type proportions
Idents(scRNA) <- "Celltype_Major"
Cellratio <- prop.table(table(Idents(scRNA), scRNA$new_id_4), margin = 2)#calculate cell type proportions per sample
sample.group<-read.csv("sample_group.csv",header=T)
head(sample.group)
sample.group<-sample.group[sample.group$new_id_4 %in% colnames(Cellratio),]
sample.order<-sample.group[,"new_id_4"]
Cellratio.sorted<-Cellratio[mixedsort(rownames(Cellratio)),sample.order]
Cellratio.t<-t(Cellratio.sorted)
cellratio.list<-data.frame(sample.group[,"new_group"],Cellratio.t[,1:length(colnames(Cellratio.t))])
cellratio.list<-cellratio.list[,c(1,3,4)]
names(cellratio.list)<-c("group","cell","ratio")

cellratio.list.sub<-cellratio.list[cellratio.list$group=="pre-pCR"|cellratio.list$group=="pre-RD",]
cellratio.list.sub$group<-factor(x=cellratio.list.sub$group,levels = c("pre-pCR","pre-RD"))
cellratio.list.sub %>% sample_n_by(cell, group, size = 1)#整理数据
cellratio.list.sub$cell<-factor(cellratio.list.sub$cell,levels=c("B cell","Plasma cell","T cell","Endothelial cell","Epithelial cell","Mesenchymal cell","Myeloid cell","pDC"))

#sccomp comparing proportions across groups
scRNA.sub<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
head(scRNA.sub@meta.data)

scRNA.sub$new_group<-factor(scRNA.sub$new_group,levels=c("pre-pCR","pre-RD"))

## ISSUES WITH HARDCODED PATH 
scRNA.sub$new_id_4<-factor(scRNA.sub$new_id_4,levels=names(table(scRNA.sub$new_id_4))[table(scRNA.sub$new_id_4)>0])
local_cache <- path.expand("~/.sccomp_models/1.10.0")
dir.create(local_cache, recursive = TRUE, showWarnings = FALSE)

# 2. Force the global environment variable to point to your new local folder
Sys.setenv(SCCOMP_MODELS_DIR = local_cache)

my_valid_path <- path.expand("~/.sccomp_models/1.10.0")
dir.create(my_valid_path, recursive = TRUE, showWarnings = FALSE)

# 2. Force BOTH the environment variable AND the global option
Sys.setenv(SCCOMP_MODELS_DIR = my_valid_path)
options(sccomp_models_dir = my_valid_path)

# 3. Double check that R sees your path (this should print your home path)
print(getOption("sccomp_models_dir"))

sccomp_result = 
  scRNA.sub |>
  sccomp_estimate( 
    formula_composition = ~ new_group, 
    .sample =  new_id_4, 
    .cell_group = Celltype_Major, 
    cores = 8
  ) |> 
  sccomp_remove_outliers(cores = 8) |> # Optional, remove let p-value be more sensitive but omit patient-specific subsets
  sccomp_test()

colnames(sccomp_result)

sccomp_result
saveRDS(sccomp_result,"output/sccomp_major_result_pre-pCRvsRD.rds")



plot_ready_data <- plot_ready_data |>
  group_by(new_id_4) |>
  mutate(proportion = count / sum(count)) |>
  ungroup()


p <- ggplot(plot_ready_data, aes(x = new_group, y = proportion, fill = new_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  facet_wrap(~ Celltype_Major, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Group",
    y = "Cell Type Proportion",
    fill = "Group"
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# DOESNT WORK
#p<-sccomp_result |> 
# sccomp_boxplot(factor = "new_group",significance_threshold = 0.1)
#ggsave("output/sccomp_major_boxplot_pre.pdf",p,width=4,height=8)

# substract FDR
sccomp.fdr.df <- sccomp_result[sccomp_result$parameter == "new_grouppre-RD", 
                               c("Celltype_Major", "c_FDR")]
colnames(sccomp.fdr.df) <- c("cell", "fdr")

# THIS NEEDED TO BE CONVERTED TO A TIBBLE
cell_fdr_df <- cellratio.list.sub %>%
  # Convert the list to a tibble so dplyr can read it
  dplyr::as_tibble() %>% 
  dplyr::select(cell) %>%
  dplyr::distinct() %>%
  dplyr::left_join(sccomp.fdr.df, by = "cell") %>%
  dplyr::mutate(
    label = paste0("FDR=", signif(fdr, 2)),
    y_label = 0.98,  
    y_line_start = 0.96, 
    y_line_end = 0.90  
  )

# boxplot
p1 <- ggboxplot(
  cellratio.list.sub, x = "cell", y = "ratio",
  color = "group", palette = c("#CC0C00FF", "#5C88DAFF"),add.params = list(size = 0.5),
  add = "jitter", ylab = "Cell ratio"
) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(
    data = cell_fdr_df,
    aes(x = cell, y = y_label, label = label),
    vjust = 0,
    size = 3
  ) +
  geom_segment(
    data = cell_fdr_df,
    aes(x = as.numeric(factor(cell)) - 0.35,
        xend = as.numeric(factor(cell)) + 0.35,
        y = y_line_start,
        yend = y_line_start),
    inherit.aes = FALSE
  )
ggsave("output/boxplot_major_sccomp.pdf",p1,width=6,height=4.5)

# Fig.S2C TNBCsubtype composition across groups
scRNA<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD","post-RD"))
head(scRNA@meta.data)
colnames(scRNA@meta.data)

# 1. Extract metadata
meta <- scRNA@meta.data[, c("new_id_4", "Celltype_Major", "new_group", "Allcell_TNBCtype")]

meta$new_group <- factor(meta$new_group,
                         levels = c("pre-pCR", "pre-RD", "post-RD"))

# 3. Summarize sample counts per TNBC subtype per new_group
tnbc_summary <- meta %>%
  dplyr::distinct(new_id_4, new_group, Allcell_TNBCtype) %>%  # one row per sample
  dplyr::group_by(new_group, Allcell_TNBCtype) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(new_group) %>%
  dplyr::mutate(freq = n / sum(n)) %>%  # fraction of samples per new_group
  ungroup()

# 2. Define colors for TNBC subtypes
tnbc_colors <- c("BLIA"="skyblue","BLIS"="tomato",
                 "MES"="seagreen","LAR"="orange")

# 4. Make stacked barplot
p2<-ggplot(tnbc_summary, aes(x = new_group, y = freq, fill = Allcell_TNBCtype)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = tnbc_colors, name = "TNBC subtype") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  labs(x = "Treatment Group", y = "Fraction of Samples") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5))
p2

# --- Chi-square test between pre-pCR and pre-RD ---
tab <- meta %>%
  dplyr::filter(new_group %in% c("pre-pCR","pre-RD")) %>%
  dplyr::distinct(new_id_4, new_group, Allcell_TNBCtype) %>% 
  dplyr::count(new_group, Allcell_TNBCtype) %>%
  tidyr::pivot_wider(names_from = Allcell_TNBCtype, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("new_group")

chi_res <- chisq.test(tab)
p_value <- chi_res$p.value

fisher_res <- fisher.test(tab)
p_value <- fisher_res$p.value #since not all cells >5, use fisher test

# Add p-value on figure
p2 <- p2 +
  annotate("text",
           x = 1.5, y = 1.05,
           label = paste0("p = ", signif(p_value, 3)),
           size = 5) +
  coord_cartesian(clip = "off")
p2

ggsave("stackBar_group_proportion_subtype.pdf", plot = p2, width = 5, height = 5)
