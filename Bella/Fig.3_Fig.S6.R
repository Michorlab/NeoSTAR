################################################################################
# Figure 3 + Supplementary Figure S6:
#   Transcriptional profiling of TME cell subsets in relation to treatment response (Fig. 3) and annotation of all
#   cell subsets (Fig. S6)
################################################################################

# ── Shared palette (36 colors) ────────────────────────────────────────────────
umapColor <- c(
  '#E5D2DD','#53A85F','#F1BB72','#F3B1A0','#D6E7A3','#57C3F3','#476D87',
  '#E95C59','#E59CC4','#AB3282','#23452F','#BD956A','#8C549C','#585658',
  '#9FA3A8','#E0D4CA','#5F3D69','#C5DEBA','#58A4C3','#E4C755','#F7F398',
  '#AA9A59','#E63863','#E39A35','#C1E6F3','#6778AE','#91D0BE','#B53E2B',
  '#712820','#DCC1DD','#CCE0F5','#CCC9E6','#625D9E','#68A180','#3A6963',
  '#968175'
)

# Group colors (pCR = red, RD = blue) — used consistently across all panels
group_colors <- c("pre-pCR" = "#CC0C00FF", "pre-RD" = "#5C88DAFF")

# Shared ggplot2 theme for publication
theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text        = element_text(size = base_size, color = "black"),
      axis.title       = element_text(size = base_size + 1, color = "black"),
      plot.title       = element_text(size = base_size + 1, hjust = 0.5,
                                      face = "bold"),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size, face = "bold"),
      legend.key.size  = unit(0.35, "cm"),
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size, face = "bold"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

################################################################################
# Fig.3A — sccomp bubble/dot plot (ES vs –log10 FDR)
################################################################################

# select TME cell subsets
scRNA.origin<-readRDS("All_cell_seurat.rds") #Seurat containing all cells
names<-names(table(scRNA.origin$Celltype_Subset))
scRNA.TME<-subset(scRNA.origin,Celltype_Subset %in% 
                names[grepl("CD4T_c",names)|grepl("CD8T_c",names)|grepl("ILC_c",names)|grepl("B_c",names)|grepl("M_c",names)|grepl("F_c",names)|grepl("EC_c",names)]
)
scRNA<-scRNA.TME

## calculate transcriptional similarity across cell subsets within each major cell type
# generate pseudobulk for each cell subset
a<-PseudobulkExpression(
  scRNA,
  group.by = "Celltype_Subset"
)
#cluster within each major cell type
matrix<-a$SCT[,grepl("CD4T-c",colnames(a$SCT))]
out.dist=dist(t(matrix),method="euclidean")    
out.hclust=hclust(out.dist,method="complete") 
pdf("output/hclust/hclust.CD4T.pdf",width=5,height=7)
plot(out.hclust, hang = -1)
dev.off()
out.1<-data.frame("cluster"=out.hclust$labels[out.hclust$order])

matrix<-a$SCT[,grepl("CD8T-c",colnames(a$SCT))]
out.dist=dist(t(matrix),method="euclidean") 
out.hclust=hclust(out.dist,method="complete") 
pdf("output/hclust/hclust.CD8T.pdf",width=5,height=7)
plot(out.hclust, hang = -1)
dev.off()
out.2<-data.frame("cluster"=out.hclust$labels[out.hclust$order])

matrix<-a$SCT[,grepl("ILC-c",colnames(a$SCT))]
out.dist=dist(t(matrix),method="euclidean") 
out.hclust=hclust(out.dist,method="complete")
pdf("output/hclust/hclust.ILC.pdf",width=5,height=7)
plot(out.hclust, hang = -1)
dev.off()
out.3<-data.frame("cluster"=out.hclust$labels[out.hclust$order])

matrix<-a$SCT[,grepl("B-c",colnames(a$SCT))]
out.dist=dist(t(matrix),method="euclidean") 
out.hclust=hclust(out.dist,method="complete")
pdf("output/hclust/hclust.B.pdf",width=5,height=7)
plot(out.hclust, hang = -1)
dev.off()
out.4<-data.frame("cluster"=out.hclust$labels[out.hclust$order])

matrix<-a$SCT[,grepl("M-c",colnames(a$SCT))]
out.dist=dist(t(matrix),method="euclidean") 
out.hclust=hclust(out.dist,method="complete") 
pdf("output/hclust/hclust.M.pdf",width=5,height=7)
plot(out.hclust, hang = -1) 
dev.off()
out.5<-data.frame("cluster"=out.hclust$labels[out.hclust$order])

matrix<-a$SCT[,grepl("F-c",colnames(a$SCT))]
out.dist=dist(t(matrix),method="euclidean")
out.hclust=hclust(out.dist,method="complete")
pdf("output/hclust/hclust.F.pdf",width=5,height=7)
plot(out.hclust, hang = -1) 
dev.off()
out.6<-data.frame("cluster"=out.hclust$labels[out.hclust$order])

matrix<-a$SCT[,grepl("EC-c",colnames(a$SCT))]
out.dist=dist(t(matrix),method="euclidean") 
out.hclust=hclust(out.dist,method="complete")
pdf("output/hclust/hclust.EC.pdf",width=5,height=7)
plot(out.hclust, hang = -1) 
dev.off()
out.7<-data.frame("cluster"=out.hclust$labels[out.hclust$order])

out<-rbind(out.2,out.1,out.3,out.4,out.5,out.6,out.7)
write.csv(out,"output/hclust_bulk_celltype.csv")

## sccomp compare between groups
# scRNA.sub<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
# scRNA.sub$new_group<-factor(scRNA.sub$new_group,levels=c("pre-pCR","pre-RD"))
# scRNA.sub$new_id_4 <- droplevels(scRNA.sub$new_id_4)
# scRNA.sub$Celltype_Subset <- droplevels(scRNA.sub$Celltype_Subset)
# 
# sccomp_result = 
#   scRNA.sub |>
#   sccomp_estimate( 
#     formula_composition = ~ new_group, 
#     .sample =  new_id_4, 
#     .cell_group = Celltype_Subset, 
#     cores = 8
#   ) |> 
#   sccomp_remove_outliers(cores = 8) |> # Optional, remove let p-value be more sensitive but omit patient-specific subsets
#   sccomp_test()
# saveRDS(sccomp_result,"output/sccomp_result_preNR.vs.R.rds")

# p<-sccomp_result |> 
#   sccomp_boxplot(factor = "group_37sample",significance_threshold = 0.05)
# ggsave("sccomp_cellchat_major_TME_boxplot_2025.10_pre.pdf",p,width=4,height=6)

# draw dotplot
sccomp.result<-read.csv("output/sccomp_ES_FDR_summary.csv")
rownames(sccomp.result)<-sccomp.result$Celltype_subset
sccomp.result[,"ES_postR.vs.preR"]<--sccomp.result[,"ES_preR.vs.postR"]
sccomp.result[,"ES_postNR.vs.preNR"]<--sccomp.result[,"ES_preNR.vs.postNR"]

hclust<- read.csv("output/hclust_bulk_celltype.csv", header = TRUE)

#set exp matrix & percent matrix
exp_mat<-sccomp.result[hclust$cluster,c("ES_preR.vs.NR","ES_postR.vs.preR","ES_postNR.vs.preNR")]
colnames(exp_mat)<-c("pre-pCR vs. pre-RD","post-pCR vs. pre-pCR","post-RD vs. pre-RD")

percent_mat<--log10(sccomp.result[hclust$cluster,c("FDR_preR.vs.NR","FDR_preR.vs.postR","FDR_preNR.vs.postNR")])
colnames(percent_mat)<-c("pre-pCR vs. pre-RD","post-pCR vs. pre-pCR","post-RD vs. pre-RD")

# annotate comparison groups
cluster_anno<-  c("pre-pCR vs. pre-RD","pre-pCR vs. post-pCR","pre-RD vs. post-RD")
column_ha<- HeatmapAnnotation(
  cluster_anno = cluster_anno,
  col = list(cluster_anno = setNames(pal_startrek()(4)[1:3], unique(cluster_anno))
  ),
  na_col = "grey"
)
#get a sense
quantile(as.matrix(exp_mat), c(0.1, 0.5, 0.9, 0.99))
col_fun = circlize::colorRamp2(c(-0.5, 0, 0.7), c("navy", "white", "firebrick3"))
quantile(as.matrix(percent_mat), c(0.1, 0.5, 0.9, 0.99))
cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h,
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/5 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}
# add asterics to indicate significance
cell_fun = function(j, i, x, y, width, height, fill){
  if(percent_mat[i, j] > 1.30103){
    grid.text(sprintf("%s", "*"), 
              x, 
              y, 
              gp = gpar(fontsize = 20,
                        col="red",fontface="bold"))
  }
}
# add circles
layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= pindex(as.matrix(percent_mat), i, j)/5 * unit(5, "mm"),
              gp = gpar(fill = col_fun(pindex(as.matrix(exp_mat), i, j)), col = NA))}
# add label for circles
lgd_list = list(
  Legend( labels = c(0,0.5,1,1.5,2), title = "-Log10(FDR)",
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0 * unit(5, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.1 * unit(5, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.2 * unit(5, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.3 * unit(5, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.4 * unit(5, "mm"),
                                             gp = gpar(fill = "black")))
  ))
# draw
pdf("output/Fig3A_sccomp_bubble.pdf",width=6,height=18)
hp<-ComplexHeatmap::Heatmap(exp_mat,
                            name = "hp",
                            heatmap_legend_param = list(title = "ES"),
                            col = col_fun,
                            rect_gp = gpar(type = "none"),
                            cell_fun = cell_fun,
                            layer_fun = layer_fun,
                            row_names_gp = gpar(fontsize = 10),cluster_rows = FALSE,cluster_columns = FALSE,
                            #heatmap_width = unit(15, "npc"),
                            #width = 5,
                            #heatmap_height = unit(15, "npc")
                            #height = 15,
                            #row_km = 4,
                            border = "black",top_annotation = column_ha
)
draw(hp, annotation_legend_list = lgd_list)
dev.off()


################################################################################
#  Fig.3B — CD8+ T cell UMAP (annotated left) + density by group (right) + dotplot of annotation markers (Fig.S6A)
################################################################################

scRNA_CD8 <- readRDS("data/CD8T_seurat.rds")#this version comes from 20240429
head(scRNA_CD8@meta.data)
scRNA <- scRNA_CD8

#re-organize sample and group ID
sample.group<-read.csv("data/sample_group.csv",header=T)
group_tmp <-sample.group$new_id_4
names(group_tmp) <- sample.group$orig.ident
Idents(scRNA) <- "orig.ident"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_id_4"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

group_tmp <-sample.group$new_group
names(group_tmp) <- sample.group$new_id_4
Idents(scRNA) <- "new_id_4"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_group"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

scRNA@meta.data<-scRNA@meta.data[,c("new_id_4","nCount_RNA","nFeature_RNA","percent.mt","percent.HB","nCount_SCT","nFeature_SCT",
                          "seq_batch","new_group","Celltype_subset","Celltype_subset_short","seurat_clusters")]
saveRDS(scRNA,"data/CD8T_seurat.rds")

## Left: annotated UMAP
pB_left <- DimPlot(scRNA,
                   group.by = "Celltype_subset",
                   label    = FALSE,
                   pt.size  = 0.01,
                   raster   = FALSE) +
  scale_color_manual(values = umapColor) +
  labs(title = bquote("CD8"^"+" * " T cells (n=24,836)")) +
  theme_void(base_size = 7) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 7, face = "bold"),
    legend.text     = element_text(size = 4.5),
    legend.key.size = unit(0.25, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
ggsave("output/Fig3B.left_CD8T_UMAP.pdf", pB_left, width = 4.5, height = 3.2)

## Right: side-by-side density plots
tempName<-"CD8T"
tempObj<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
tempObj$new_group<-factor(tempObj$new_group,levels=c("pre-pCR","pre-RD"))
#Faceted Density Plot for groups
Idents(tempObj) <- "new_group"
coord <- Embeddings(object = tempObj, reduction = "umap")[, 1:2]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord <- data.frame(ID = rownames(coord), coord)

# Metadata
meta <- tempObj@meta.data %>%
  rownames_to_column("ID") %>%
  left_join(coord, by = "ID")

# Optional: downsample equally across groups
# May consider downsample equally across samples as well
minNum <- min(table(meta$new_group))
meta <- meta %>%
  group_by(new_group) %>%
  slice_sample(n = minNum) %>%
  ungroup()

# Faceted density plot (with fixed axes)
p<-ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE
  ) +
  geom_point(color = "white", size = 0.01
  ) +
  facet_wrap(~ new_group, scales = "fixed") +  # fixed = shared axes
  scale_fill_viridis(name = "Density", option = "A") +
  theme_void() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle(paste("UMAP Density by group -", tempName))
ggsave(
  filename = paste0("output/Fig3B.right_", tempName, "_faceted_density_group_umap.pdf"),
  plot = p,
  width = 16,
  height = 7
)

## Fig.S6A dotplot showing annotation markers
# set subsets order
genegroup=read.table("data/CD8T_cluster_annot.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
seurat_annotations <-genegroup$subset_short
names(seurat_annotations) <- rownames(genegroup)
scRNA$Celltype_subset_short<-factor(scRNA$Celltype_subset_short,levels=c(seurat_annotations))

# read markers
genegroup=read.table("data/CD8T_annot_markers.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
annot.markerGenes<-rownames(genegroup)
# test the expression of these markers in assays
allGenes = row.names(GetAssayData(scRNA, slot="data"))
annot.markerGenes %in% allGenes

# draw
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="Celltype_subset_short"#"seurat_clusters","Celltype_Subset"
)+#coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=1))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/Fig.S6A_CD8T_annot_markers_Dotplot.pdf", p, width = 12, height = 4)


################################################################################
#  Fig.3C — CD8+ T cell signature violin plots (Naïve / Cytokine / Chemokine)
################################################################################

# Load and score signatures
TNBC_signature <- read.csv("CD8T_function_signature.csv", header = TRUE)
# score using Addmodulescore
colnames(TNBC_signature)
geneSets.list<-as.list(TNBC_signature[,])
for (i in 1:length(geneSets.list)){
  geneSets.list[[i]]<-geneSets.list[[i]][geneSets.list[[i]]!=""]
}
scRNA <- AddModuleScore(
  object = scRNA,
  features = geneSets.list,
  ctrl = 100, #by default
  name = names(geneSets.list)
)
colnames(scRNA@meta.data)

# Violin plot comparing cell-level scores across groups
scRNA.sub<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
signature.matrix<-scRNA.sub@meta.data[,c("Naive1","Cytokine.Cytokine.receptor6","Chemokine.Chemokine.receptor7")]#CD8T
group<-scRNA.sub$new_group
cell.type<-scRNA.sub$Celltype_Subset
plots<-list()
for(i in 1:ncol(signature.matrix)){
  exp.matrix<-data.frame("expression"=signature.matrix[,i],"group"=group,"cell"=cell.type)
  plots[[i]]<-ggviolin(exp.matrix,x="group",y="expression",fill="group",
                       palette = "nejm",
                       #add="boxplot",
                       add = "mean_sd", error.plot = "crossbar",
                       add.params = list(fill="white")
  )+ylab("Score")+
    ggtitle(colnames(signature.matrix)[i])+theme(plot.title = element_text(hjust = 0.5))+
    stat_compare_means(label = "p.format",method="wilcox.test")
  #facet_wrap(~cell)+
  #ylim(c(-2,2))
  #coord_cartesian(ylim = boxplot.stats(gene.matrix)$stats[c(1, 5)]*20)#adjust y axis range
}
signature.score.plot <- wrap_plots(plots = plots, ncol=3)
ggsave("output/Fig3C_CD8T_signature_violin.pdf", signature.score.plot, width = 10, height = 6)

################################################################################
# Fig.3D — Monocle3 pseudotime UMAP (trajectory overlay)
################################################################################

### Run trajectory analysis using monocle3
# set file path
dir.create("output/monocle")
# read expression and phenotypic data
expr_matrix <- as(as.matrix(scRNA@assays$SCT@counts), 'sparseMatrix')
p_data <- scRNA@meta.data[,c("new_id_4","new_group","Celltype_subset","Celltype_subset_short")]
f_data <- data.frame(gene_short_name = row.names(scRNA),row.names = row.names(scRNA))
# convert p_data and f_data from data.frame to AnnotatedDataFrame object
cds <- new_cell_data_set(expr_matrix,
                         cell_metadata = p_data,
                         gene_metadata = f_data)
# NormalizeData+ScaleData+RunPCA+umap by monocle
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="Celltype_subset_short") + ggtitle('cds.umap')
# integrate umap coordinates from seurat
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="Celltype_subset_short") + ggtitle('int.umap')
# save umap figures
plots = list(p1,p2)
UMAP.plot <- wrap_plots(plots = plots, ncol=2)
ggsave("output/monocle/UMAP_cluster_compare.pdf", plot = UMAP.plot, width = 10, height = 5)
## Monocle3 partition
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
ggsave("output/monocle/UMAP_partition.pdf", plot = p, width = 10, height = 5)
## identify branches
cds <- learn_graph(cds,learn_graph_control = list(euclidean_distance_ratio=0.2,#默认1
                                                  geodesic_distance_ratio=0.3,#默认1/3
                                                  minimal_branch_len=20#默认10
))
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = TRUE, 
               label_branch_points =FALSE)
## order cells
# set colors
scales::viridis_pal(alpha = 1, begin = 0, end = 1, direction = 1, option = "A")(10)
c1=viridis(10, alpha = 1, begin = 0, end = 1, direction = 1, option = "A")
cds <- order_cells(cds) # manually choose root
p<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
              label_leaves = TRUE,  label_branch_points = FALSE, label_roots = FALSE)+ scale_colour_gradient2(low="#451077FF",mid="#9F2F7FFF",high="#FCFDBFFF")
ggsave("trajectory_pseudotime.pdf", plot = p, width = 5.5, height = 4)
# save data
saveRDS(cds,"output/monocle/monocle_cds.rds")

cds <- readRDS("output/monocle/monocle_cds.rds")
# Pseudotime UMAP with trajectory graph
pD <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                 label_leaves = TRUE,  label_branch_points = FALSE, label_roots = FALSE)+ scale_colour_gradient2(low="#451077FF",mid="#9F2F7FFF",high="#FCFDBFFF")
ggsave("output/Fig3D_CD8T_pseudotime_UMAP.pdf", pD, width = 4, height = 3)


################################################################################
# Fig.3E — Pseudotime density plots: Tpex & T NK-like trajectories
################################################################################

# manually select 2 trajectories

pseudotime_vec <- monocle3::pseudotime(cds)
scRNA <- AddMetaData(scRNA, metadata = pseudotime_vec,
                         col.name = "Pseudotime")
scRNA_sub_E <- subset(scRNA,
                      new_group %in% c("pre-pCR", "pre-RD"))

cds_pex <- readRDS("output/monocle/Tpex_cds.rds")
cds_nk  <- readRDS("output/monocle/Tnk_cds.rds")

make_density_plot <- function(cell_ids, traj_label) {
  df <- scRNA_sub_E@meta.data[
    intersect(cell_ids, colnames(scRNA_sub_E)),
    c("new_group", "Pseudotime")
  ] %>%
    mutate(Group = recode(new_group, !!!group_labels))
  
  ggplot(df, aes(x = Pseudotime, color = Group)) +
    geom_density(linewidth = 0.7) +
    scale_color_manual(values = unname(group_colors),
                       labels = unname(group_labels)) +
    labs(title = traj_label,
         x = "Pseudotime", y = "Density",
         color = NULL) +
    theme_pub(base_size = 7) +
    theme(legend.position = c(0.65, 0.85),
          legend.background = element_blank())
}

pE1 <- make_density_plot(colnames(cds_pex), "CD8⁺ Tpex trajectory")
pE2 <- make_density_plot(colnames(cds_nk),  "CD8⁺ T NK-like trajectory")

pE  <- pE1 + pE2 + plot_layout(ncol = 2)
ggsave("output/Fig3E_pseudotime_density.pdf", pE, width = 5.5, height = 2.8)

################################################################################
# Fig.S6B — Line plot showing marker expression along Tpex & T NK-like trajectories
################################################################################

annot.markerGenes<-c("PDCD1","MKI67","IFNG","TCF7","GNLY")
gene.matrix<-GetAssayData(scRNA, assay="SCT", slot="scale.data")[annot.markerGenes,]
for (i in 1:length(annot.markerGenes)){
  scRNA <- AddMetaData(scRNA, metadata = gene.matrix[i,], col.name = annot.markerGenes[i])
}

cell.ID.Tpex<-colnames(cds_pex)
matrix<-scRNA@meta.data[cell.ID.Tex,c("Pseudotime",annot.markerGenes)]
matrix$trajectory<-"Tpex"

cell.ID.Tnk<-colnames(cds_nk)
matrix.3<-scRNA@meta.data[cell.ID.Tnk,c("Pseudotime",annot.markerGenes)]
matrix.3$trajectory<-"Tnk"

matrix.all<-rbind(matrix,matrix.3)
table(matrix.all$trajectory)

plots = list()
for (i in 1:length(annot.markerGenes)){
  matrix.sub<-matrix.all[,c("Pseudotime",annot.markerGenes[i],"trajectory")]
  colnames(matrix.sub)[2]<-"gene"
  plots[[i]]<-ggplot(matrix.sub,aes(x = Pseudotime, y = gene, group = trajectory))+
    #geom_line(lwd = 0.8)+theme_light()+#未拟合折线图
    geom_smooth(aes(linetype=trajectory),method = "loess",colour="black")+#曲线拟合
    scale_linetype_manual(values=c("solid","dashed","dotted"))+
    labs(x = "Pseudotime",
         y = "Scaled expression",
         title = annot.markerGenes[i])+
    theme_minimal()+ theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(legend.position = 'top')
}
curve.plot <- wrap_plots(plots = plots, ncol=5)    
ggsave("output/FigS6B_Marker_exp_trajectory.pdf", plot = curve.plot, width = 15, height = 5)


################################################################################
# Panel H — Violin: MKI67 scaled expression in CD8+ T cells (scRNA)
################################################################################

scRNA_sub_H <- subset(scRNA,
                      new_group %in% c("pre-pCR", "pre-RD"))
scRNA_sub_H$new_group <- factor(scRNA_sub_H$new_group,
                                     levels = c("pre-pCR", "pre-RD"))

mki67 <- GetAssayData(scRNA_sub_H, assay = "SCT",
                      slot = "scale.data")["MKI67", ]
df_H  <- data.frame(
  Score = mki67,
  Group = recode(scRNA_sub_H$new_group, !!!group_labels)
)

pH <- ggviolin(df_H, x = "Group", y = "Score",
               fill    = "Group",
               palette = unname(group_colors),
               add     = "mean_sd", error.plot = "crossbar",
               add.params = list(fill = "white", width = 0.12, size = 0.4)) +
  stat_compare_means(method  = "wilcox.test",
                     label   = "p.format",
                     size    = 2.5,
                     label.y = max(df_H$Score, na.rm = TRUE) * 1.08) +
  labs(title = bquote("CD8"^"+" * " T cells in scRNA"),
       y = "Scaled expression", x = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_pub(base_size = 7) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylim(c(0,35))

ggsave("output/Fig3H_MKI67_violin.pdf", pH, width = 2, height = 3)


################################################################################
#  Fig.S6C — CD4+ T cell UMAP (upper left) + density by group (upper right) + dotplot of annotation markers (lower)
################################################################################

scRNA_CD4 <- readRDS("data/CD4T_seurat.rds")#this version comes from 20240412
scRNA <- scRNA_CD4
head(scRNA@meta.data)

#re-organize sample and group ID
sample.group<-read.csv("data/sample_group.csv",header=T)
group_tmp <-sample.group$new_id_4
names(group_tmp) <- sample.group$orig.ident
Idents(scRNA) <- "orig.ident"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_id_4"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

group_tmp <-sample.group$new_group
names(group_tmp) <- sample.group$new_id_4
Idents(scRNA) <- "new_id_4"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_group"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

scRNA@meta.data<-scRNA@meta.data[,c("new_id_4","nCount_RNA","nFeature_RNA","percent.mt","percent.HB","nCount_SCT","nFeature_SCT",
                                    "seq_batch","new_group","Celltype_subset","Celltype_subset_short","seurat_clusters")]
saveRDS(scRNA,"data/CD4T_seurat.rds")

## Left: annotated UMAP
pB_left <- DimPlot(scRNA,
                   group.by = "Celltype_subset",
                   label    = FALSE,
                   pt.size  = 0.01,
                   raster   = FALSE) +
  scale_color_manual(values = umapColor) +
  labs(title = bquote("CD4"^"+" * " T cells")) +
  theme_void(base_size = 7) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 7, face = "bold"),
    legend.text     = element_text(size = 4.5),
    legend.key.size = unit(0.25, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
ggsave("output/FigS6C.left_CD4T_UMAP.pdf", pB_left, width = 4.5, height = 3.2)

## Right: side-by-side density plots
tempName<-"CD4T"
tempObj<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
tempObj$new_group<-factor(tempObj$new_group,levels=c("pre-pCR","pre-RD"))
#Faceted Density Plot for groups
Idents(tempObj) <- "new_group"
coord <- Embeddings(object = tempObj, reduction = "umap")[, 1:2]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord <- data.frame(ID = rownames(coord), coord)

# Metadata
meta <- tempObj@meta.data %>%
  rownames_to_column("ID") %>%
  left_join(coord, by = "ID")

# Optional: downsample equally across groups
# May consider downsample equally across samples as well
minNum <- min(table(meta$new_group))
meta <- meta %>%
  group_by(new_group) %>%
  slice_sample(n = minNum) %>%
  ungroup()

# Faceted density plot (with fixed axes)
p<-ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE
  ) +
  geom_point(color = "white", size = 0.01
  ) +
  facet_wrap(~ new_group, scales = "fixed") +  # fixed = shared axes
  scale_fill_viridis(name = "Density", option = "A") +
  theme_void() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle(paste("UMAP Density by group -", tempName))
ggsave(
  filename = paste0("output/FigS6C.right_", tempName, "_faceted_density_group_umap.pdf"),
  plot = p,
  width = 16,
  height = 7
)

## lower: dotplot showing annotation markers
# set subsets order
genegroup=read.table("data/CD4T_cluster_annot.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
seurat_annotations <-genegroup$subset_short
names(seurat_annotations) <- rownames(genegroup)
scRNA$Celltype_subset_short<-factor(scRNA$Celltype_subset_short,levels=unique(c(seurat_annotations)))

# read markers
genegroup=read.table("data/CD4T_annot_markers.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
annot.markerGenes<-rownames(genegroup)
# test the expression of these markers in assays
allGenes = row.names(GetAssayData(scRNA, slot="data"))
annot.markerGenes %in% allGenes

# draw
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="Celltype_subset_short"#"seurat_clusters","Celltype_Subset"
)+#coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=1))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/FigS6C.lower_CD4T_annot_markers_Dotplot.pdf", p, width = 12, height = 4)


################################################################################
#  Fig.S6D — ILC cell UMAP (upper left) + density by group (upper right) + dotplot of annotation markers (lower)
################################################################################

scRNA_NK.ILC <- readRDS("data/NK_seurat.rds")#this version comes from 20240429
scRNA <- scRNA_NK.ILC
head(scRNA@meta.data)

#re-organize sample and group ID
sample.group<-read.csv("data/sample_group.csv",header=T)
group_tmp <-sample.group$new_id_4
names(group_tmp) <- sample.group$orig.ident
Idents(scRNA) <- "orig.ident"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_id_4"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

group_tmp <-sample.group$new_group
names(group_tmp) <- sample.group$new_id_4
Idents(scRNA) <- "new_id_4"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_group"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

scRNA@meta.data<-scRNA@meta.data[,c("new_id_4","nCount_RNA","nFeature_RNA","percent.mt","percent.HB","nCount_SCT","nFeature_SCT",
                                    "seq_batch","new_group","Celltype_subset","Celltype_subset_short","seurat_clusters")]
saveRDS(scRNA,"data/NK_seurat.rds")

## Left: annotated UMAP
pB_left <- DimPlot(scRNA,
                   group.by = "Celltype_subset",
                   label    = FALSE,
                   pt.size  = 0.01,
                   raster   = FALSE) +
  scale_color_manual(values = umapColor) +
  labs(title = bquote("Innate lymphoid cells")) +
  theme_void(base_size = 7) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 7, face = "bold"),
    legend.text     = element_text(size = 4.5),
    legend.key.size = unit(0.25, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
ggsave("output/FigS6D.left_NK_UMAP.pdf", pB_left, width = 4.5, height = 3.2)

## Right: side-by-side density plots
tempName<-"NK"
tempObj<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
tempObj$new_group<-factor(tempObj$new_group,levels=c("pre-pCR","pre-RD"))
#Faceted Density Plot for groups
Idents(tempObj) <- "new_group"
coord <- Embeddings(object = tempObj, reduction = "umap")[, 1:2]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord <- data.frame(ID = rownames(coord), coord)

# Metadata
meta <- tempObj@meta.data %>%
  rownames_to_column("ID") %>%
  left_join(coord, by = "ID")

# Optional: downsample equally across groups
# May consider downsample equally across samples as well
minNum <- min(table(meta$new_group))
meta <- meta %>%
  group_by(new_group) %>%
  slice_sample(n = minNum) %>%
  ungroup()

# Faceted density plot (with fixed axes)
p<-ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE
  ) +
  geom_point(color = "white", size = 0.01
  ) +
  facet_wrap(~ new_group, scales = "fixed") +  # fixed = shared axes
  scale_fill_viridis(name = "Density", option = "A") +
  theme_void() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle(paste("UMAP Density by group -", tempName))
ggsave(
  filename = paste0("output/FigS6D.right_", tempName, "_faceted_density_group_umap.pdf"),
  plot = p,
  width = 16,
  height = 7
)

## lower: dotplot showing annotation markers
# set subsets order
genegroup=read.table("data/NK_cluster_annot.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
seurat_annotations <-genegroup$subset_short
names(seurat_annotations) <- rownames(genegroup)
scRNA$Celltype_subset_short<-factor(scRNA$Celltype_subset_short,levels=unique(c(seurat_annotations)))

# read markers
genegroup=read.table("data/NK_annot_markers.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
annot.markerGenes<-rownames(genegroup)
# test the expression of these markers in assays
allGenes = row.names(GetAssayData(scRNA, slot="data"))
annot.markerGenes %in% allGenes

# draw
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="Celltype_subset_short"#"seurat_clusters","Celltype_Subset"
)+#coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=1))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/FigS6D.lower_NK_annot_markers_Dotplot.pdf", p, width = 12, height = 4)

### Fig.3I dotplot showing group diff gene expression in NK cells
allGenes = row.names(GetAssayData(scRNA, slot="data"))
annot.markerGenes<-c("DNAJB1","EGR1","HSP90AA1","HSP90AB1","HSPA1A","HSPA1B","JUN",
                     "CCL3","CCL4","CCL4L2","CCL5",#"CXCR2","CX3CR1","IFNG","TNF",
                     "ITGA1","ITGAE","RGS1"#"CD69",
                     #"KIR2DL1","KIR2DL3","KIR3DL1","KIR3DL2","LILRB1","KLRC1"
)
annot.markerGenes %in% allGenes
# dotplot
scRNA<-subset(scRNA, new_group %in% c("pre-pCR","pre-RD"))
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="new_group"#"seurat_clusters","Celltype_Subset"
)+coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/Fig3I_NK_markers_group_Dotplot.pdf", p, width = 4.5, height = 6)


### Fig.S6E heatmap showing functional signatures across TNK cell subsets
head(scRNA.origin@meta.data)
scRNA<-subset(scRNA.origin, Celltype_Major %in% c("T cell"))

# calculate signature scores
TNBC_signature_1<-read.csv("data/literature_signature/IFN_TNF_TGF_signature.csv",header=T)
TNBC_signature_2<-read.csv("data/literature_signature/M_function_signature.csv",header=T)
TNBC_signature_3<-read.csv("data/literature_signature/KEGG_metabolism.csv",header=T)
TNBC_signature_4<-read.csv("data/literature_signature/CD8T_function_signature.csv",header=T)
geneSets.list.1<-as.list(TNBC_signature_1)
geneSets.list.2<-as.list(TNBC_signature_2)
geneSets.list.3<-as.list(TNBC_signature_3[,c("Glycolysis...Gluconeogenesis","Oxidative.phosphorylation","Fatty.acid.biosynthesis")])
geneSets.list.4<-as.list(TNBC_signature_4)
geneSets.list<-c(geneSets.list.1,geneSets.list.2,geneSets.list.3,geneSets.list.4)

length(geneSets.list)
for (i in 1:length(geneSets.list)){
  geneSets.list[[i]]<-geneSets.list[[i]][geneSets.list[[i]]!=""]
}# remove emtpy if gene number is unequal
scRNA <- AddModuleScore(#check parameters if necessary
  object = scRNA,
  features = geneSets.list,
  ctrl = 100, #default
  name = names(geneSets.list)
)

# calculate mean score for each cell subset
colnames(scRNA@meta.data)
signature.idx<-c(13:15,
                 27,28,31,29,
                 #25,24,26,
                 41:43
                 )#for TNK cells
names(table(scRNA$Celltype_Subset))
c<-matrix(,nrow=28,ncol=10)
for (i in 1:10){
  a<-aggregate(scRNA@meta.data[,signature.idx[i]],by=list(scRNA$Celltype_Subset),mean)
  c[,i]<-a$x
  rownames(c)<-a$Group.1
}
colnames(c)<-colnames(scRNA@meta.data)[signature.idx]

# read clustering info for subsets
hclust<-read.csv("output/hclust_bulk_celltype.csv",header=T)
hclust$cluster[1:28]#select TNK cells
# check and correct cell subsets names
setdiff(hclust$cluster[1:28], rownames(c))  # names in hclust but not in c
setdiff(rownames(c), hclust$cluster[1:28])  # names in c but not in hclust
hclust$cluster[1:23]<-c("CD8T_c10_Tcyc_MKI67","CD8T_c3_Tpex_CXCL13","CD8T_c5_Tpex_IFNG","CD8T_c9_Tpex_CCL4",     
                        "CD8T_c1_Tem_GZMK","CD8T_c6_Tem_HSPA8","CD8T_c0_Tnaive_CCR7","CD8T_c2_Tem_RNF125",
                        "CD8T_c7_Tem_ZNF683","CD8T_c4_Tnk_PRF1","CD8T_c8_Tnk_XCL2","CD4T_c4_Tem_FYN",
                        "CD4T_c0_Tem_GPR183","CD4T_c1_Tnaive_CCR7","CD4T_c3_Th1_JUN","CD4T_c9_Tem_CCL5",
                        "CD4T_c2_Treg_FOXP3","CD4T_c6_Tisr_IFIT3","CD4T_c7_Tem_GZMB","CD4T_c10_Tcyc_MKI67",
                        "CD4T_c8_Tem_NFKB1","CD4T_c11_Tem_S100A4","CD4T_c5_Tfh_CXCL13"
                        )

matrix<-c[hclust$cluster[1:28],]#select TNK cells
head(matrix)
# draw heatmap
cc = colorRampPalette(c("navy", "white", "firebrick3"))
pdf("output/FigS6E_heatmap_TNK_signature.pdf",width=25,height=25)
pheatmap(matrix,color = cc(100),
         #annotation_row = gene_info,annotation_colors = ann_colors,
         main="heatmap",
         #fontsize = 2,#gene name size
         scale="column",
         border_color = "white",
         cluster_rows = F,cluster_cols = F,#clustering_distance_rows = 'maximum',
         show_rownames = T,show_colnames = T,
         #treeheight_row = 30,treeheight_col = 30,
         cellheight = 15,cellwidth = 15,
         #cutree_row=1,cutree_col=2,
         gaps_row = c(11,23),
         gaps_col = c(3,7),
         display_numbers = F,legend = T,
         angle_col = "45"
)
dev.off()

################################################################################
#  Fig.S6F — B cell UMAP (upper left) + density by group (upper right) + dotplot of annotation markers (lower)
################################################################################

scRNA_B <- readRDS("data/B_seurat.rds")#this version comes from 20240910
scRNA <- scRNA_B
head(scRNA@meta.data)

#re-organize sample and group ID
sample.group<-read.csv("data/sample_group.csv",header=T)
group_tmp <-sample.group$new_id_4
names(group_tmp) <- sample.group$orig.ident
Idents(scRNA) <- "orig.ident"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_id_4"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

group_tmp <-sample.group$new_group
names(group_tmp) <- sample.group$new_id_4
Idents(scRNA) <- "new_id_4"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_group"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

scRNA@meta.data<-scRNA@meta.data[,c("new_id_4","nCount_RNA","nFeature_RNA","percent.mt","percent.HB","nCount_SCT","nFeature_SCT",
                                    "seq_batch","new_group","Celltype_subset","Celltype_subset_short","seurat_clusters")]
saveRDS(scRNA,"data/B_seurat.rds")

## Left: annotated UMAP
pB_left <- DimPlot(scRNA,
                   group.by = "Celltype_subset",
                   label    = FALSE,
                   pt.size  = 0.01,
                   raster   = FALSE) +
  scale_color_manual(values = umapColor) +
  labs(title = bquote("B cells")) +
  theme_void(base_size = 7) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 7, face = "bold"),
    legend.text     = element_text(size = 4.5),
    legend.key.size = unit(0.25, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
ggsave("output/FigS6F.left_B_UMAP.pdf", pB_left, width = 4.5, height = 3.2)

## Right: side-by-side density plots
tempName<-"B"
tempObj<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
tempObj$new_group<-factor(tempObj$new_group,levels=c("pre-pCR","pre-RD"))
#Faceted Density Plot for groups
Idents(tempObj) <- "new_group"
coord <- Embeddings(object = tempObj, reduction = "umap")[, 1:2]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord <- data.frame(ID = rownames(coord), coord)

# Metadata
meta <- tempObj@meta.data %>%
  rownames_to_column("ID") %>%
  left_join(coord, by = "ID")

# Optional: downsample equally across groups
# May consider downsample equally across samples as well
minNum <- min(table(meta$new_group))
meta <- meta %>%
  group_by(new_group) %>%
  slice_sample(n = minNum) %>%
  ungroup()

# Faceted density plot (with fixed axes)
p<-ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE
  ) +
  geom_point(color = "white", size = 0.01
  ) +
  facet_wrap(~ new_group, scales = "fixed") +  # fixed = shared axes
  scale_fill_viridis(name = "Density", option = "A") +
  theme_void() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle(paste("UMAP Density by group -", tempName))
ggsave(
  filename = paste0("output/FigS6F.right_", tempName, "_faceted_density_group_umap.pdf"),
  plot = p,
  width = 16,
  height = 7
)

## lower: dotplot showing annotation markers
# set subsets order
genegroup=read.table("data/B_cluster_annot.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
seurat_annotations <-genegroup$subset_short
names(seurat_annotations) <- rownames(genegroup)
scRNA$Celltype_subset_short<-factor(scRNA$Celltype_subset_short,levels=unique(c(seurat_annotations)))

# read markers
genegroup=read.table("data/B_annot_markers.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
annot.markerGenes<-rownames(genegroup)
# test the expression of these markers in assays
allGenes = row.names(GetAssayData(scRNA, slot="data"))
annot.markerGenes %in% allGenes

# draw
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="Celltype_subset_short"#"seurat_clusters","Celltype_Subset"
)+#coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=1))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/FigS6F.lower_B_annot_markers_Dotplot.pdf", p, width = 12, height = 4)


################################################################################
#  Fig.S6G — Myeloid cell UMAP (upper left) + density by group (upper right) + dotplot of annotation markers (lower)
################################################################################
## here ##
  
scRNA_M <- readRDS("M_DC_seurat.rds")#this version comes from 20250423
scRNA <- scRNA_M
head(scRNA@meta.data)

table(scRNA@meta.data$Celltype_subset)

#re-organize sample and group ID
sample.group<-read.csv("data/sample_group.csv",header=T)
group_tmp <-sample.group$new_id_4
names(group_tmp) <- sample.group$orig.ident
Idents(scRNA) <- "orig.ident"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_id_4"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

group_tmp <-sample.group$new_group
names(group_tmp) <- sample.group$new_id_4
Idents(scRNA) <- "new_id_4"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_group"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

scRNA@meta.data<-scRNA@meta.data[,c("new_id_4","nCount_RNA","nFeature_RNA","percent.mt","percent.HB","nCount_SCT","nFeature_SCT",
                                    "seq_batch","new_group","Celltype_subset","Celltype_subset_short","seurat_clusters")]
saveRDS(scRNA,"data/M_DC_seurat.rds")

## Left: annotated UMAP
pB_left <- DimPlot(scRNA,
                   group.by = "Celltype_subset",
                   label    = FALSE,
                   pt.size  = 0.01,
                   raster   = FALSE) +
  scale_color_manual(values = umapColor) +
  labs(title = bquote("Myeloid cells")) +
  theme_void(base_size = 7) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 7, face = "bold"),
    legend.text     = element_text(size = 4.5),
    legend.key.size = unit(0.25, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
ggsave("output/FigS6G.left_M_UMAP.pdf", pB_left, width = 4.5, height = 3.2)

## Right: side-by-side density plots
tempName<-"M"
tempObj<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
tempObj$new_group<-factor(tempObj$new_group,levels=c("pre-pCR","pre-RD"))
#Faceted Density Plot for groups
Idents(tempObj) <- "new_group"
coord <- Embeddings(object = tempObj, reduction = "umap")[, 1:2]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord <- data.frame(ID = rownames(coord), coord)

# Metadata
meta <- tempObj@meta.data %>%
  rownames_to_column("ID") %>%
  left_join(coord, by = "ID")

# Optional: downsample equally across groups
# May consider downsample equally across samples as well
minNum <- min(table(meta$new_group))
meta <- meta %>%
  group_by(new_group) %>%
  slice_sample(n = minNum) %>%
  ungroup()

# Faceted density plot (with fixed axes)
p<-ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE
  ) +
  geom_point(color = "white", size = 0.01
  ) +
  facet_wrap(~ new_group, scales = "fixed") +  # fixed = shared axes
  scale_fill_viridis(name = "Density", option = "A") +
  theme_void() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle(paste("UMAP Density by group -", tempName))
ggsave(
  filename = paste0("output/FigS6G.right_", tempName, "_faceted_density_group_umap.pdf"),
  plot = p,
  width = 16,
  height = 7
)

## lower: dotplot showing annotation markers
# set subsets order
genegroup=read.table("data/M_cluster_annot.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
seurat_annotations <-genegroup$subset_short
names(seurat_annotations) <- rownames(genegroup)
scRNA$Celltype_subset_short<-factor(scRNA$Celltype_subset_short,levels=unique(c(seurat_annotations)))

# read markers
genegroup=read.table("data/M_annot_markers.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
annot.markerGenes<-rownames(genegroup)
# test the expression of these markers in assays
allGenes = row.names(GetAssayData(scRNA, slot="data"))
annot.markerGenes %in% allGenes

# draw
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="Celltype_subset_short"#"seurat_clusters","Celltype_Subset"
)+#coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=1))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/FigS6G.lower_M_annot_markers_Dotplot.pdf", p, width = 12, height = 4)


### Fig.S6H heatmap showing functional signatures across myeloid cell subsets
# calculate signature scores
TNBC_signature_1<-read.csv("data/literature_signature/IFN_TNF_TGF_signature.csv",header=T)
TNBC_signature_2<-read.csv("data/literature_signature/M_function_signature.csv",header=T)
TNBC_signature_3<-read.csv("data/literature_signature/KEGG_metabolism.csv",header=T)
#TNBC_signature_4<-read.csv("data/literature_signature/CD8T_function_signature.csv",header=T)
geneSets.list.1<-as.list(TNBC_signature_1)
geneSets.list.2<-as.list(TNBC_signature_2)
geneSets.list.3<-as.list(TNBC_signature_3[,c("Glycolysis...Gluconeogenesis","Oxidative.phosphorylation","Fatty.acid.biosynthesis")])
#geneSets.list.4<-as.list(TNBC_signature_4)
geneSets.list<-c(geneSets.list.1,geneSets.list.2,geneSets.list.3)

length(geneSets.list)
for (i in 1:length(geneSets.list)){
  geneSets.list[[i]]<-geneSets.list[[i]][geneSets.list[[i]]!=""]
}# remove emtpy if gene number is unequal
scRNA <- AddModuleScore(#check parameters if necessary
  object = scRNA,
  features = geneSets.list,
  ctrl = 100, #default
  name = names(geneSets.list)
)

# calculate mean score for each cell subset
colnames(scRNA@meta.data)
signature.idx<-c(13:15,
                 16:17,
                 18,20,23,
                 24:26
)#for myeloid cells
names(table(scRNA$Celltype_subset))
c<-matrix(,nrow=15,ncol=11)
for (i in 1:11){
  a<-aggregate(scRNA@meta.data[,signature.idx[i]],by=list(scRNA$Celltype_subset),mean)
  c[,i]<-a$x
  rownames(c)<-a$Group.1
}
colnames(c)<-colnames(scRNA@meta.data)[signature.idx]

# read clustering info for subsets
hclust<-read.csv("output/hclust_bulk_celltype.csv",header=T)
hclust$cluster[40:54]#select myeloid cells
# check and correct cell subsets names
setdiff(hclust$cluster[40:54], rownames(c))  # names in hclust but not in c
setdiff(rownames(c), hclust$cluster[40:54])  # names in c but not in hclust
hclust$cluster[48]<-c("M_c8_Macro_SIGLEC15"
)

matrix<-c[hclust$cluster[40:54],]#select myeloid cells
head(matrix)
# draw heatmap
cc = colorRampPalette(c("navy", "white", "firebrick3"))
pdf("output/FigS6H_heatmap_M_signature.pdf",width=25,height=25)
pheatmap(matrix,color = cc(100),
         #annotation_row = gene_info,annotation_colors = ann_colors,
         main="heatmap",
         #fontsize = 2,#gene name size
         scale="column",
         border_color = "white",
         cluster_rows = F,cluster_cols = F,#clustering_distance_rows = 'maximum',
         show_rownames = T,show_colnames = T,
         #treeheight_row = 30,treeheight_col = 30,
         cellheight = 15,cellwidth = 15,
         #cutree_row=1,cutree_col=2,
         #gaps_row = c(11,23),
         gaps_col = c(3,5,8),
         display_numbers = F,legend = T,
         angle_col = "45"
)
dev.off()


################################################################################
# Fig.3K — Macrophage M1/M2 signature violin plots
################################################################################

# Score M1 / M2 signatures (assumed already in meta.data as M14, M25)
scRNA_M_sub <- subset(scRNA,
                      Celltype_subset %in% c("M_c2_Macro_APOE",
                                             "M_c4_Macro_CCL4",
                                             "M_c6_Macro_SPP1",
                                             "M_c7_Macro_SELENOP",
                                             "M_c12_Macro_ISG15",
                                             "M_c11_Cycling_MKI67",
                                             "M_c8_Macro_SIGLEC15") &
                        new_group %in% c("pre-pCR", "pre-RD"))
scRNA_M_sub$new_group <- factor(scRNA_M_sub$new_group,
                                     levels = c("pre-pCR", "pre-RD"))

make_violin_K <- function(score_col, sig_name) {
  df <- data.frame(
    Score = scRNA_M_sub@meta.data[[score_col]],
    Group = recode(scRNA_M_sub$new_group, !!!group_labels)
  )
  ggviolin(df, x = "Group", y = "Score",
           fill    = "Group",
           palette = unname(group_colors),
           add     = "mean_sd", error.plot = "crossbar",
           add.params = list(fill = "white", width = 0.12, size = 0.4)) +
    stat_compare_means(method  = "wilcox.test",
                       label   = "p.format",
                       size    = 2.5,
                       label.y = max(df$Score, na.rm = TRUE) * 1.08) +
    labs(title = sig_name, y = "Score", x = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    theme_pub(base_size = 7) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}

pK1 <- make_violin_K("M14", "M1-like signature")
pK2 <- make_violin_K("M25", "M2-like signature")
pK  <- pK1 + pK2 + plot_layout(ncol = 2)
ggsave("output/Fig3K_Macrophage_M1M2_violin.pdf", pK, width = 3.2, height = 3)


################################################################################
# Fig.3L — Macrophage Fc-receptor dot plot
################################################################################

fcr_genes <- c("FCGR3A", "FCGR2B", "FCGR2A", "FCGR1A")
p<-DotPlot(scRNA_M_sub,
           features =fcr_genes, 
           cols = c("lightgrey", "red"),
           group.by="new_group"#"seurat_clusters","Celltype_Subset"
)+coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/Fig3L_Macrophage_FcR_dotplot.pdf", p, width = 4.5, height = 6) ###looks like FCGR2A has some problem!!!


################################################################################
#  Fig.S6I — Mesenchymal cell UMAP (upper left) + density by group (upper right) + dotplot of annotation markers (lower)
################################################################################

scRNA_S <- readRDS("data/CAF_seurat.rds")#this version comes from 20240910
scRNA <- scRNA_S
head(scRNA@meta.data)

#re-organize sample and group ID
sample.group<-read.csv("data/sample_group.csv",header=T)
group_tmp <-sample.group$new_id_4
names(group_tmp) <- sample.group$orig.ident
Idents(scRNA) <- "orig.ident"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_id_4"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

group_tmp <-sample.group$new_group
names(group_tmp) <- sample.group$new_id_4
Idents(scRNA) <- "new_id_4"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_group"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

scRNA@meta.data<-scRNA@meta.data[,c("new_id_4","nCount_RNA","nFeature_RNA","percent.mt","percent.HB","nCount_SCT","nFeature_SCT",
                                    "seq_batch","new_group","Celltype_subset","Celltype_subset_short","seurat_clusters")]
saveRDS(scRNA,"data/CAF_seurat.rds")

## Left: annotated UMAP
pB_left <- DimPlot(scRNA,
                   group.by = "Celltype_subset",
                   label    = FALSE,
                   pt.size  = 0.01,
                   raster   = FALSE) +
  scale_color_manual(values = umapColor) +
  labs(title = bquote("Mesenchymal cells")) +
  theme_void(base_size = 7) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 7, face = "bold"),
    legend.text     = element_text(size = 4.5),
    legend.key.size = unit(0.25, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
ggsave("output/FigS6I.left_CAF_UMAP.pdf", pB_left, width = 4.5, height = 3.2)

## Right: side-by-side density plots
tempName<-"CAF"
tempObj<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
tempObj$new_group<-factor(tempObj$new_group,levels=c("pre-pCR","pre-RD"))
#Faceted Density Plot for groups
Idents(tempObj) <- "new_group"
coord <- Embeddings(object = tempObj, reduction = "umap")[, 1:2]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord <- data.frame(ID = rownames(coord), coord)

# Metadata
meta <- tempObj@meta.data %>%
  rownames_to_column("ID") %>%
  left_join(coord, by = "ID")

# Optional: downsample equally across groups
# May consider downsample equally across samples as well
minNum <- min(table(meta$new_group))
meta <- meta %>%
  group_by(new_group) %>%
  slice_sample(n = minNum) %>%
  ungroup()

# Faceted density plot (with fixed axes)
p<-ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE
  ) +
  geom_point(color = "white", size = 0.01
  ) +
  facet_wrap(~ new_group, scales = "fixed") +  # fixed = shared axes
  scale_fill_viridis(name = "Density", option = "A") +
  theme_void() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle(paste("UMAP Density by group -", tempName))
ggsave(
  filename = paste0("output/FigS6I.right_", tempName, "_faceted_density_group_umap.pdf"),
  plot = p,
  width = 16,
  height = 7
)

## lower: dotplot showing annotation markers
# set subsets order
genegroup=read.table("data/S_cluster_annot.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
seurat_annotations <-genegroup$subset_short
names(seurat_annotations) <- rownames(genegroup)
scRNA$Celltype_subset_short<-factor(scRNA$Celltype_subset_short,levels=unique(c(seurat_annotations)))

# read markers
genegroup=read.table("data/S_annot_markers.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
annot.markerGenes<-rownames(genegroup)
# test the expression of these markers in assays
allGenes = row.names(GetAssayData(scRNA, slot="data"))
annot.markerGenes %in% allGenes

# draw
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="Celltype_subset_short"#"seurat_clusters","Celltype_Subset"
)+#coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=1))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/FigS6I.lower_CAF_annot_markers_Dotplot.pdf", p, width = 12, height = 4)


################################################################################
#  Fig.S6J — Endothelial cell UMAP (upper left) + density by group (upper right) + dotplot of annotation markers (lower)
################################################################################

scRNA_Endo <- readRDS("data/EC_seurat.rds")#this version comes from 20240917
scRNA <- scRNA_Endo
head(scRNA@meta.data)

#re-organize sample and group ID
sample.group<-read.csv("data/sample_group.csv",header=T)
group_tmp <-sample.group$new_id_4
names(group_tmp) <- sample.group$orig.ident
Idents(scRNA) <- "orig.ident"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_id_4"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

group_tmp <-sample.group$new_group
names(group_tmp) <- sample.group$new_id_4
Idents(scRNA) <- "new_id_4"
scRNA <- RenameIdents(scRNA, group_tmp)#set annot as current idents
scRNA[["new_group"]] <- Idents(object = scRNA)
head(scRNA@meta.data)

scRNA@meta.data<-scRNA@meta.data[,c("new_id_4","nCount_RNA","nFeature_RNA","percent.mt","percent.HB","nCount_SCT","nFeature_SCT",
                                    "seq_batch","new_group","Celltype_subset","Celltype_subset_short","seurat_clusters")]
saveRDS(scRNA,"data/EC_seurat.rds")

## Left: annotated UMAP
pB_left <- DimPlot(scRNA,
                   group.by = "Celltype_subset",
                   label    = FALSE,
                   pt.size  = 0.01,
                   raster   = FALSE) +
  scale_color_manual(values = umapColor) +
  labs(title = bquote("Endothelial cells")) +
  theme_void(base_size = 7) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 7, face = "bold"),
    legend.text     = element_text(size = 4.5),
    legend.key.size = unit(0.25, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))
ggsave("output/FigS6J.left_EC_UMAP.pdf", pB_left, width = 4.5, height = 3.2)

## Right: side-by-side density plots
tempName<-"EC"
tempObj<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
tempObj$new_group<-factor(tempObj$new_group,levels=c("pre-pCR","pre-RD"))
#Faceted Density Plot for groups
Idents(tempObj) <- "new_group"
coord <- Embeddings(object = tempObj, reduction = "umap")[, 1:2]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord <- data.frame(ID = rownames(coord), coord)

# Metadata
meta <- tempObj@meta.data %>%
  rownames_to_column("ID") %>%
  left_join(coord, by = "ID")

# Optional: downsample equally across groups
# May consider downsample equally across samples as well
minNum <- min(table(meta$new_group))
meta <- meta %>%
  group_by(new_group) %>%
  slice_sample(n = minNum) %>%
  ungroup()

# Faceted density plot (with fixed axes)
p<-ggplot(meta, aes(x = UMAP_1, y = UMAP_2)) +
  stat_density_2d(
    aes(fill = after_stat(density)),
    geom = "raster",
    contour = FALSE
  ) +
  geom_point(color = "white", size = 0.01
  ) +
  facet_wrap(~ new_group, scales = "fixed") +  # fixed = shared axes
  scale_fill_viridis(name = "Density", option = "A") +
  theme_void() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle(paste("UMAP Density by group -", tempName))
ggsave(
  filename = paste0("output/FigS6J.right_", tempName, "_faceted_density_group_umap.pdf"),
  plot = p,
  width = 16,
  height = 7
)

## lower: dotplot showing annotation markers
# set subsets order
genegroup=read.table("data/EC_cluster_annot.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
seurat_annotations <-genegroup$subset_short
names(seurat_annotations) <- rownames(genegroup)
scRNA$Celltype_subset_short<-factor(scRNA$Celltype_subset_short,levels=unique(c(seurat_annotations)))

# read markers
genegroup=read.table("data/EC_annot_markers.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
annot.markerGenes<-rownames(genegroup)
# test the expression of these markers in assays
allGenes = row.names(GetAssayData(scRNA, slot="data"))
annot.markerGenes %in% allGenes

# draw
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="Celltype_subset_short"#"seurat_clusters","Celltype_Subset"
)+#coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=1))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/FigS6J.lower_EC_annot_markers_Dotplot.pdf", p, width = 12, height = 4)


################################################################################
# CODE TO REPLACE 3C Violin Plots with Sample Wise Box Plots
################################################################################

library(ggrepel)
library(stringr)

# Sample-wise boxplot comparing per-sample mean scores across groups
scRNA.sub <- subset(scRNA, new_group %in% c("pre-pCR", "pre-RD"))

signature.matrix <- scRNA.sub@meta.data[, c("Naive1", "Cytokine.Cytokine.receptor6", "Chemokine.Chemokine.receptor7")]  # CD8T
group <- scRNA.sub$new_group
sample <- scRNA.sub$new_id_4

clean_title <- function(x) {
  gsub("\\.", " ", x) %>%  # replace dots with spaces
    stringr::str_wrap(width = 20)  # wrap at 20 characters
}
plots <- list()
for (i in 1:ncol(signature.matrix)) {
  
  exp.matrix <- data.frame(
    expression = signature.matrix[, i],
    group = group,
    sample = sample
  )
  
  sample.means <- exp.matrix %>%
    group_by(sample, group) %>%
    summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop")
  
  n_labels <- sample.means %>%
    group_by(group) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(label = paste0(group, "\n(n=", n, ")"))
  label_map <- setNames(n_labels$label, n_labels$group)
  
  # Separate left (pre-pCR) and right (pre-RD) for directional nudging
  sample.means <- sample.means %>%
    mutate(nudge = ifelse(group == "pre-pCR", -0.6, 0.6))
  
  plots[[i]] <- ggboxplot(
    sample.means,
    x = "group",
    y = "expression",
    fill = "group",
    palette = "nejm",
    add = "jitter",
    add.params = list(size = 1.5, alpha = 0.8)
  ) +
    geom_text_repel(
      data = sample.means,
      aes(x = group, y = expression, label = sample),
      size = 1.8,
      max.overlaps = Inf,
      box.padding = 0.2,
      point.padding = 0.2,
      segment.size = 0.25,
      segment.color = "grey60",
      segment.alpha = 0.6,
      direction = "y",
      nudge_x = sample.means$nudge,
      min.segment.length = 0,
      force = 2
    ) +
    scale_x_discrete(labels = label_map) +
    ylab("Score") +
    xlab(NULL) +
    ggtitle(colnames(signature.matrix)[i]) +
    theme_classic(base_size = 11) 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_text(size = 9, face = "bold"),
      legend.position = "none",
      plot.margin = margin(5, 40, 5, 40)  # extra left/right margin for labels
    ) +
    stat_compare_means(
      label = "p.format",
      method = "wilcox.test",
      label.x.npc = 0.65,   # push p-value to right
      label.y.npc = 0.97,   # top of plot
      size = 3.5
    )
}

signature.score.plot <- wrap_plots(plots = plots, ncol = 3)
ggsave("output/Fig3C_CD8T_boxplots.pdf", signature.score.plot, width = 10, height = 6)




