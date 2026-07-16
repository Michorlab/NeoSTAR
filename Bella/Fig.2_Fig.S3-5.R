# macrophage object in figure 3

scRNA.origin<-readRDS("All_cell_seurat.rds") #Seurat containing all cells
scRNA.epi<-readRDS("epi_seurat.rds") #Seurat containing epithelial cells
scRNA.epi.ref<-readRDS("epi_normal1500_ref500.rds") #Seurat containing epithelial cells, external normal epithelial cells (500 cells per AV, HS, BA), and 500 immune cells as reference for inferCNV
scRNA.cancer<-readRDS("cancer_seurat.rds") #Seurat containing annotated cancer cells, coming from 20250811 version

# ======== Fig.2A TROP2 expression in all cells UMAP ========
# SLOT IS DEPRECATED - SWITCH TO LAYER
scRNA<-scRNA.origin
allGenes = row.names(GetAssayData(scRNA, layer ="data"))
# define annotated genes
annot.markerGenes<-c("TACSTD2")
# check gene existence
annot.markerGenes %in% allGenes
#作图
p<-FeaturePlot(scRNA, 
               features =annot.markerGenes, 
               cols = c("grey", "red"),
               #pt.size=0.01,
               ncol=1)
ggsave("output/TROP2_featureplot_allUMAP.pdf", plot = p, width = 7, height = 6)

# ======== Fig.2B & S3A epithelial cell UMAP ========
scRNA<-scRNA.epi
#draw UMAP
umapColor <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175','#FF7F0E', # vivid orange
               '#1F77B4', # blue
               '#2CA02C', # green
               '#D62728', # red
               '#9467BD'  # purple
)#36 colors
scRNA$Celltype_Subset_short<-factor(scRNA$Celltype_Subset_short,levels=mixedsort(levels(scRNA$Celltype_Subset_short)))
scRNA$Celltype_Subset<-factor(scRNA$Celltype_Subset,levels=mixedsort(levels(scRNA$Celltype_Subset)))
scRNA$new_id_4<-factor(scRNA$new_id_4,levels=mixedsort(levels(scRNA$new_id_4)))
scRNA$new_group<-factor(scRNA$new_group,levels=c("pre-pCR","pre-RD","post-RD"))
#plot
p1<-DimPlot(scRNA, group.by = "Celltype_Subset", label = T,pt.size=0.01,raster=FALSE)+scale_color_manual(values = umapColor)
p2<-DimPlot(scRNA, group.by = "Celltype_Subset_short", label = T,pt.size=0.01,raster=FALSE)+scale_color_manual(values = umapColor)
p3<- DimPlot(scRNA, group.by = "new_id_4", label = F,pt.size=0.01,raster=FALSE)+scale_color_manual(values = umapColor)
p4<- DimPlot(scRNA, group.by = "new_group", label = F,pt.size=0.01,raster=FALSE)+scale_color_manual(values = c("#CC0C00FF","#5C88DAFF","#FFCD00FF"))
plots = list(p1,p2,p3,p4)
UMAP.plot <- wrap_plots(plots = plots, nrow=4)
ggsave("output/epi_UMAP_annotations.pdf", plot = UMAP.plot, width = 7, height = 18)

# ======== Fig.S3B dotplot of annotation for epithelial clusters ========
allGenes = row.names(GetAssayData(scRNA, layer="data"))
# read markers
genegroup=read.table("E_annot_markers.txt",sep="\t",header=T,row.names=1,check.names=F,quote="")
annot.markerGenes<-rownames(genegroup)
# check gene existence
annot.markerGenes %in% allGenes
# dotplot
p<-DotPlot(scRNA,
           features =annot.markerGenes, 
           cols = c("lightgrey", "red"),
           group.by="Celltype_Subset"#"seurat_clusters","Celltype_Subset"
)+#coord_flip()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave("output/epi_annot_markers_Dotplot.pdf", p, width = 14, height = 5)

# ======== Fig.S3C barplots for epithelial clusters ========
# barplots showing sample origin, group, and patient occupancy per cluster
Idents(scRNA) <- "Celltype_Subset"
a<-table(scRNA$new_id_4,Idents(scRNA))#cell ratio across samples
Cellratio <- prop.table(a, margin = 2)
max<-apply(Cellratio,2,max)
max<-sort(max,decreasing=F)
max<-data.frame("cluster"=names(max),"occupancy"=max)
max$cluster<-factor(max$cluster,levels=max$cluster)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2<-factor(Cellratio$Var2,levels=max$cluster)#order by occupancy
# plot
p1 = ggplot( Cellratio, aes( x = Var2, weight = Freq, fill = Var1))+
  geom_bar(position = "stack",width=1)+
  scale_fill_manual(values=umapColor)+
  theme_classic()+
  labs(x='Cluster',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x=element_text(angle = 0,hjust = 1, vjust = 1),
        axis.ticks.x = element_line())+
  scale_y_continuous(limits=c(0, 1),
                     expand = c(0, 0),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))

p2<-ggplot(max,aes( x = cluster, y=occupancy, fill=occupancy))+
  geom_bar(stat="identity",width=1)+
  #scale_fill_startrek()+
  theme_classic()+
  labs(x='Cluster',y = 'Patient occupancy')+
  coord_flip()+
  scale_fill_gradient(low = "black", high="gray")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x=element_text(angle = 0,hjust = 1, vjust = 1),
        axis.ticks.x = element_line())+
  scale_y_continuous(limits=c(0, 1),
                     expand = c(0, 0),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))

a<-table(scRNA$new_group,Idents(scRNA))#cell ratio across groups
Cellratio <- prop.table(a, margin = 2)
Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2<-factor(Cellratio$Var2,levels=max$cluster)
p3<-ggplot( Cellratio, aes( x = Var2, weight = Freq, fill = Var1))+
  geom_bar(position = "stack",width=1)+
  scale_fill_manual(values=c("#CC0C00FF","#5C88DAFF","#FFCD00FF"))+
  theme_classic()+
  labs(x='Cluster',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x=element_text(angle = 0,hjust = 1, vjust = 1),
        axis.ticks.x = element_line())+
  scale_y_continuous(limits=c(0, 1),
                     expand = c(0, 0),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))
plots = list(p1,p2,p3)
ratio.plot <- wrap_plots(plots = plots, ncol=3)
ggsave("output/stackBar_epi_cluster_cellratio.pdf", plot = ratio.plot, width = 16, height = 6)

# ======== Fig.S3D-E & S4 Identifying cancer vs normal epi based on CNV and transcriptional signatures ========
## merge external normal epi and immune cells as reference
# seurat_object.1<-readRDS("data/literature_signature/GSE180878_Li_normal_epi.rds")
# seurat_object.1$dataset<-"External epi"
# set.seed(123)
# meta<-seurat_object.1@meta.data
# cells_use <- unlist(lapply(c("AV","HS","BA"), function(x){
#   cells <- rownames(meta)[meta$Major.subtype == x]
#   sample(cells, 500)
# }))
# seurat_object.1_sub <- subset(seurat_object.1, cells = cells_use)
# 
# seurat_object.2<-subset(scRNA.origin, Celltype_Major != "Epithelial cell")
# seurat_object.2$dataset<-"NeoSTAR TME"
# cells_use<-sample(colnames(seurat_object.2),500)
# seurat_object.2_sub <- subset(seurat_object.2, cells = cells_use)
# 
# seurat_object.1_sub <- CreateSeuratObject(
#   counts = GetAssayData(seurat_object.1_sub, assay = "RNA", slot = "counts"),
#   meta.data = seurat_object.1_sub@meta.data
# )
# seurat_object.2_sub <- CreateSeuratObject(
#   counts = GetAssayData(seurat_object.2_sub, assay = "RNA", slot = "counts"),
#   meta.data = seurat_object.2_sub@meta.data
# )
# scRNA.epi_sub <- CreateSeuratObject(
#   counts = GetAssayData(scRNA.epi, assay = "RNA", slot = "counts"),
#   meta.data = scRNA.epi@meta.data
# )
# scRNA <- merge(
#   x = seurat_object.1_sub,
#   y = list(seurat_object.2_sub, scRNA.epi_sub),
#   add.cell.ids = c("NormalEpi", "Ref", "NeoSTAR")
# )
# 
# scRNA$infercnv_id <- scRNA$new_id_4
# scRNA@meta.data$infercnv_id <- as.character(scRNA@meta.data$infercnv_id)
# scRNA@meta.data$infercnv_id[scRNA$dataset == "External epi"] <- 
#   scRNA@meta.data$Major.subtype[scRNA$dataset == "External epi"]
# scRNA@meta.data$infercnv_id[scRNA$dataset == "NeoSTAR TME"] <- 
#   "normal_TME_ref"
# saveRDS(scRNA,"data/epi_normal1500_ref500.rds")

# Fig.S4B compute CNV score and transcriptional signature score
# scRNA.cnv <- scRNA.epi.ref
# ## prepare for inferCNV
# counts <- GetAssayData(scRNA.cnv, slot = 'counts',assay = "RNA")# not use SCT
# anno.df<-data.frame(scRNA.cnv$infercnv_id)
# gene_order <- "data/literature_signature/data.broadinstitute.org_Trinity_CTAT_cnv_hg38_gencode_v27.txt"
# infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
#                                     annotations_file = anno.df,
#                                     delim="\t",
#                                     gene_order_file = gene_order,
#                                     min_max_counts_per_cell = c(100, +Inf),
#                                     ref_group_names = c("normal_TME_ref"))
# # run infercnv
# infercnv_obj = infercnv::run(infercnv_obj,
#                              cutoff = 0.1,# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics, or no thresholds
#                              out_dir = "output/inferCNV_RNA_test+normal_HMM_leiden_res0.01/", 
#                              cluster_by_groups = T,
#                              #k_obs_groups = 15,#default=1
#                              HMM = FALSE, 
#                              analysis_mode='subclusters',
#                              tumor_subcluster_partition_method="leiden",#default by leiden
#                              leiden_resolution=0.01,#lower resolution brings less clusters
#                              tumor_subcluster_pval=0.05,#default by 0.1
#                              denoise = TRUE, 
#                              BayesMaxPNormal = 0,#default by 0.5, time consuming
#                              num_threads = 8#multiple threads on server
# )
# #calculate CNV score：mean of the absolute values of the CNV matrix for each cell
# scRNA.cnv@meta.data[["cnv_score"]] = vector(mode="double", length=ncol(infercnv_obj@expr.data))
# scRNA.cnv@meta.data[["cnv_score"]] =  mean(abs(infercnv_obj@expr.data[ , , drop=FALSE]-1))
# #project CNV scores back to epithelial seurat
# cnv.matrix<-scRNA.cnv$cnv_score
# scRNA <- AddMetaData(object = scRNA,                #seurat object for epithelial cells
#                      metadata = cnv.matrix,               
#                      col.name = "CNV_score_cell")

#plot scores on 2D
scRNA<-scRNA.epi.ref
scRNA.sub<-subset(scRNA,dataset %in% c("External epi","NeoSTAR epi"))
scRNA.sub@meta.data[scRNA.sub$dataset=="NeoSTAR epi","Major.subtype"]<-"NeoSTAR epi"

data<-scRNA.sub@meta.data[,c("Major.subtype","Luminal_AV_SCT1","Luminal_HS_SCT1","Basal_Myepi_SCT1","CNV_score_cell")]
scatter <- ggplot(data=data,aes(x=CNV_score_cell,y=Luminal_AV_SCT1,colour=Major.subtype,fill=Major.subtype)) + 
  geom_point(aes(fill=Major.subtype),shape=21,size=1)+#,colour="black")+
  scale_fill_manual(values= c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2"))+
  scale_colour_manual(values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2"))+
  scale_x_continuous(limits = c(0, 0.2)) +
  theme_minimal()+
  theme(
    #text=element_text(size=15,face="plain",color="black"),
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=13,face="plain",color="black"),
    legend.text= element_text(size=13,face="plain",color="black"),
    legend.title=element_text(size=12,face="plain",color="black"),
    legend.background=element_blank(),
    legend.position = c(0.9,0.15)
  )
p1<-ggMarginal(scatter,type="density",color="black",groupColour = FALSE,groupFill = TRUE)
scatter <- ggplot(data=data,aes(x=CNV_score_cell,y=Luminal_HS_SCT1,colour=Major.subtype,fill=Major.subtype)) + 
  geom_point(aes(fill=Major.subtype),shape=21,size=1)+#,colour="black")+
  scale_fill_manual(values= c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2"))+
  scale_colour_manual(values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2"))+
  scale_x_continuous(limits = c(0, 0.2))+
  theme_minimal()+
  theme(
    #text=element_text(size=15,face="plain",color="black"),
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=13,face="plain",color="black"),
    legend.text= element_text(size=13,face="plain",color="black"),
    legend.title=element_text(size=12,face="plain",color="black"),
    legend.background=element_blank(),
    legend.position = c(0.9,0.15)
  )
p2<-ggMarginal(scatter,type="density",color="black",groupColour = FALSE,groupFill = TRUE)
scatter <- ggplot(data=data,aes(x=CNV_score_cell,y=Basal_Myepi_SCT1,colour=Major.subtype,fill=Major.subtype))+ 
  geom_point(aes(fill=Major.subtype),shape=21,size=1)+#,colour="black")+
  scale_fill_manual(values= c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2"))+
  scale_colour_manual(values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2"))+
  scale_x_continuous(limits = c(0, 0.2))+
  theme_minimal()+
  theme(
    #text=element_text(size=15,face="plain",color="black"),
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=13,face="plain",color="black"),
    legend.text= element_text(size=13,face="plain",color="black"),
    legend.title=element_text(size=12,face="plain",color="black"),
    legend.background=element_blank(),
    legend.position = c(0.9,0.15)
  )
p3<-ggMarginal(scatter,type="density",color="black",groupColour = FALSE,groupFill = TRUE)
plots = list(p1,p2,p3)
scatter.plot <- wrap_plots(plots = plots, ncol=3)
ggsave("output/scatter_normal_signature_CNV_score.pdf", plot = scatter.plot, width = 18, height = 5)

# #set threshold to differentiate cancer vs. normal
# scRNA.sub$class.merge<-"NA"
# scRNA.sub@meta.data[scRNA.sub$Luminal_AV_SCT1 < 0.04 & scRNA.sub$Luminal_HS_SCT1 < 0.0 & scRNA.sub$Basal_Myepi_SCT1 < 0.15,]$class.merge<-"cancer"
# scRNA.sub@meta.data[scRNA.sub$CNV_score_cell>0.03,]$class.merge<-"cancer"
# scRNA.sub@meta.data[scRNA.sub$class.merge=="NA",]$class.merge<-"normal"

# Fig.S4E UMAP plot of CNV and transcriptional signature scores
scRNA<-scRNA.epi
# 
# #calculate scores for 3 normal breast epi subtypes
# TNBC_signature<-read.csv("data/literature_signature/normal_breast_epi_signature.csv",header=T)
# geneSets.list<-as.list(TNBC_signature)
# 
# for (i in 1:length(geneSets.list)){
#   geneSets.list[[i]]<-geneSets.list[[i]][geneSets.list[[i]]!=""]
# }# remove emtpy if gene number is unequal
# scRNA <- AddModuleScore(#check parameters if necessary
#   object = scRNA,
#   features = geneSets.list,
#   ctrl = 100, #default
#   name = names(geneSets.list)
# )
# head(scRNA@meta.data)#check meta.data

# draw UMAP for scores
colnames(scRNA@meta.data)
signature.matrix<-scRNA@meta.data[,c(14:17)]
plots<-list()
for(i in 1:ncol(signature.matrix)){
  mydata<- FetchData(scRNA,vars = c("umap_1","umap_2",colnames(signature.matrix)[i]))
  colnames(mydata)<-c("UMAP_1","UMAP_2","Score")
  if (i==1){
    mydata$Score[mydata$Score>0.06]<-0.06
    plots[[i]] <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,color = Score))+
      geom_point( size=0.1) + scale_color_viridis(option="A", limits = c(0, 0.06)) + theme_light(base_size = 15)+
      labs(title = colnames(signature.matrix)[i])+
      theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }else{
    plots[[i]] <- ggplot(mydata,aes(x = UMAP_1, y = UMAP_2, color = Score))+
      geom_point(size=0.1) + scale_color_viridis(option="A") + theme_light(base_size = 15)+
      labs(title = colnames(signature.matrix)[i])+
      theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }
}
signature.score.plot <- wrap_plots(plots = plots, ncol=2)
ggsave("output/UMAP_cnv_lineage3_score.pdf", plot = signature.score.plot, width = 10, height = 8)

# # Calculate cellcycle scores
# cc.genes
# # Get marker genes for S phase 
# s.genes <- cc.genes$s.genes 
# # Get marker genes for G2M phase 
# g2m.genes <- cc.genes$g2m.genes 
# # Compute cell cycle signature scores
# s.score   <- Seurat::AddModuleScore(scRNA, features = list(s.genes), name = "S.Score")$S.Score1
# g2m.score <- Seurat::AddModuleScore(scRNA, features = list(g2m.genes), name = "G2M.Score")$G2M.Score1
# # Store in meta.data
# scRNA$S.Score   <- s.score
# scRNA$G2M.Score <- g2m.score
# # default phase is S
# scRNA$Phase.cat <- "S"
# # if G2M is higher than S → G2M
# scRNA$Phase.cat[scRNA$G2M.Score > scRNA$S.Score] <- "G2M"
# # if both scores are negative → G1
# scRNA$Phase.cat[scRNA$S.Score < 0 & scRNA$G2M.Score < 0] <- "G1"
# # turn into factor for consistency
# #scRNA$Phase.cat <- factor(scRNA$Phase.cat, levels = c("G1", "S", "G2M"))
# # check distribution
# table(scRNA$Phase.cat)
# save data
# head(scRNA@meta.data)
# saveRDS(scRNA,"data/epi_seurat.rds")

## Fig.S4C and S4E UMAP for classifications and cellcycle phase
# draw UMAP
p1 <-DimPlot(scRNA, group.by = "Celltype_Subset_short",label = T,pt.size=0.05,raster=FALSE)+scale_color_manual(values = umapColor)
p2 <- DimPlot(scRNA, group.by = "class.merge",label = T,pt.size=0.05,raster=FALSE)+scale_color_nejm()
p3 <- DimPlot(scRNA, group.by = "Phase.cat", label = T,pt.size=0.05,raster=FALSE)+
  scale_colour_manual(values = c("G1" = "#4E79A7", "S" = "#F28E2B", "G2M" = "#E15759"))
plots = list(p1,p2,p3)
UMAP.plot <- wrap_plots(plots = plots, nrow=3)    
ggsave("output/UMAP_class_Cellcycle_Epicell.pdf", plot = UMAP.plot, width = 7, height = 15)

# Fig. S4D Draw barplot for classifications
meta <- scRNA@meta.data
meta <- meta[order(meta$new_group, meta$new_id_4), ]
ordered_samples <- unique(meta$new_id_4)
# merge new_group with sample ID
sample_group <- meta %>%
  dplyr::select(new_id_4, new_group) %>%
  dplyr::distinct()
# calculate proportions of class.merge per sample
df <- scRNA@meta.data %>%
  dplyr::count(new_id_4, class.merge) %>%
  dplyr::group_by(new_id_4) %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  left_join(sample_group, by = "new_id_4")
# fix order for sample IDs 
df$new_id_4 <- factor(df$new_id_4, levels = ordered_samples)
# draw
p <- ggplot(df, aes(x = new_id_4, y = freq, fill = class.merge)) +
  geom_bar(stat = "identity", width = 0.9) +
  facet_grid(~ new_group, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_nejm(name = "Class") +
  labs(x = "Sample ID", y = "Cell Composition (%)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave("output/stackbarplot_sample_Epi_class.pdf", p, width = 8.5, height = 4)

# Fig. S4F Draw barplot for cell cycle phase
df <- scRNA@meta.data %>%
  dplyr::count(class.merge, Phase.cat) %>%   
  dplyr::group_by(class.merge) %>%
  dplyr::mutate(freq = n / sum(n))
# draw
p<-ggplot(df, aes(x = class.merge, y = freq, fill = Phase.cat)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("G1" = "#4E79A7", "S" = "#F28E2B", "G2M" = "#E15759")) +
  labs(x = "Classification", y = "Cell Composition (%)", fill = "Phase") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())
ggsave("output/stackbarplot_phase_class.pdf", p, width = 4, height = 4)


# Fig.2C & S4G TROP2/TOP1 level and TROP2 intratumoral heterogeneity within annotated cancer cells
scRNA<-scRNA.cancer
gene.matrix<-GetAssayData(scRNA, layer="data")[c("TACSTD2"),]
a<-data.frame(gene.matrix,"sample"=scRNA$new_id_4,"group_37sample"=scRNA$new_group)
a$gene.cat<-ifelse(a$gene.matrix>0,"yes","no")
df1<-aggregate(a$gene.cat,by=list(sample=a$sample,group=a$group_37sample),function(x) length(which(x=="yes"))/length(x))
colnames(df1)<-c("sample","group","Ratio")
df2<-aggregate(a$gene.matrix,by=list(sample=a$sample,group=a$group_37sample),mean)
colnames(df2)<-c("sample","group","Mean_expression")
df<-merge(df1,df2,by="sample")
colnames(df)[2]<-"group"

df<-subset(df,group %in% c("pre-pCR", "pre-RD", "post-RD"))
my_comparisons <- list(c("pre-pCR", "pre-RD"),c("pre-RD", "post-RD"))
p1<-ggboxplot(
  df, x = "group", y = "Ratio",
  color="group", palette = "jco",add = "jitter",
  #facet.by = "cell",
  short.panel.labs = FALSE
)+scale_color_manual(values=c("#CC0C00FF","#5C88DAFF","#FFCD00FF"))+
  stat_compare_means(comparisons=my_comparisons,label = "p.format",method="wilcox.test")+
  ggtitle("TROP2 ratio")+theme(plot.title = element_text(hjust = 0.5))
p2<-ggboxplot(
  df, x = "group", y = "Mean_expression",
  color="group", palette = "jco",add = "jitter",
  #facet.by = "cell",
  short.panel.labs = FALSE
)+scale_color_manual(values=c("#CC0C00FF","#5C88DAFF","#FFCD00FF"))+
  stat_compare_means(comparisons=my_comparisons,label = "p.format",method="wilcox.test")+
  ggtitle("TROP2 mean")+theme(plot.title = element_text(hjust = 0.5))
plots = list(p1,p2)
p <- wrap_plots(plots = plots, ncol=2)    
ggsave("output/Boxplot_TROP2_cancer_group.pdf", plot = p, width = 8, height = 5)

gene.matrix<-GetAssayData(scRNA, layer="data")[c("TOP1"),]
a<-data.frame(gene.matrix,"sample"=scRNA$new_id_4,"group_37sample"=scRNA$new_group)
a$gene.cat<-ifelse(a$gene.matrix>0,"yes","no")
df1<-aggregate(a$gene.cat,by=list(sample=a$sample,group=a$group_37sample),function(x) length(which(x=="yes"))/length(x))
colnames(df1)<-c("sample","group","Ratio")
df2<-aggregate(a$gene.matrix,by=list(sample=a$sample,group=a$group_37sample),mean)
colnames(df2)<-c("sample","group","Mean_expression")
df<-merge(df1,df2,by="sample")
colnames(df)[2]<-"group"

df<-subset(df,group %in% c("pre-pCR", "pre-RD", "post-RD"))
my_comparisons <- list(c("pre-pCR", "pre-RD"),c("pre-RD", "post-RD"))
p1<-ggboxplot(
  df, x = "group", y = "Ratio",
  color="group", palette = "jco",add = "jitter",
  #facet.by = "cell",
  short.panel.labs = FALSE
)+scale_color_manual(values=c("#CC0C00FF","#5C88DAFF","#FFCD00FF"))+
  stat_compare_means(comparisons=my_comparisons,label = "p.format",method="wilcox.test")+
  ggtitle("TOP1 ratio")+theme(plot.title = element_text(hjust = 0.5))
p2<-ggboxplot(
  df, x = "group", y = "Mean_expression",
  color="group", palette = "jco",add = "jitter",
  #facet.by = "cell",
  short.panel.labs = FALSE
)+scale_color_manual(values=c("#CC0C00FF","#5C88DAFF","#FFCD00FF"))+
  stat_compare_means(comparisons=my_comparisons,label = "p.format",method="wilcox.test")+
  ggtitle("TOP1 mean")+theme(plot.title = element_text(hjust = 0.5))
plots = list(p1,p2)
p <- wrap_plots(plots = plots, ncol=2)    
ggsave("output/Boxplot_TOP1_cancer_group.pdf", plot = p, width = 8, height = 5)

# ITH test using poisson
rna_counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")
# Make a list: each element is a sample with tumor cells
samples_list <- split(colnames(scRNA), scRNA$new_id_4)
# Prepare a named list per sample
sample_data <- lapply(names(samples_list), function(sample_id) {
  cells <- samples_list[[sample_id]]
  trop2_counts <- rna_counts["TACSTD2", cells]
  total_umi <- colSums(rna_counts[, cells])
  list(
    trop2 = as.numeric(trop2_counts),
    total_umi = as.numeric(total_umi),
    n_cells = length(cells),
    group = unique(scRNA$new_group[cells]) # optional: for later group comparison
  )
})
names(sample_data) <- names(samples_list)

set.seed(123)
n_sim <- 10000
results <- lapply(names(samples_list), function(sample_id) {
  trop2 <- sample_data[[sample_id]]$trop2
  total_umi <- sample_data[[sample_id]]$total_umi
  N <- length(trop2)
  # Observed zero proportion
  obs_zero <- mean(trop2 == 0)
  # Expected lambda_i
  lambda_i <- mean(trop2) * (total_umi / mean(total_umi))
  # Simulate Poisson counts
  sim_zero <- replicate(n_sim, {
    mean(rpois(N, lambda_i) == 0)
  })
  # Empirical p-value (optional, for reference)
  p_empirical <- mean(sim_zero >= obs_zero)
  # Effect size: Δ
  delta <- obs_zero - mean(sim_zero)
  # Effect size: standardized Z*
  z_star <- (obs_zero - mean(sim_zero)) / sd(sim_zero)
  data.frame(
    sample_id = sample_id,
    n_cells = N,
    obs_zero = obs_zero,
    expected_zero = mean(sim_zero),
    p_empirical = p_empirical,
    delta = delta,
    z_star = z_star,
    group = sample_data[[sample_id]]$group
  )
})
results_df <- do.call(rbind, results)
write.csv(results_df,"output/Poisson_TROP2_ITH_cancer.csv")

#作图
results_df<-read.csv("output/Poisson_TROP2_ITH_cancer.csv")
data<-subset(results_df , group %in% c("pre-pCR", "pre-RD", "post-RD"))
my_comparisons <- list(c("pre-pCR", "pre-RD"),c("pre-RD", "post-RD") )
p3<-ggboxplot(
  data, x = "group", y = "z_star",
  color = "group", palette = "jco", add = "jitter",
  #facet.by = "exp.cat",
  short.panel.labs = FALSE,
  ylab = "Standardized zero inflation"
)+scale_color_manual(values=c("#CC0C00FF","#5C88DAFF","#FFCD00FF"))+
  stat_compare_means(comparisons=my_comparisons,label = "p.format",method="wilcox.test")+
  ggtitle("TACSTD2 heterogeneity")+theme(plot.title = element_text(hjust = 0.5))
p4<-ggboxplot(
  data, x = "group", y = "delta",
  color = "group", palette = "jco", add = "jitter",
  #facet.by = "exp.cat",
  short.panel.labs = FALSE,
  ylab = "Excess zero proportion"
)+scale_color_manual(values=c("#CC0C00FF","#5C88DAFF","#FFCD00FF"))+
  stat_compare_means(comparisons=my_comparisons,label = "p.format",method="wilcox.test")+
  ggtitle("TACSTD2 heterogeneity")+theme(plot.title = element_text(hjust = 0.5))
plots = list(p3,p4)
p <- wrap_plots(plots = plots, ncol=2)    
ggsave("output/Boxplot_TROP2_ITH_cancer_group.pdf", plot = p, width = 8, height = 5)


#calculate MP score
#generate MP using cNMF
scRNA<-scRNA.cancer
scRNA.pre<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
#load function
robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10) {
  # Select NMF programs based on the minimum overlap with other NMF programs from the same cell line
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
  names(nmf_sel) <- names(nmf_programs)
  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
  # ii) the minimum overlap with programs from another cell line
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs
  final_filter <- NULL 
  for(i in names(nmf_sel)) {
    a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]
    b <- sort(apply(a, 2, max), decreasing = T) # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
    if(inter_filter==T) b <- b[b>=inter_min] # selects programs with a maximum intersection of at least 10
    if(length(b) > 1) {
      c <- names(b[1]) 
      for(y in 2:length(b)) {
        if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  return(final_filter)                                                      
}
# get NMF matrix
sct.scale.all = GetAssayData(scRNA.pre, layer="scale.data",assay="SCT")
cell_counts <- table(scRNA.pre$new_id_4)
sample <- names(cell_counts[cell_counts > 0])
score.list<-list()
for (i in 1:length(sample)){
  score.df<-data.frame("gene"=rownames(sct.scale.all))
  scRNA.sub<-subset(scRNA,new_id_4 %in% sample[i])
  sct.scale.matrix = GetAssayData(scRNA.sub, layer ="scale.data",assay="SCT")
  sct.scale.matrix[sct.scale.matrix<0]<-0#non-negative
  sct.scale.matrix.filter<-sct.scale.matrix [rowSums(sct.scale.matrix)>0,]#exclude non-expressing genes
  for (j in 4:9){
    NMFs_per_sample = nmf(x = sct.scale.matrix.filter, rank = j, method="snmf/r", nrun = 10, seed=123456)
    score<-basis(NMFs_per_sample)#W matrix as co-efficient
    colnames(score)<-paste0(paste0(sample[i],"_rank4_9_nrun10.RDS.",j,"."),c(1:j))
    score.df.tmp<-as.data.frame(score)
    score.df.tmp$gene<-rownames(score.df.tmp)
    score.df<-as.data.frame(plyr::join(score.df ,score.df.tmp ,by="gene"))
  }
  rownames(score.df)<-score.df$gene
  score.list[[paste0(sample[i],"_rank4_9_nrun10.RDS")]]<-as.matrix(score.df[,-1])
}
# save data
saveRDS(score.list,"output/nmf_SCTscale3000_rank4_9_w_basis_24pre.rds") #generate 39*24 NMF program

score.list<-readRDS("output/nmf_SCTscale3000_rank4_9_w_basis_24pre.rds")
# get top 50 genes for each NMF program 
Genes_nmf_w_basis<-score.list
nmf_programs          <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_programs          <- lapply(nmf_programs,toupper) ## convert all genes to uppercase 
## Parameters 
intra_min_parameter <- 35 #robust within tumour (similar programs by multiple K values; at least 70% gene overlap (35 out of 50 genes)
intra_max_parameter <- 10 #robust across tumours (at least 20% similarity (by top 50 genes) with any NMF program in any of the other tumours analysed
inter_min_parameter <- 10 #non-redundant within the tumour (within each tumour, NMF programs ranked by their similarity (gene overlap) with NMFs from other tumours and selected in decreasing order; 
#once an NMF was selected, any other NMF within the tumour that had 20% overlap (or more) with the selected NMF was removed, to avoid redundancy).
# for each sample, select robust NMF programs (i.e. observed using different ranks in the same sample), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
nmf_filter_ccle       <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter, inter_filter=T, inter_min = inter_min_parameter)  
nmf_programs          <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs          <- do.call(cbind, nmf_programs)
dim(nmf_programs) #generate 90 robust programs
write.csv(nmf_programs,"output/NMFprogram_genelist_robust90.csv")

# calculate similarity between programs
nmf_intersect         <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y)))) 
# hierarchical clustering of the similarity matrix 
nmf_intersect_hc     <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]
# Cluster selected NMF programs to generate MPs
### Parameters for clustering
Min_intersect_initial <- 11    # the minimal intersection cutoff for defining the first NMF program in a cluster
Min_intersect_cluster <- 11    # the minimal intersection cutoff for adding a new NMF to the forming cluster 
Min_group_size        <- 3     # ,default 5, the minimal group size to consider for defining the first NMF program in a cluster 

Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)
Cluster_list              <- list()   ### Every entry contains the NMFs of a chosen cluster
MP_list                   <- list()
k                         <- 1
Curr_cluster              <- c()
nmf_intersect_original    <- nmf_intersect

while (Sorted_intersection[1]>Min_group_size) {  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  ### intersection between all remaining NMFs and Genes in MP 
  Genes_MP                    <- nmf_programs[,names(Sorted_intersection[1])] # Genes in the forming MP are first chosen to be those in the first NMF. Genes_MP always has only 50 genes and evolves during the formation of the cluster
  nmf_programs                <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]  # remove selected NMF
  Intersection_with_Genes_MP  <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  NMF_history                 <- Genes_MP  # has genes in all NMFs in the current cluster, for redefining Genes_MP after adding a new NMF 
  ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.    
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {  
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene
    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run across all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- paste( strsplit(i , "[.]")[[1]][1 : which(strsplit(i , "[.]")[[1]] == "RDS")]   , collapse = "."  )
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear 
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   
      Genes_MP_temp             <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    NMF_history     <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP        <- Genes_MP_temp[1:50]
    nmf_programs    <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]  # remove selected NMF
    Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  }
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]]           <- Genes_MP
  k <- k+1
  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # Remove current chosen cluster
  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   # Sort intersection of remaining NMFs not included in any of the previous clusters
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}

####  Sort Jaccard similarity plot according to new clusters
inds_sorted <- c()
for (j in 1:length(Cluster_list)){
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))
}
inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[inds_new,inds_new]) 
# Custom color palette
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
# Fig.2L plot heatmap without sample label
p<-ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  scale_y_discrete(limits = rev) +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), #设置label
        axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
# add sample label
# convert matrix back to wide format
sim_mat <- reshape2::acast(nmf_intersect_meltI_NEW, Var1 ~ Var2, value.var = "value")
# make annotation
annotation <- data.frame(Sample = substr(rownames(sim_mat), 1, 4))
rownames(annotation) <- rownames(sim_mat)
annotation$Program <- rownames(annotation)
program_order <- rev(unique(nmf_intersect_meltI_NEW$Var2))
annotation$Program <- factor(annotation$Program, levels = program_order)
label_plot_with_legend <- ggplot(annotation, aes(x = 1, y = Program, fill = Sample)) +
  geom_tile() +
  scale_fill_manual(values = umapColor[c(1:2,4:9,11:12,14:17,19:20,23,25:26,28:29,31:33)]) +
  theme_void() +
  theme(
    legend.position = "right", 
    plot.margin = margin(t = 5, r = 0, b = 5, l = 5)
  )
legend_sample <- get_legend(label_plot_with_legend)
label_plot <- label_plot_with_legend + theme(legend.position = "none")
main_plot <- plot_grid(label_plot, p, ncol = 2, rel_widths = c(0.08, 1))
final_plot <- plot_grid(main_plot, legend_sample, ncol = 2, rel_widths = c(1, 0.12))
ggsave("output/Heatmap_Jaccard_similarity_label.pdf",final_plot,width = 12, height = 8)

# save data
MP_list <-  do.call(cbind, MP_list)
write.csv(MP_list,"output/MPgenelist_MP5_pre_unfilter.csv")
df <- data.frame(lapply(Cluster_list, function(x) {
  x <- unlist(x)
  length(x) <- max(lengths(Cluster_list))
  return(x)
}))
write.csv(df,"output/MPclusterlist_MP5_pre_unfilter.csv")

# manually remove MP
nmf_programs<-read.csv("output/NMFprogram_genelist_robust90.csv",header = T,row.names = 1)
list<-c(df$Cluster_5[!is.na(df$Cluster_5)])
nmf_programs_filter<-nmf_programs[,!(colnames(nmf_programs) %in% list)]
write.csv(nmf_programs_filter,"output/NMFprogram_genelist_robust90_removeMP5.csv")
# regenerate heatmap for publication

# Fig.2M generate MP score cell heatmap 
# TNBC_signature<-read.csv("output/MPgenelist_MP4_pre_filter.csv",header=T)
# TNBC_signature<-TNBC_signature[,-1]
# geneSets.list<-as.list(TNBC_signature)
# 
# for (i in 1:length(geneSets.list)){
#   geneSets.list[[i]]<-geneSets.list[[i]][geneSets.list[[i]]!=""]
# }# remove emtpy if gene number is unequal
# scRNA <- AddModuleScore(#check parameters if necessary
#   object = scRNA,nbin = 12,
#   features = geneSets.list,
#   ctrl = 50, #default
#   name = names(geneSets.list)
# )
# 
# # assign subtype
# df<-scRNA@meta.data[,c(22:25)]
# df$max <- apply(df, 1, which.max)
# df$MP.type<-NA
# df[df$max==1,]$MP.type<-"MP1"
# df[df$max==2,]$MP.type<-"MP2"
# df[df$max==3,]$MP.type<-"MP3"
# df[df$max==4,]$MP.type<-"MP4"
# scRNA <- AddMetaData(scRNA, metadata = df$MP.type, col.name = "MP.type")

colnames(scRNA@meta.data)
Idents(scRNA)<-"MP.type"
n<-names(table(scRNA$MP.type))
cell_id<-c()
for (i in 1:4){
  scRNA.sub<-subset(scRNA,MP.type == n[i])
  a<-scRNA.sub[,order(scRNA.sub@meta.data[,21+i],decreasing = T)[1:50]]
  cell_id<-c(cell_id,colnames(a))
}
downsampled_seurat <- scRNA[,cell_id]
pdf("output/cellHeatmap_MP4.pdf",width = 6,height=5)
DoHeatmap(
  downsampled_seurat,
  features = unlist(geneSets.list),
  group.by = "MP.type",
  group.bar = TRUE,disp.max=6,
  #group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"),
  slot = "scale.data",
  assay = "SCT")+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))
dev.off()

# get MP mean score for each sample
scRNA<-scRNA.cancer
colnames(scRNA@meta.data)
signature.matrix<-scRNA@meta.data[,c(1,9,18:21)]
df<-data.frame(matrix(ncol = 4, nrow = 33))
for (i in 1:4){
  tmp<-aggregate(signature.matrix[,i+2],by=list(signature.matrix[,1]),mean)
  rownames(df)<-tmp$Group.1
  df[,i]<-tmp$x
  colnames(df)[i]<-colnames(signature.matrix)[i+2]
}
df$sample<-rownames(df)

sample.group<-read.csv("sample_group.csv",header=T)
sample.group<-sample.group[sample.group$new_id_4 %in% rownames(df),]
sample.df<-data.frame("sample"=sample.group$new_id_4,"group"=sample.group$new_group,"survival_group"=sample.group$survival_outcome,"RD2_4group"=sample.group$X4.group_outcome)
merge.df<-merge(df,sample.df,by="sample")

# Fig.2N-P Boxplot for sample level comparison
df<-merge.df
plots<-list()
for (i in 1:4){
  df.sub<-df[,c("sample","group",colnames(df)[i+1])]
  df.sub$group<-factor(df.sub$group,levels = c("pre-pCR","pre-RD","post-RD"))#decide x order
  colnames(df.sub)[3]<-"score"
  #plot
  my_comparisons <- list(c("pre-pCR", "pre-RD"),c("pre-RD", "post-RD"))
  plots[[i]]<-ggboxplot(
    df.sub, x = "group", y = "score",
    color = "group", palette = "jco",
    add = "jitter",ylab = "Score",
    short.panel.labs = FALSE
  )+scale_color_manual(values=c("#CC0C00FF","#5C88DAFF","#FFCD00FF"))+
    theme(axis.text.x=element_text(angle = 30,hjust = 1, vjust = 1))+
    stat_compare_means(comparisons=my_comparisons,#
                       label = "p.format",method="t.test")+#p.format,p.signif
    ggtitle(colnames(df)[i+1])+theme(plot.title = element_text(hjust = 0.5))
}
signature.score.plot <- wrap_plots(plots = plots, ncol=4)
ggsave("output/MP4_sample_boxplot.pdf", plot = signature.score.plot, width = 22, height = 7)


# Fig.2O Correlation plot between MP2 and GIS per sample level
df<-read.csv("pair_scRNA_exome_metric.csv",header=T)
# set theme
common_theme <- theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank(),   # remove grid lines
    plot.title = element_text(hjust = 0.5)
  )
custom_colors <- c("pre-pCR" = "#CC0C00FF", "pre-RD" = "#5C88DAFF")
# draw
p1 <- ggplot(df, aes(x = Ploidy, y = MP1, color = group)) +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = TRUE, color = "red") +
  stat_cor(method = "spearman", color = "red",
           label.x.npc = "left", label.y.npc = "top") +
  scale_color_manual(values = custom_colors) +
  common_theme
p2 <- ggplot(df, aes(x = Ploidy, y = MP2, color = group)) +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = TRUE, color = "red") +
  stat_cor(method = "spearman", color = "red",
           label.x.npc = "left", label.y.npc = "top") +
  scale_color_manual(values = custom_colors) +
  common_theme
p3 <- ggplot(df, aes(x = dominant_cin, y = MP1, color = group)) +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = TRUE, color = "red") +
  stat_cor(method = "spearman", color = "red",
           label.x.npc = "left", label.y.npc = "top") +
  scale_color_manual(values = custom_colors) +
  common_theme
p4 <- ggplot(df, aes(x = dominant_cin, y = MP2, color = group)) +
  geom_point(size = 2) +
  geom_smooth(method = lm, se = TRUE, color = "red") +
  stat_cor(method = "spearman", color = "red",
           label.x.npc = "left", label.y.npc = "top") +
  scale_color_manual(values = custom_colors) +
  common_theme
plots<-wrap_plots(p1, p2, p3, p4, ncol = 2)
ggsave("output/cor_CIN_MP_sample.pdf",plots,width=14,height=10.5)


# Fig. S5I correlation plot for MP scores
# Compute correlation matrix
scRNA.pre<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
scRNA<-scRNA.pre
a<-data.frame("MP1_IFN"=scRNA$MP_11,"MP2_cell_cycle"=scRNA$MP_22,"MP3_epi_diff"=scRNA$MP_33,"MP4_EMT_inflammatory"=scRNA$MP_44)
corr_matrix <- cor(a, method = "spearman")
color_palette <- colorRampPalette(c("blue", "white", "red"))(200)
# draw heatmap
pdf("output/cor_MP_spearman.pdf",width=6,height=5)
pheatmap(corr_matrix,
         color = color_palette,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 10,
         fontsize = 13,
         main = "Clustered Correlation Heatmap of MP Scores",
         angle_col = "45",           # x轴标签倾斜45度
         legend = TRUE,
         border_color = "white",   # 去除方格线
         breaks = seq(-0.5, 0.5, length.out = 200))  # 确保0对应白色
dev.off()

# Fig. S5J correlation plot for MP scores and IRDS score
scRNA.pre<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
scRNA<-scRNA.pre
#calculate substracted score for IRDS
# Define genes
up_genes <- c("STAT1","IFI44","IFIT3","OAS1","IFIT1","ISG15","MX1")
down_genes <- c("ZNF273","CTDSPL","UBE2D3","IGF1R","FLNB","GLUD1","ABCD3")
# Make sure gene symbols match the Seurat object's naming
up_genes <- up_genes[up_genes %in% rownames(scRNA)]
down_genes <- down_genes[down_genes %in% rownames(scRNA)]
# Calculate module scores
scRNA <- AddModuleScore(scRNA, features = list(up_genes), name = "IRDS_up")
scRNA <- AddModuleScore(scRNA, features = list(down_genes), name = "IRDS_down")
# Subtract to get IRDS score
scRNA$IRDS_score <- scRNA$IRDS_up1 - scRNA$IRDS_down1

#draw correlation plot
a<-data.frame("CNV_score"=scRNA$CNV_score_cell,"MP1_IFN"=scRNA$MP_11,"MP2_cell_cycle"=scRNA$MP_22,"IRDS_score"=scRNA$IRDS_score)
p1 <- ggplot(a, aes(x = IRDS_score, y = MP1_IFN)) +
  geom_point(size = 0.3, shape = 15) +
  geom_smooth(method = lm, se = TRUE, color = "red") +
  stat_cor(method = "spearman", 
           color = "red", 
           label.x.npc = "left", 
           label.y.npc = "top") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(a, aes(x = IRDS_score, y = MP2_cell_cycle)) +
  geom_point(size = 0.3, shape = 15) +
  geom_smooth(method = lm, se = TRUE, color = "red") +
  stat_cor(method = "spearman", 
           color = "red", 
           label.x.npc = "left", 
           label.y.npc = "top") +
  theme(plot.title = element_text(hjust = 0.5))
plots <- list(p1, p2)
SCATTER.plot <- wrap_plots(plots = plots, ncol = 2)
ggsave("output/cor_IRDS_MP_spearman.pdf",SCATTER.plot,width=7,height=3.5)


# Fig. S5G boxplot comparing MP score for 4 TNBC subtypes (just pre)
scRNA.pre<-subset(scRNA,new_group %in% c("pre-pCR","pre-RD"))
scRNA<-scRNA.pre
# calculate mean MP score for each sample
scRNA$new_id_4 <- droplevels(scRNA$new_id_4)
colnames(scRNA@meta.data)
signature.matrix<-scRNA@meta.data[,c(1,9,18:21)]
df<-data.frame(matrix(ncol = 4, nrow = 24))
for (i in 1:4){
  tmp<-aggregate(signature.matrix[,i+2],by=list(signature.matrix[,1]),mean)
  rownames(df)<-tmp$Group.1
  df[,i]<-tmp$x
  colnames(df)[i]<-colnames(signature.matrix)[i+2]
}
df$sample<-rownames(df)

sample.group<-read.csv("sample_group.csv",header=T)
sample.group<-sample.group[sample.group$new_id_4 %in% rownames(df),]
sample.df<-data.frame("sample"=sample.group$new_id_4,"group"=sample.group$new_group,"subtype"=sample.group$Allcell_TNBCtype)
merge.df<-merge(df,sample.df,by="sample")

# compare sample level mean score across TNBC subtypes
plots<-list()
for (i in 1:4){
  df.sub<-merge.df[,c("sample","subtype",colnames(merge.df)[i+1])]
  df.sub$subtype<-factor(df.sub$subtype,levels = c("BLIA","BLIS","MES","LAR"))#decide x order
  colnames(df.sub)[3]<-"score"
  #plot
  my_comparisons <- list(c("BLIA","BLIS"),c("BLIA","MES"),c("BLIA","LAR"))
  plots[[i]]<-ggboxplot(
    df.sub, x = "subtype", y = "score",
    color = "subtype", palette = "jco",
    add = "jitter",ylab = "Score",
    short.panel.labs = FALSE
  )+scale_color_manual(values=c("skyblue","tomato",
                                "seagreen","orange"))+
    theme(axis.text.x=element_text(angle = 30,hjust = 1, vjust = 1))+
    stat_compare_means(comparisons=my_comparisons,#,
                       label = "p.format",method="wilcox.test")+#p.format,p.signif
    ggtitle(colnames(merge.df)[i+1])+theme(plot.title = element_text(hjust = 0.5))
}
signature.score.plot <- wrap_plots(plots = plots, ncol=4)
signature.score.plot
ggsave("output/MP4_sample_boxplot_TNBCsubtype.pdf", plot = signature.score.plot, width = 22, height = 7)



########################################
#
# Isabella's additions, June-July 2026
#
# Figure + statistics code for the TNBC scRNA-seq manuscript revision.
# Every section is self-labeled with the output file it writes.
#
# Objects assumed already loaded in the session:
#   scRNA.origin  - all cells, annotated (Celltype_Minor / Celltype_Subset, UPPERCASE "S")
#   scRNA.cancer  - epithelial/cancer cells, carries the MP metaprogram scores
#   M_DC_seurat.rds - myeloid/DC compartment, loaded below (Celltype_subset, lowercase "s")
#
# Recurring gotchas this script guards against:
#   - AnnotationDbi masks dplyr::select / dplyr::filter  -> conflicted + explicit reassignment
#   - Celltype_Subset (all-cell object) vs Celltype_subset (myeloid object) differ in case
#   - new_group has FOUR levels (pre-pCR, pre-RD, post-pCR, post-RD); pre-treatment must be
#     filtered EXPLICITLY, not just by setting factor levels (that silently makes post- cells NA)
#   - subsetting a factor keeps unused levels -> wilcox.test errors on >2 groups; droplevels() it
#
########################################


### Setup: libraries, dplyr conflict handling, output folder ###
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scales)           # hue_pal(), percent_format()
library(readxl)
library(logistf)          # Firth-penalized logistic regression (small-n safe)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggbeeswarm)
library(ggpubr)           # ggboxplot() + stat_compare_means() for the Fig 3C-style panels
library(ggrepel)

# Force dplyr verbs to win over AnnotationDbi / S4 masking.
# Without this, select()/filter() silently dispatch to the wrong method and throw
# "unable to find an inherited method" errors mid-pipeline.
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count",  "dplyr")
select <- dplyr::select
filter <- dplyr::filter

dir.create("output", showWarnings = FALSE)

setwd("/Users/isabellapabon/Documents/Ting Data")

# Shared palette used across every figure: pCR / M1 = red, RD / M2 = blue
group_colors <- c("pre-pCR" = "#CC0C00FF", "pre-RD" = "#5C88DAFF")


### CD8+ shared data objects (metadata + expression, pre-treatment only) ###
# Built once here, reused by every CD8 figure below. Pre-treatment only throughout:
# the post-treatment groups are a different question and are not part of these figures.

# CD8+ cell metadata
cd8_meta <- scRNA.origin@meta.data %>%
  filter(Celltype_Minor == "CD8+T cell",
         new_group %in% c("pre-pCR", "pre-RD")) %>%
  mutate(new_group = droplevels(factor(new_group)))   # drop post- levels or wilcox.test sees 4 groups

# Marker expression per cell. new_id_4 (patient/sample ID) is carried through so the
# patient-level views can aggregate without re-fetching.
genes <- c("IFNG", "TNF", "GZMK", "GZMB", "GZMH", "MKI67")

expr_df <- FetchData(scRNA.origin,
                     vars = c("Celltype_Minor", "Celltype_Subset",
                              "new_group", "new_id_4", genes)) %>%
  rownames_to_column("cell_id") %>%
  filter(Celltype_Minor == "CD8+T cell",
         new_group %in% c("pre-pCR", "pre-RD")) %>%
  mutate(new_group = droplevels(factor(new_group)))

# Long form, one row per cell x gene. The CD8T_ prefix is stripped so cluster labels
# fit on an axis; downstream regexes are written to tolerate its presence or absence.
expr_long <- expr_df %>%
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expression") %>%
  mutate(cluster = gsub("CD8T_", "", Celltype_Subset),
         cluster = factor(cluster, levels = sort(unique(cluster))))

# Patient-level means, one row per patient x gene.
# This is the pseudoreplication fix: tests run on n = patients, not n = cells.
expr_long_patient <- expr_df %>%
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expression") %>%
  group_by(new_id_4, new_group, gene) %>%
  summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop")

# Pooled per-cluster summary for the dot plot: mean expression, % of cells expressing,
# and a per-gene z-score across clusters (so colors are comparable down a column).
dotplot_df <- expr_long %>%
  group_by(gene, cluster) %>%
  summarise(mean_expr   = mean(expression, na.rm = TRUE),
            pct_express = mean(expression > 0, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(gene) %>%
  mutate(z_score = as.numeric(scale(mean_expr))) %>%
  ungroup()

# Per-response subset proportions + the paper's UMAP palette (feeds pie + stacked bar)
cd8_clusters <- unique(cd8_meta$Celltype_Subset)
umapColor <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175', '#FF7F0E', '#1F77B4', '#2CA02C', '#D62728', '#9467BD')
cluster_colors <- setNames(umapColor[seq_along(cd8_clusters)], sort(cd8_clusters))

cd8_prop <- cd8_meta %>%
  group_by(new_group, Celltype_Subset) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(new_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()


### CD8+ subset composition: pie chart -> output/CD8T_subset_piechart.pdf ###
# Pooled (cell-level) proportions, one pie per response group. Descriptive only --
# the statistics for composition live in the dumbbell figure further down, which
# tests per-patient proportions rather than pooled cell counts.
p_pie <- ggplot(cd8_prop, aes(x = "", y = prop, fill = Celltype_Subset)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar("y") +
  scale_fill_manual(values = cluster_colors) +
  facet_wrap(~ new_group) +
  labs(title = "CD8+ T cell subset composition by response", fill = "Cluster") +
  theme_void() +
  theme(legend.position = "right",
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("output/CD8T_subset_piechart.pdf", plot = p_pie, width = 12, height = 6)


### CD8+ subset composition: stacked bar -> output/CD8T_subset_stackedbar.pdf ###
# Same data as the pie chart, easier to read across groups. Also descriptive.
p_stacked <- ggplot(cd8_prop, aes(x = new_group, y = prop, fill = Celltype_Subset)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(title = "CD8+ T cell subset composition by response",
       x = "Response Group", y = "Proportion", fill = "Cluster") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title  = element_text(size = 12),
        plot.title  = element_text(hjust = 0.5, face = "bold", size = 13),
        panel.grid.major.y = element_line(linetype = "dashed", color = "grey85"))

ggsave("output/CD8T_subset_stackedbar.pdf", plot = p_stacked, width = 7, height = 6)


### CD8+ subset compositional shift: per-patient Wilcoxon (BH) dumbbell -> output/CD8T_compositional_shift.pdf ###
# Answers the reviewer's pseudoreplication concern for composition: the DOTS are pooled
# proportions (what the reader wants to see), but the STARS come from per-patient
# proportions, so the test has n = patients rather than n = cells.
# Result: nothing survives BH correction (all ns) -- see da_stats.

# Lineage-based cluster order, shared with the heatmap + stacked bar so the three figures
# list clusters identically. Regex tolerates the CD8T_ prefix being present or stripped.
raw_clusters  <- unique(as.character(cd8_meta$Celltype_Subset))
lineage_of    <- sub("^.*c[0-9]+_([A-Za-z0-9]+)_.*$", "\\1", raw_clusters)
lin_levels    <- c("Tnaive", "Tem", "Tpex", "Tnk", "Tcyc")
cluster_order <- raw_clusters[order(factor(lineage_of, levels = lin_levels),
                                    as.integer(sub("^.*c([0-9]+)_.*", "\\1", raw_clusters)))]

# Pooled proportions (the plotted points)
pooled <- cd8_meta %>%
  group_by(new_group, Celltype_Subset) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(new_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Per-patient proportions -> the statistics.
# complete() 0-fills subsets a patient lacks entirely: a true 0% is data, not a missing row,
# and dropping those rows would bias the test toward patients who happen to have the subset.
per_pt <- cd8_meta %>%
  group_by(new_group, new_id_4, Celltype_Subset) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(new_group, new_id_4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(nesting(new_group, new_id_4), Celltype_Subset,
           fill = list(n = 0, prop = 0))

# Wilcoxon per subset, BH-corrected across the 11 subsets.
# tryCatch guards subsets where one group is entirely 0 (test undefined -> NA, not an error).
da_stats <- per_pt %>%
  group_by(Celltype_Subset) %>%
  summarise(
    p = tryCatch(wilcox.test(prop ~ new_group)$p.value, error = function(e) NA_real_),
    .groups = "drop") %>%
  mutate(
    padj = p.adjust(p, method = "BH"),
    star = case_when(is.na(padj) ~ "",
                     padj < 0.001 ~ "***", padj < 0.01 ~ "**",
                     padj < 0.05  ~ "*",   TRUE ~ "ns"))

print(da_stats)   # exact FDRs for the manuscript text

# Wide frame for the dumbbell: one row per subset, pCR and RD proportions side by side
plot_df <- pooled %>%
  select(new_group, Celltype_Subset, prop) %>%
  pivot_wider(names_from = new_group, values_from = prop, values_fill = 0) %>%
  left_join(da_stats, by = "Celltype_Subset") %>%
  mutate(
    Celltype_Subset = factor(Celltype_Subset, levels = rev(cluster_order)),  # rev: ggplot y-axis is bottom-up
    direction = ifelse(`pre-RD` >= `pre-pCR`, "up in RD", "down in RD"),
    xmax = pmax(`pre-pCR`, `pre-RD`))                                        # star sits past the further dot

p_shift <- ggplot(plot_df, aes(y = Celltype_Subset)) +
  # Arrow runs pCR -> RD, so its direction reads as "what happens in non-responders"
  geom_segment(aes(x = `pre-pCR`, xend = `pre-RD`,
                   yend = Celltype_Subset, color = direction),
               linewidth = 1.1,
               arrow = arrow(length = unit(0.10, "cm"), type = "closed")) +
  geom_point(aes(x = `pre-pCR`), color = "#CC0C00FF", size = 3.2) +
  geom_point(aes(x = `pre-RD`),  color = "#5C88DAFF", size = 3.2) +
  geom_text(aes(x = xmax + 0.012, label = star),
            hjust = 0, fontface = "bold", size = 4.2, na.rm = TRUE) +
  scale_color_manual(values = c("up in RD" = "#5C88DAFF",
                                "down in RD" = "#CC0C00FF"), name = "Shift") +
  scale_x_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0.02, 0.12))) +   # right headroom for the stars
  labs(
    title    = "CD8+ subset compositional shift: pre-pCR -> pre-RD",
    subtitle = paste0("Dots = pooled proportion (red pre-pCR, blue pre-RD); arrow = direction of shift\n",
                      "Stars = per-patient Wilcoxon, BH-corrected (* <.05  ** <.01  *** <.001)"),
    x = "Proportion of CD8+ compartment", y = NULL) +
  theme_bw() +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.5, face = "italic", lineheight = 1.1),
        axis.text.y   = element_text(size = 9),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave("output/CD8T_compositional_shift.pdf", plot = p_shift, width = 8.5, height = 6)


### CD8+ marker expression by cluster: boxplot -> output/CD8T_gene_expression_boxplot.pdf ###
# Exploratory, cell-level. One panel per gene, clusters on x, split by response.
# Cell-level so the boxes show real spread -- NOT the figure to quote p-values from.
p_boxcluster <- ggplot(expr_long, aes(x = cluster, y = expression, fill = new_group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  scale_fill_manual(values = group_colors) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  labs(title = "CD8+ T cell gene expression by cluster and response",
       x = "Cluster", y = "Normalized expression", fill = "Response") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text  = element_text(face = "bold"),
        plot.title  = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x))   # wrap long cluster names at underscores

ggsave("output/CD8T_gene_expression_boxplot.pdf", plot = p_boxcluster, width = 14, height = 10)


### CD8+ marker expression: dot plot -> output/CD8T_marker_dotplot.pdf ###
# Standard Seurat-style dot plot: size = % of cells expressing, color = z-scored mean.
# Pooled across response groups (the by-response version is the ComplexHeatmap below).
p_dot <- ggplot(dotplot_df, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct_express, color = z_score)) +
  scale_color_gradient2(low = "#5C88DAFF", mid = "white", high = "#CC0C00FF",
                        midpoint = 0, name = "Z-score\n(mean expr)") +
  scale_size_continuous(range = c(0.5, 7), labels = scales::percent_format(),
                        name = "% Expressed") +
  labs(title = "CD8+ T cell marker gene expression by subtype",
       x = "Marker Gene", y = "CD8+ Subtype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
        axis.text.y = element_text(size = 10),
        plot.title  = element_text(hjust = 0.5, face = "bold", size = 13),
        panel.grid.major = element_line(color = "#EFEFEF", linewidth = 0.2),
        panel.grid.minor = element_blank())

ggsave("output/CD8T_marker_dotplot.pdf", plot = p_dot, width = 10, height = 6)


### CD8+ marker expression, patient-level (Fig 6A style): beeswarm -> output/CD8T_pretreatment_expression_6Astyle.pdf ###
# The pseudoreplication-safe view of marker expression: one point per patient, labeled,
# so a reviewer can see exactly how many independent observations there are (9 pCR, 16 RD).
p_beeswarm <- ggplot(expr_long_patient, aes(x = new_group, y = expression)) +
  geom_boxplot(aes(color = new_group), fill = NA, outlier.shape = NA,
               width = 0.45, linewidth = 0.7) +
  geom_beeswarm(aes(color = new_group), size = 1.8, cex = 2.5, alpha = 0.9) +
  geom_text(aes(label = new_id_4, color = new_group), size = 2,
            hjust = -0.2, check_overlap = TRUE) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  labs(title = "CD8+ T cell gene expression by response (pre-treatment)",
       x = NULL, y = "Mean expression") +
  theme_classic() +
  theme(strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(hjust = 0.5, size = 11),
        legend.position = "none")

ggsave("output/CD8T_pretreatment_expression_6Astyle.pdf", plot = p_beeswarm, width = 8, height = 10)


### CD8+ marker dot-heatmap by response -> output/CD8T_dotheatmap_by_response.pdf ###
# Publication version of the dot plot: two panels (pCR | RD) sharing rows, so the same
# cluster can be compared across responses by eye. Built with ComplexHeatmap rather than
# ggplot so the rows can be split by lineage and annotated with cluster size.

# Per-group summary: mean expr, % expressed, within-group z-score
dotplot_grouped <- expr_long %>%
  group_by(gene, cluster, new_group) %>%
  summarise(
    mean_expr   = mean(expression, na.rm = TRUE),
    pct_express = mean(expression > 0, na.rm = TRUE),
    n_cells     = dplyr::n(),
    .groups = "drop"
  ) %>%
  # z-score scaled WITHIN each response group. Colors therefore mean "high for a cluster
  # in this group", i.e. compare down a panel, not across the two panels.
  group_by(gene, new_group) %>%
  mutate(z_score = scale(mean_expr)[, 1]) %>%
  ungroup() %>%
  mutate(new_group = factor(new_group, levels = c("pre-pCR", "pre-RD")))

# One n_cells per cluster for the left barplot (max across group/gene; the two groups
# have different totals, and we want a single representative size per row)
cluster_sizes <- dotplot_grouped %>%
  group_by(cluster) %>%
  summarise(n_cells = max(n_cells), .groups = "drop")

# Gene order curated so categories are contiguous (cytotoxic -> cytokine -> proliferation),
# not first-appearance order
gene_levels   <- c("GZMB", "GZMH", "GZMK", "IFNG", "TNF", "MKI67")
gene_category <- c(GZMB = "Cytotoxic", GZMH = "Cytotoxic", GZMK = "Cytotoxic",
                   IFNG = "Cytokine",  TNF  = "Cytokine",  MKI67 = "Proliferation")

# Cluster rows sorted numerically by cXX, then split into lineage blocks.
# "^.*c([0-9]+)_" tolerates the CD8T_ prefix; an earlier "^c([0-9]+)_" version returned
# NAs and silently scrambled the row order.
raw_clusters_h <- unique(as.character(dotplot_grouped$cluster))
cluster_levels <- raw_clusters_h[order(as.integer(sub("^.*c([0-9]+)_.*", "\\1", raw_clusters_h)))]
lineage_vec    <- sub("^.*c[0-9]+_([A-Za-z0-9]+)_.*$", "\\1", cluster_levels)
lineage_order  <- c("Tnaive", "Tem", "Tpex", "Tnk", "Tcyc")
row_split      <- factor(lineage_vec, levels = lineage_order)

groups <- c("pre-pCR", "pre-RD")

# Aligned z-score + % matrices, one pair per response group.
# Both are re-indexed to [cluster_levels, gene_levels] so the two panels and the
# cell_fun lookups all address the same cell in the same order.
make_mats <- function(g) {
  d <- filter(dotplot_grouped, new_group == g)
  to_mat <- function(val) d %>%
    select(cluster, gene, dplyr::all_of(val)) %>%
    pivot_wider(names_from = gene, values_from = dplyr::all_of(val)) %>%
    column_to_rownames("cluster") %>% as.matrix()
  list(z = to_mat("z_score")[cluster_levels, gene_levels, drop = FALSE],
       p = to_mat("pct_express")[cluster_levels, gene_levels, drop = FALSE])
}
mats <- setNames(lapply(groups, make_mats), groups)

# Colors + dot geometry
col_fun <- colorRamp2(c(-2, 0, 2), c("#5C88DAFF", "white", "#CC0C00FF"))
cell_mm <- 7; max_r <- 2.45          # max_r ~= 0.35*cell_mm so dots never touch at 100%

category_cols <- c(Cytotoxic = "#41725B", Cytokine = "#B5651D", Proliferation = "#6B4C9A")
lineage_cols  <- c(Tnaive = "#C6C6C6", Tem = "#8FBC8F", Tpex = "#E8998D",
                   Tnk = "#9FB1D4", Tcyc = "#D9C47E")

# One dot per cell: outline rectangle + circle sized by % expressed.
# rect_gp = gpar(type = "none") in the Heatmap() call suppresses the default fill so this
# function draws everything. i/j are ORIGINAL matrix indices even under row_split.
make_cellfun <- function(pmat) function(j, i, x, y, width, height, fill) {
  grid.rect(x, y, width, height, gp = gpar(col = "grey92", fill = NA, lwd = 0.4))
  pe <- pmat[i, j]
  if (!is.na(pe) && pe > 0)
    grid.circle(x, y, r = unit(pe * max_r, "mm"),
                gp = gpar(fill = fill, col = "grey35", lwd = 0.3))
}

# Top annotation: "Mean % expressing" barplot + gene-category strip.
# Each bar = mean fraction of cells expressing that gene, averaged across clusters within
# the group. axis/name shown on the left panel only, so the two panels don't duplicate them.
make_top <- function(pmat, first = FALSE) HeatmapAnnotation(
  `Mean % expressing` = anno_barplot(colMeans(pmat, na.rm = TRUE),
                                     gp = gpar(fill = "grey55", col = NA),
                                     height = unit(1, "cm"),
                                     axis = first,
                                     axis_param = list(gp = gpar(fontsize = 6))),
  Category = gene_category[gene_levels],
  col = list(Category = category_cols),
  simple_anno_size     = unit(3.5, "mm"),
  show_annotation_name = if (first) c(TRUE, FALSE) else FALSE,  # label barplot on left panel only
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 7),
  gap = unit(1, "mm")
)

# Left annotation (first panel only): lineage color strip + N-cells barplot.
# The N-cells bars are the reader's caveat -- a tiny cluster's z-score is noisy.
n_cells_vec <- cluster_sizes$n_cells[match(cluster_levels, cluster_sizes$cluster)]
left_anno <- rowAnnotation(
  Lineage = lineage_vec,
  `N cells` = anno_barplot(n_cells_vec, gp = gpar(fill = "grey45", col = NA),
                           width = unit(1.6, "cm"),
                           axis_param = list(gp = gpar(fontsize = 6))),
  col = list(Lineage = lineage_cols),
  simple_anno_size    = unit(3.5, "mm"),
  annotation_name_gp  = gpar(fontsize = 7),
  annotation_name_rot = 90,
  gap = unit(1, "mm")
)

# Heatmap builder, called once per response group.
# `first` controls everything that should appear only once across the two panels:
# row names, left annotation, the color legend, and the barplot axis.
mk_ht <- function(g, first = FALSE) {
  z <- mats[[g]]$z
  Heatmap(
    z, col = col_fun,
    rect_gp  = gpar(type = "none"),          # cell_fun draws the cells instead
    cell_fun = make_cellfun(mats[[g]]$p),
    cluster_rows = FALSE, cluster_columns = FALSE,   # order is curated, not clustered
    row_split = row_split, row_gap = unit(2, "mm"),
    row_title_gp = gpar(fontsize = 9, fontface = "bold"), row_title_rot = 0,
    top_annotation  = make_top(mats[[g]]$p, first = first),
    left_annotation = if (first) left_anno else NULL,
    width  = unit(ncol(z) * cell_mm, "mm"),  # fixed cell size -> predictable square cells
    height = unit(nrow(z) * cell_mm, "mm"),
    column_title    = g,
    column_title_gp = gpar(fontface = "bold", fontsize = 11),
    column_names_gp = gpar(fontsize = 9, fontface = "italic"),
    row_names_side  = "left", row_names_gp = gpar(fontsize = 9),
    show_row_names  = first,
    # distinct internal name per panel or ComplexHeatmap merges the two legends
    name = if (first) "Z-score\n(mean expr)" else paste0("z_", g),
    show_heatmap_legend = first,
    border = TRUE
  )
}
ht_pcr <- mk_ht("pre-pCR", first = TRUE)
ht_rd  <- mk_ht("pre-RD",  first = FALSE)

# Manual "% Expressed" size legend. ComplexHeatmap can't auto-generate this because the
# dot sizes come from cell_fun, so the geometry is reproduced by hand: diameter = 2 * radius.
pct_breaks <- c(0.25, 0.5, 0.75, 1.0)
lgd_size <- Legend(
  title = "% Expressed", labels = paste0(pct_breaks * 100, "%"),
  type = "points", pch = 21, size = unit(pct_breaks * max_r * 2, "mm"),
  legend_gp = gpar(fill = "grey55", col = "grey35"), background = "white",
  title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)
)

# Draw + save, with a caption spelling out what the top barplot means
pdf("output/CD8T_dotheatmap_by_response.pdf", width = 13, height = 7)
draw(ht_pcr + ht_rd,
     ht_gap = unit(6, "mm"),
     column_title = "CD8+ T cell marker expression by subtype and response",
     column_title_gp = gpar(fontface = "bold", fontsize = 13),
     annotation_legend_list = list(lgd_size),
     merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right")
grid.text(
  paste('Top bars ("Mean % expressing") = mean fraction of cells expressing each gene,',
        "averaged across clusters within each response group."),
  x = unit(0.5, "npc"), y = unit(0.02, "npc"),
  gp = gpar(fontsize = 8, fontface = "italic")
)
dev.off()


### MP -> pCR: Firth model, per-MP Wilcoxon (BH), MP_22 boxplot -> output/MP22_by_response.pdf ###
# scRNA.cancer holds the epithelial/cancer cells with the NMF metaprogram scores (MP_11..MP_44).
# Everything here is PATIENT-level: MP scores are averaged within a patient first, so both the
# regression and the Wilcoxon have n = 24 patients (15 RD, 9 pCR), not n = cells.

# Clinical covariates (summary form only -- the sheet is a formatted table, not tidy data,
# so it's loaded for reference rather than joined)
clinical <- read_excel("~/Documents/Ting Data/Table S1_baseline clinical characteristics.xlsx")

# Pre-treatment cells -> one mean MP score per patient -> binary outcome
patient_df <- scRNA.cancer@meta.data %>%
  filter(new_group %in% c("pre-pCR", "pre-RD")) %>%
  group_by(new_id_4, new_group) %>%
  summarise(
    MP_11 = mean(MP_11, na.rm = TRUE),
    MP_22 = mean(MP_22, na.rm = TRUE),
    MP_33 = mean(MP_33, na.rm = TRUE),
    MP_44 = mean(MP_44, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pCR = ifelse(new_group == "pre-pCR", 1, 0))

table(patient_df$pCR)   # 15 RD / 9 pCR

# Firth-penalized logistic regression rather than plain glm(): with 24 patients, 4 predictors
# and an imbalanced outcome, standard ML can separate and blow the coefficients up to infinity.
# Firth's penalty keeps the estimates finite and the profile-likelihood CIs usable.
model_firth <- logistf(pCR ~ MP_11 + MP_22 + MP_33 + MP_44, data = patient_df)
summary(model_firth)

# Univariate check on the same patient-level values, BH-corrected across the 4 MPs.
# MP_22 (cell cycle) is the one that survives (p_adj ~ 0.02).
mp_cols <- c("MP_11", "MP_22", "MP_33", "MP_44")
pvals <- sapply(mp_cols, function(mp) wilcox.test(patient_df[[mp]] ~ patient_df$pCR)$p.value)
mp_results <- data.frame(MP = mp_cols, p_raw = pvals,
                         p_adj = p.adjust(pvals, method = "BH"))
write.csv(mp_results, "output/MP_wilcox_FDR.csv", row.names = FALSE)

# MP_22 by response. One point per patient; n printed in the subtitle so the reader
# can see the sample size without hunting for it.
p_mp22 <- ggplot(patient_df, aes(x = factor(pCR, labels = c("RD", "pCR")),
                                 y = MP_22,
                                 fill = factor(pCR, labels = c("RD", "pCR")))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +      # outliers hidden: jitter already shows every point
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("RD" = "#5C88DAFF", "pCR" = "#CC0C00FF")) +
  labs(title = "MP_22 score by treatment response",
       subtitle = paste0("RD n=", sum(patient_df$pCR == 0),
                         ", pCR n=", sum(patient_df$pCR == 1)),
       x = "Response", y = "MP_22 score (patient mean)", fill = "Response") +
  theme_classic() +
  theme(legend.position = "none")

ggsave("output/MP22_by_response.pdf", plot = p_mp22, width = 4, height = 5)


########################################
#
# MYELOID / MACROPHAGE SECTIONS
#
# Separate Seurat object (M_DC_seurat.rds) with its own metadata schema:
#   - cluster column is Celltype_subset (lowercase "s"), NOT Celltype_Subset
#   - sample ID is new_id_4, response is new_group (four levels)
#   - cluster names carry an "M_" prefix (M_c4_Macro_CCL4, not c4_Macro_CCL4)
#
# NOTE ON `mac`: it is assigned the FULL object below -- all 15 myeloid/DC clusters,
# including monocytes, cDC/pDC/mregDC, mast and neutrophils. Despite the variable name it
# is not a macrophage subset. The four-panel figure and the per-cluster grid therefore
# cover the whole myeloid compartment; only the by-sample sections at the end filter down
# to macrophages proper. If those two figures should be macrophage-only, subset `mac` to
# `macro_levels` (defined below) before building them.
#
########################################


### Macrophage M1/M2 four-panel figure -> output/Macrophage_M1_M2_UMAP.pdf ###
# Panels: M1 signature UMAP | M2 signature UMAP | annotated clusters | c4-vs-c11 boxplot.
# The boxplot is scoped to c4 and c11 because those are the two clusters the reviewer named.
scRNA_M <- readRDS("M_DC_seurat.rds")
mac <- scRNA_M        # full myeloid/DC object -- see note above

# M1-like / M2-like module scores. AddModuleScore appends a "1" to the name argument,
# so these land in the metadata as M1_score1 / M2_score1 (not M1_score).
m1_genes <- c("IL1B", "NLRP3", "CCL4", "FCGR2A")
m2_genes <- c("CD163", "C1QC", "APOE", "SELENOP", "PLTP", "SPP1", "FCGR3A")

mac <- AddModuleScore(
  object = mac,
  features = list(M1_score = m1_genes),
  ctrl = 100,
  name = "M1_score"
)

mac <- AddModuleScore(
  object = mac,
  features = list(M2_score = m2_genes),
  ctrl = 100,
  name = "M2_score"
)

scRNA_M_sub$new_group

# c4-vs-c11 module scores (feeds the boxplot).
# droplevels() is essential: subsetting a factor keeps all 15 levels, and wilcox.test
# then errors with "grouping factor must have exactly 2 levels".
score_df <- mac@meta.data %>%
  filter(Celltype_subset %in% c("M_c4_Macro_CCL4", "M_c11_Cycling_MKI67")) %>%
  select(Celltype_subset, M1_score1, M2_score1) %>%
  pivot_longer(c(M1_score1, M2_score1), names_to = "signature", values_to = "score") %>%
  mutate(signature = recode(signature, "M1_score1" = "M1-like", "M2_score1" = "M2-like"),
         cluster   = droplevels(Celltype_subset))

# Wilcoxon per signature (c4 vs c11), BH-corrected across the two tests.
# Cell-level: this compares two clusters within the same tissue, so it is not the
# pseudoreplication question the by-sample figures address.
m1_p <- wilcox.test(M1_score1 ~ Celltype_subset,
                    data = filter(mac@meta.data, Celltype_subset %in% c("M_c4_Macro_CCL4", "M_c11_Cycling_MKI67")))$p.value
m2_p <- wilcox.test(M2_score1 ~ Celltype_subset,
                    data = filter(mac@meta.data, Celltype_subset %in% c("M_c4_Macro_CCL4", "M_c11_Cycling_MKI67")))$p.value
padj_c4c11 <- p.adjust(c(m1_p, m2_p), method = "BH")   # renamed: `padj` collides with other sections
fdr_labels <- data.frame(signature = c("M1-like", "M2-like"),
                         label = paste0("FDR = ", signif(padj_c4c11, 2)))

# y position for each FDR label, just above the tallest box in that signature's group
y_pos <- score_df %>%
  group_by(signature) %>%
  summarise(y = max(score, na.rm = TRUE) * 1.08, .groups = "drop") %>%
  left_join(fdr_labels, by = "signature")

# M1 signature on the UMAP. order = TRUE draws high-scoring cells last so they aren't
# buried under grey low-scoring ones.
p1 <- FeaturePlot(mac, features = "M1_score1", order = TRUE, pt.size = 0.5) +
  scale_color_gradientn(colors = c("lightgrey", "#CC0C00FF")) +
  labs(title = "M1-like signature") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# M2 signature on the UMAP
p2 <- FeaturePlot(mac, features = "M2_score1", order = TRUE, pt.size = 0.5) +
  scale_color_gradientn(colors = c("lightgrey", "#5C88DAFF")) +
  labs(title = "M2-like signature") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Reproduce Seurat's default UMAP palette so the boxplot fills can match the UMAP exactly
clust_levels <- levels(factor(mac$Celltype_subset))
mac_cols <- setNames(hue_pal()(length(clust_levels)), clust_levels)

# Annotated clusters (the color source for the boxplot)
p3 <- DimPlot(mac, group.by = "Celltype_subset", cols = mac_cols,
              label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.5) +
  labs(title = "Macrophage clusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Read the exact color the UMAP actually drew for each cluster, rather than assuming
# hue_pal() and DimPlot agree on ordering. ggplot_build()$data[[1]]$group is indexed by
# factor level, so match() maps level position -> drawn color.
umap_layer   <- ggplot_build(p3)$data[[1]]
idx          <- match(seq_along(clust_levels), umap_layer$group)
umap_palette <- setNames(umap_layer$colour[idx], clust_levels)

# c4 vs c11 module scores, filled with the same cluster colors as the UMAP.
# x = signature, fill = cluster: the two clusters sit side by side within each signature.
p4 <- ggplot(score_df, aes(x = signature, y = score, fill = cluster)) +
  geom_boxplot(outlier.size = 0.3, alpha = 1, width = 0.6,
               position = position_dodge(0.75)) +
  geom_text(data = y_pos, aes(x = signature, y = y, label = label),
            inherit.aes = FALSE, size = 3, fontface = "italic") +
  scale_fill_manual(values = umap_palette[c("M_c4_Macro_CCL4", "M_c11_Cycling_MKI67")]) +
  labs(title = "M1 vs M2 scores: c4 vs c11", x = NULL, y = "Module score", fill = "Cluster") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "right")

pdf("output/Macrophage_M1_M2_UMAP.pdf", width = 22, height = 6)
p1 | p2 | p3 | p4
dev.off()


### Per-cluster M1 vs M2 module-score boxplots -> output/Macrophage_M1_M2_c2.pdf, output/Macrophage_M1_M2_per_cluster.pdf ###
# Follow-up to the c4-vs-c11 boxplot, which compared the two clusters the reviewer named.
# Here the comparison is *within* each cluster instead: for a given cluster, are its cells
# more M1-like or more M2-like? Each cluster gets its own panel with the same layout as the
# c4/c11 boxplot -- M1-like and M2-like on the x-axis, boxes filled with that cluster's color
# taken straight from the UMAP (p3), so a panel can be matched back to the UMAP by eye.
# Clusters are never tested against each other here.
#
# Stats: paired Wilcoxon per cluster (M1_score1 vs M2_score1 in the SAME cells), BH-corrected
# across clusters. Reads as relative polarization within a cluster, not absolute expression --
# the two module scores come from different gene sets, so their magnitudes aren't interchangeable.
#
# Outputs: c2 alone (standalone panel) and all clusters (one panel each, patchwork grid).
# Requires M1_score1 / M2_score1 from the AddModuleScore step above, and p3 from the 4-panel figure.

cl_levels <- levels(factor(mac$Celltype_subset))

# Same trick as the c4/c11 boxplot: read the exact colors the UMAP drew, so the fills match.
# Falls back to hue_pal() if p3 isn't in the environment (e.g. running this section alone).
if (exists("p3")) {
  umap_layer   <- ggplot_build(p3)$data[[1]]
  idx          <- match(seq_along(cl_levels), umap_layer$group)
  umap_palette <- setNames(umap_layer$colour[idx], cl_levels)
} else {
  umap_palette <- setNames(hue_pal()(length(cl_levels)), cl_levels)
}

# Paired Wilcoxon (M1 vs M2 in the same cells), one per cluster, BH across clusters.
# Guard on n < 3: a cluster too small to test returns NA rather than erroring out the sapply.
cl_p <- sapply(cl_levels, function(cl) {
  d <- dplyr::filter(mac@meta.data, Celltype_subset == cl)
  if (nrow(d) < 3) return(NA_real_)
  wilcox.test(d$M1_score1, d$M2_score1, paired = TRUE)$p.value
})
cl_padj <- setNames(p.adjust(cl_p, method = "BH"), cl_levels)

# One M1-vs-M2 boxplot for one cluster, filled with that cluster's UMAP color.
# Both boxes share the cluster's color (identity), and M1 vs M2 is read off the x-axis --
# same encoding as the c4/c11 panel.
plot_cluster <- function(cl) {
  df <- mac@meta.data %>%
    dplyr::filter(Celltype_subset == cl) %>%
    dplyr::select(Celltype_subset, M1_score1, M2_score1) %>%
    pivot_longer(c(M1_score1, M2_score1), names_to = "signature", values_to = "score") %>%
    mutate(signature = recode(signature, M1_score1 = "M1-like", M2_score1 = "M2-like"),
           cluster   = droplevels(factor(Celltype_subset, levels = cl_levels)))
  
  lab <- data.frame(x = 1.5,                                   # centered between the two boxes
                    y = max(df$score, na.rm = TRUE) * 1.08,
                    label = paste0("FDR = ", signif(cl_padj[[cl]], 2)))
  
  ggplot(df, aes(x = signature, y = score, fill = cluster)) +
    geom_boxplot(outlier.size = 0.3, alpha = 1, width = 0.6) +
    geom_text(data = lab, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 3, fontface = "italic") +
    scale_fill_manual(values = umap_palette[cl]) +
    labs(title = cl, x = NULL, y = "Module score") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          legend.position = "none")   # cluster identity is already in the panel title
}

# c2 alone
c2_name <- grep("^M_c2_", cl_levels, value = TRUE)   # hardcode if grep is ever ambiguous
p_c2 <- plot_cluster(c2_name)

pdf("output/Macrophage_M1_M2_c2.pdf", width = 4, height = 5)
print(p_c2); dev.off()

# All clusters, one panel each.
# Y-axes are free per panel, so spread is NOT comparable across panels -- add a shared
# coord_cartesian(ylim = ...) inside plot_cluster if that matters (and pin the label y).
plots    <- lapply(cl_levels, plot_cluster)
n_col    <- 5
combined <- wrap_plots(plots, ncol = n_col)

pdf("output/Macrophage_M1_M2_per_cluster.pdf",
    width  = 3.0 * n_col,
    height = 3.6 * ceiling(length(plots) / n_col))
print(combined); dev.off()


########################################
#   1. Macrophage_M1_M2_by_sample.pdf         - THE deliverable. Same layout as the Fig 3C CD8
#                                               signature boxplots: per-sample means, labeled
#                                               points, Wilcoxon. This is what "overall ... by
#                                               sample, comparing pCR to RD" asks for.
#   2. Macrophage_M1_M2_dotheatmap_by_sample.pdf - supplementary, per-gene detail. Shows that
#                                               samples vary; does not itself carry the test.
#
# Shared setup for both is defined once immediately below.
#
########################################


### Shared: macrophage-only, pre-treatment-only cell set + sample-level scores ###

sample_col      <- "new_id_4"     # sample / patient ID
response_col    <- "new_group"    # FOUR levels; pre-treatment filtered explicitly below
min_cells       <- 10             # min macrophages for a sample's mean to be trustworthy
include_cycling <- TRUE           # c11_Cycling_MKI67 counted as macrophage (reviewer treated it
# as one in the c4-vs-c11 ask). FALSE = strict _Macro_ labels only.
min_pct         <- 0.01           # heatmap: drop genes expressed in <1% of macrophages (e.g. ARG1)

# Macrophage clusters only: monocytes (c1, c9), cDC/pDC/mregDC (c3, c13, c14, c10),
# mast (c5) and neutrophils (c0) are all excluded, per "not other myeloid populations".
cl_levels    <- levels(factor(mac$Celltype_subset))
macro_levels <- grep("_Macro_", cl_levels, value = TRUE)
if (include_cycling) macro_levels <- c(macro_levels, grep("_Cycling_", cl_levels, value = TRUE))
macro_levels <- cl_levels[cl_levels %in% macro_levels]   # keep the object's cluster ordering
message("Macrophage clusters used: ", paste(macro_levels, collapse = ", "))

# Macrophage cells, pre-treatment only, with sample + response under stable names.
# The pre-treatment filter must be explicit: new_group carries post-pCR / post-RD, and
# relying on factor(levels = c("pre-pCR","pre-RD")) would convert those cells to NA
# rather than removing them.
macro_meta <- mac@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::filter(Celltype_subset %in% macro_levels,
                .data[[response_col]] %in% c("pre-pCR", "pre-RD")) %>%
  dplyr::mutate(sample   = as.character(.data[[sample_col]]),
                response = factor(as.character(.data[[response_col]]),
                                  levels = c("pre-pCR", "pre-RD")))

# One M1 and one M2 value per sample = the pseudoreplication fix. Tests below have
# n = samples, not n = cells, matching the CD8 patient-level reanalysis.
sample_scores <- macro_meta %>%
  dplyr::group_by(sample, response) %>%
  dplyr::summarise(n_mac = dplyr::n(),
                   M1    = mean(M1_score1, na.rm = TRUE),
                   M2    = mean(M2_score1, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::filter(n_mac >= min_cells) %>%
  dplyr::mutate(polarization = M1 - M2)   # >0 = M1-skewed sample

# Sample accounting, for the response letter:
#   - 25 pre-treatment samples exist
#   - P02R and P06R contain ZERO macrophages (2 and 5 myeloid cells respectively)
#   - P17R contains 2 macrophages -> dropped by min_cells = 10
#   - leaves 7 pCR vs 15 RD. The macrophage counts jump 2 -> 12 -> 16 -> 20, so a threshold
#     anywhere in 3..12 gives the same result; 10 is not a knife-edge choice.
message("Samples retained: ",
        sum(sample_scores$response == "pre-pCR"), " pCR, ",
        sum(sample_scores$response == "pre-RD"),  " RD")
write.csv(sample_scores, "output/Macrophage_M1_M2_sample_scores.csv", row.names = FALSE)

# Restrict the cell table to the samples that survived the count filter
macro_meta <- dplyr::filter(macro_meta, sample %in% sample_scores$sample)


### PRIMARY: Macrophage M1/M2 by sample, pCR vs RD -> output/Macrophage_M1_M2_by_sample.pdf ###
# Structurally identical to Fig3C_CD8T_boxplots.pdf: signature.matrix -> per-sample means ->
# ggboxplot + labeled jitter + Wilcoxon, one panel per signature. Only the inputs differ
# (M1/M2 module scores instead of the CD8 function signatures; macrophage cells only).

mac.sub <- dplyr::rename(macro_meta, group = response)   # `group` matches the Fig 3C loop's naming

signature.matrix <- mac.sub[, c("M1_score1", "M2_score1")]
sig_titles       <- c(M1_score1 = "M1-like signature", M2_score1 = "M2-like signature")
group  <- mac.sub$group
sample <- mac.sub$sample

plots <- list()
raw_p <- numeric(ncol(signature.matrix))

for (i in 1:ncol(signature.matrix)) {
  
  exp.matrix <- data.frame(expression = signature.matrix[, i],
                           group = group, sample = sample)
  
  # One value per sample: mean module score across that sample's macrophages
  sample.means <- exp.matrix %>%
    group_by(sample, group) %>%
    summarise(expression = mean(expression, na.rm = TRUE), .groups = "drop")
  
  # Captured separately so it can be BH-corrected across the two signatures below;
  # stat_compare_means() on the panel prints the RAW p, as Fig 3C did.
  raw_p[i] <- wilcox.test(expression ~ group, data = sample.means)$p.value
  
  # Axis labels carry the per-group n, so the reader sees the sample size on the figure
  n_labels <- sample.means %>%
    group_by(group) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    mutate(label = paste0(group, "\n(n=", n, ")"))
  label_map <- setNames(n_labels$label, n_labels$group)
  
  # Push pCR labels left and RD labels right so ggrepel doesn't stack them over the boxes
  sample.means <- sample.means %>%
    mutate(nudge = ifelse(group == "pre-pCR", -0.6, 0.6))
  
  plots[[i]] <- ggboxplot(
    sample.means, x = "group", y = "expression", fill = "group",
    palette = "nejm", add = "jitter",
    add.params = list(size = 1.5, alpha = 0.8)
  ) +
    # Label every point with its sample ID -- at n = 7 vs 15, individual samples are
    # identifiable and a reviewer will want to trace outliers back to a patient
    geom_text_repel(
      data = sample.means,
      aes(x = group, y = expression, label = sample),
      size = 1.8, max.overlaps = Inf,
      box.padding = 0.2, point.padding = 0.2,
      segment.size = 0.25, segment.color = "grey60", segment.alpha = 0.6,
      direction = "y", nudge_x = sample.means$nudge,
      min.segment.length = 0, force = 2
    ) +
    scale_x_discrete(labels = label_map) +
    ylab("Score") + xlab(NULL) +
    ggtitle(sig_titles[colnames(signature.matrix)[i]]) +
    theme_classic(base_size = 11) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text.x = element_text(size = 9, face = "bold"),
          legend.position = "none",
          plot.margin = margin(5, 40, 5, 40)) +   # extra L/R margin so nudged labels aren't clipped
    stat_compare_means(label = "p.format", method = "wilcox.test",
                       label.x.npc = 0.65, label.y.npc = 0.97, size = 3.5)
}

signature.score.plot <- wrap_plots(plots = plots, ncol = 2)
ggsave("output/Macrophage_M1_M2_by_sample.pdf", signature.score.plot,
       width = 7, height = 5)

# BH across the two signatures. Quote THESE in the text -- the on-panel numbers are raw p,
# and every other test in the manuscript is FDR-corrected.
mac_stats <- data.frame(signature = c("M1-like", "M2-like"),
                        p_raw = raw_p,
                        p_adj = p.adjust(raw_p, method = "BH"))
print(mac_stats)
write.csv(mac_stats, "output/Macrophage_M1_M2_by_sample_stats.csv", row.names = FALSE)


### SUPPLEMENTARY: Macrophage M1/M2 dot-heatmap by sample -> output/Macrophage_M1_M2_dotheatmap_by_sample.pdf ###

# Drop genes essentially undetected in this compartment
genes_all <- c(m1_genes, m2_genes)
pct_all   <- FetchData(mac, vars = genes_all)[macro_meta$cell, , drop = FALSE] %>%
  summarise(across(everything(), ~ mean(.x > 0))) %>%
  unlist()

dropped <- names(pct_all)[pct_all < min_pct]
if (length(dropped))
  message("Dropped (expressed in <", min_pct * 100, "% of macrophages): ",
          paste(dropped, collapse = ", "))

genes_use     <- genes_all[!genes_all %in% dropped]
gene_category <- setNames(ifelse(genes_use %in% m1_genes, "M1-like", "M2-like"), genes_use)

# Sanity check: rows (now columns) sorted by polarization
ct <- suppressWarnings(cor.test(sample_scores$n_mac, sample_scores$polarization,
                                method = "spearman"))
message("n_mac vs polarization: rho = ", signif(ct$estimate, 2),
        ", p = ", signif(ct$p.value, 2))

# Per-sample, per-gene summary
dot_by_sample <- FetchData(mac, vars = genes_use) %>%
  rownames_to_column("cell") %>%
  dplyr::inner_join(dplyr::select(macro_meta, cell, sample, response), by = "cell") %>%
  pivot_longer(dplyr::all_of(genes_use), names_to = "gene", values_to = "expression") %>%
  dplyr::group_by(gene, sample, response) %>%
  dplyr::summarise(mean_expr   = mean(expression, na.rm = TRUE),
                   pct_express = mean(expression > 0, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(z_score = scale(mean_expr)[, 1]) %>%
  dplyr::ungroup()

# Align sample levels and split vectors
sample_levels <- sample_scores %>%
  dplyr::arrange(response, dplyr::desc(polarization)) %>%
  dplyr::pull(sample)
col_split_s <- factor(sample_scores$response[match(sample_levels, sample_scores$sample)],
                      levels = c("pre-pCR", "pre-RD"))   

# Original matrices [samples, genes]
to_mat <- function(val) dot_by_sample %>%
  dplyr::select(sample, gene, dplyr::all_of(val)) %>%
  pivot_wider(names_from = gene, values_from = dplyr::all_of(val)) %>%
  column_to_rownames("sample") %>% as.matrix()
zmat <- to_mat("z_score")[sample_levels, genes_use, drop = FALSE]
pmat <- to_mat("pct_express")[sample_levels, genes_use, drop = FALSE]

# --- Transpose matrices so genes are rows, samples are columns ---
zmat_t <- t(zmat)
pmat_t <- t(pmat)

# Colors + dot geometry 
col_fun_s <- colorRamp2(c(-2, 0, 2), c("#5C88DAFF", "white", "#CC0C00FF"))
cell_mm <- 7; max_r <- 2.45   

# PI FIX: Change M1/M2 to neutral tones (Charcoal & Light Grey) 
# so Red/Blue can be strictly reserved for patient response groups.
signature_cols <- c("M1-like" = "#4D4D4D", "M2-like" = "#C0C0C0")

# Updated cell function indexing the transposed pmat_t [gene_row, sample_col]
cell_fun_s <- function(j, i, x, y, width, height, fill) {
  grid.rect(x, y, width, height, gp = gpar(col = "grey92", fill = NA, lwd = 0.4))
  pe <- pmat_t[i, j]
  if (!is.na(pe) && pe > 0)
    grid.circle(x, y, r = unit(pe * max_r, "mm"),
                gp = gpar(fill = fill, col = "grey35", lwd = 0.3))
}

# --- AXIS SWAP: Top Annotation (now tracks SAMPLES along the columns) ---
# PI FIX: Added "Response" bar colored with your shared red/blue group colors
n_mac_vec <- sample_scores$n_mac[match(sample_levels, sample_scores$sample)]
top_anno_s <- HeatmapAnnotation(
  Response = col_split_s,
  `N macrophages` = anno_barplot(n_mac_vec, 
                                 gp = gpar(fill = "grey45", col = NA),
                                 height = unit(1.6, "cm"),
                                 axis_param = list(gp = gpar(fontsize = 6))),
  col = list(Response = group_colors), # Automatically maps red to pre-pCR and blue to pre-RD
  annotation_name_gp   = gpar(fontsize = 7),
  annotation_name_side = "left",
  gap = unit(1.5, "mm")
)

# --- AXIS SWAP: Left Annotation (now tracks GENES along the rows) ---
left_anno_s <- rowAnnotation(
  `Mean % expressing` = anno_barplot(rowMeans(pmat_t, na.rm = TRUE),
                                     gp = gpar(fill = "grey55", col = NA),
                                     width = unit(1, "cm"),
                                     axis_param = list(gp = gpar(fontsize = 6))),
  Signature = gene_category[rownames(zmat_t)],
  col = list(Signature = signature_cols), # Uses the new neutral colors
  simple_anno_size     = unit(3.5, "mm"),
  show_annotation_name = c(TRUE, FALSE),   
  annotation_name_side = "bottom",
  annotation_name_gp   = gpar(fontsize = 7),
  gap = unit(1, "mm")
)

# Render Heatmap with swapped configurations
ht_s <- Heatmap(
  zmat_t, col = col_fun_s,
  rect_gp  = gpar(type = "none"),          
  cell_fun = cell_fun_s,
  cluster_rows = FALSE, cluster_columns = FALSE,   
  column_split = col_split_s, 
  column_gap = unit(2, "mm"),
  column_title_gp = gpar(fontsize = 9, fontface = "bold"), 
  row_split = factor(gene_category[rownames(zmat_t)], levels = c("M1-like", "M2-like")),
  row_gap = unit(2, "mm"),
  row_title_gp = gpar(fontface = "bold", fontsize = 10), row_title_rot = 0,
  top_annotation  = top_anno_s,
  left_annotation = left_anno_s,
  width  = unit(ncol(zmat_t) * cell_mm, "mm"),
  height = unit(nrow(zmat_t) * cell_mm, "mm"),
  column_names_gp = gpar(fontsize = 8),
  row_names_side  = "left", row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  name = "Z-score\n(mean expr)",
  border = TRUE
)

# Legend
pct_breaks <- c(0.25, 0.5, 0.75, 1.0)
lgd_size_s <- Legend(
  title = "% Expressed", labels = paste0(pct_breaks * 100, "%"),
  type = "points", pch = 21, size = unit(pct_breaks * max_r * 2, "mm"),
  legend_gp = gpar(fill = "grey55", col = "grey35"), background = "white",
  title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8)
)

# Swapped layout dimensions (wider rather than taller, scaling with sample columns)
pdf("output/Macrophage_M1_M2_dotheatmap_by_sample.pdf",
    width = 4 + 0.30 * ncol(zmat_t), height = 6.5)

draw(ht_s,
     column_title = "Macrophage M1/M2 marker expression by sample and response",
     column_title_gp = gpar(fontface = "bold", fontsize = 13),
     annotation_legend_list = list(lgd_size_s),
     merge_legend = TRUE,
     heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(12, 12, 2, 2), "mm"))   # bottom, left, top, right

grid.text(
  paste("Columns = samples (macrophage clusters only), sorted by M1-M2 polarization within response.",
        "\nZ-scores computed across all samples per gene, so pCR and RD blocks are directly comparable."),
  x = unit(0.5, "npc"), y = unit(0.025, "npc"),
  gp = gpar(fontsize = 8, fontface = "italic")
)
dev.off()

### SUPPLEMENTARY: Macrophage M1 vs M2 Scatterplot -> output/Macrophage_M1_M2_scatterplot.pdf ###
library(ggplot2)
library(ggpubr) # Optional, for adding correlation coefficients easily

# 1. Calculate the sample-level average expression for M1 and M2 signatures
# (Using the same pre-filtered 'genes_use' from your heatmap script)
m1_use <- genes_use[genes_use %in% m1_genes]
m2_use <- genes_use[genes_use %in% m2_genes]

scatter_data <- FetchData(mac, vars = genes_use) %>%
  rownames_to_column("cell") %>%
  dplyr::inner_join(dplyr::select(macro_meta, cell, sample, response), by = "cell") %>%
  pivot_longer(dplyr::all_of(genes_use), names_to = "gene", values_to = "expression") %>%
  # Determine which signature each gene belongs to
  dplyr::mutate(signature = ifelse(gene %in% m1_use, "M1_score", "M2_score")) %>%
  # Average expression per sample, per signature
  dplyr::group_by(sample, response, signature) %>%
  dplyr::summarise(mean_val = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = signature, values_from = mean_val) %>%
  # Ensure response is ordered correctly
  dplyr::mutate(response = factor(response, levels = c("pre-pCR", "pre-RD")))

# 2. Build the ggplot
# We use custom colors matching your previous theme: red for pCR, blue for RD (or vice versa)
response_cols <- c("pre-pCR" = "#CC0C00FF", "pre-RD" = "#5C88DAFF")

p_scatter <- ggplot(scatter_data, aes(x = M1_score, y = M2_score)) +
  # Add a light background grid for readability
  theme_bw(base_size = 11) +
  # Add overall trend line (linear regression)
  geom_smooth(method = "lm", formula = y ~ x, color = "grey40", fill = "grey90", alpha = 0.5, linetype = "dashed") +
  # Plot sample points colored/shaped by clinical response
  geom_point(aes(fill = response, shape = response), size = 4.5, color = "grey30", stroke = 0.5) +
  # Style shapes: 21 (circle) and 24 (triangle) allow both border and fill coloring
  scale_shape_manual(values = c("pre-pCR" = 21, "pre-RD" = 24)) +
  scale_fill_manual(values = response_cols) +
  # Labels
  labs(
    title = "Macrophage M1 vs. M2 Marker Expression",
    subtitle = "Pre-treatment samples by clinical response",
    x = "Mean M1-like Expression",
    y = "Mean M2-like Expression",
    fill = "Response",
    shape = "Response"
  ) +
  # Clean up theme elements to match publication-quality layout
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", size = 10, hjust = 0.5, color = "grey30"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "grey80", size = 0.5)
  ) +
  # Add statistical annotation (Pearson or Spearman correlation)
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.5, fontface = "italic")

# 3. Save the plot
pdf("output/Macrophage_M1_M2_scatterplot.pdf", width = 5.5, height = 4.5)
print(p_scatter)
dev.off()