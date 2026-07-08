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
sct.scale.all = GetAssayData(scRNA.pre, slot="scale.data",assay="SCT")
cell_counts <- table(scRNA.pre$new_id_4)
sample <- names(cell_counts[cell_counts > 0])
score.list<-list()
for (i in 1:length(sample)){
  score.df<-data.frame("gene"=rownames(sct.scale.all))
  scRNA.sub<-subset(scRNA,new_id_4 %in% sample[i])
  sct.scale.matrix = GetAssayData(scRNA.sub, slot="scale.data",assay="SCT")
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

#### LOGESTIC REGRESSION ####
str(scRNA.origin, max.level = 2)
colnames(scRNA.origin@meta.data)
colnames(scRNA.cancer@meta.data)

library(readxl)
clinical <- read_excel("~/Documents/Ting Data/Table S1_baseline clinical characteristics.xlsx")
head(clinical)
colnames(clinical)

unique(scRNA.cancer@meta.data$new_group)
table(scRNA.cancer@meta.data$new_group)

unique(scRNA.cancer@meta.data$new_id_4) |> head(20)

# 1. Pull metadata and filter to pre-treatment only
meta <- scRNA.cancer@meta.data %>%
  filter(new_group %in% c("pre-pCR", "pre-RD"))

# 2. Aggregate MP scores per patient (mean across cells)
patient_df <- meta %>%
  group_by(new_id_4, new_group) %>%
  summarise(
    MP_11 = mean(MP_11, na.rm = TRUE),
    MP_22 = mean(MP_22, na.rm = TRUE),
    MP_33 = mean(MP_33, na.rm = TRUE),
    MP_44 = mean(MP_44, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pCR = ifelse(new_group == "pre-pCR", 1, 0))

# 3. Sanity Check
table(patient_df$pCR)

# 4. Logistic regression
library(logistf)
model_firth <- logistf(pCR ~ MP_11 + MP_22 + MP_33 + MP_44,
                       data = patient_df)
summary(model_firth)
head(patient_df)
table(patient_df$pCR)

#possibly need to be corrected (p values)
# 5. testing correlation 
mp_cols <- c("MP_11", "MP_22", "MP_33", "MP_44")

pvals <- sapply(mp_cols, function(mp) {
  wilcox.test(patient_df[[mp]] ~ patient_df$pCR)$p.value
})

pvals

library(ggplot2)

ggplot(patient_df, aes(x = factor(pCR, labels = c("RD", "pCR")), 
                       y = MP_22, 
                       fill = factor(pCR, labels = c("RD", "pCR")))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("RD" = "#4575b4", "pCR" = "#d73027")) +
  labs(
    title = "MP_22 score by treatment response",
    subtitle = paste0("RD n=", sum(patient_df$pCR==0), ", pCR n=", sum(patient_df$pCR==1)),
    x = "Response",
    y = "MP_22 score (patient mean)",
    fill = "Response"
  ) +
  theme_classic() +
  theme(legend.position = "none")
ggsave("output/MP22_by_response.pdf", width = 4, height = 5)



#### PIE CHART AND BOX PLOTS #####
ls()[grep("CD8|T.cell|Tcell|immune", ls(), ignore.case = TRUE)]
unique(scRNA.origin@meta.data$Celltype_Minor)


# Check for separate CD8 object
ls()[grep("CD8|cd8", ls(), ignore.case = TRUE)]

# Check Celltype_Subset for CD8 cluster labels (c0, c1, c2 etc)
scRNA.origin@meta.data %>%
  filter(Celltype_Minor == "CD8+T cell") %>%
  pull(Celltype_Subset) %>%
  unique()

# Subset to CD8T cells, pre-treatment only
cd8_meta <- scRNA.origin@meta.data %>%
  filter(Celltype_Minor == "CD8+T cell",
         new_group %in% c("pre-pCR", "pre-RD"))

# Get cluster order and assign colors
cd8_clusters <- unique(cd8_meta$Celltype_Subset)
n_clusters <- length(cd8_clusters)

umapColor <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175', '#FF7F0E', '#1F77B4', '#2CA02C', '#D62728', '#9467BD')

cluster_colors <- setNames(umapColor[1:n_clusters], sort(cd8_clusters))

# Calculate proportions per group
cd8_prop <- cd8_meta %>%
  group_by(new_group, Celltype_Subset) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(new_group) %>%
  mutate(prop = n / sum(n))

# --- Pie chart ---
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

# --- Stacked bar ---
p_bar_comp <- ggplot(cd8_prop, aes(x = new_group, y = prop, fill = Celltype_Subset)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(title = "CD8+ T cell subset composition by response",
       x = "Response Group", y = "Proportion", fill = "Cluster") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        panel.grid.major.y = element_line(linetype = "dashed", color = "grey85"))


# --- Boxplot by cluster ---
p_boxcluster <- ggplot(expr_long, aes(x = cluster, y = expression, fill = new_group)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  scale_fill_manual(values = c("pre-pCR" = "#CC0C00FF", "pre-RD" = "#5C88DAFF")) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  labs(title = "CD8+ T cell gene expression by cluster and response",
       x = "Cluster", y = "Normalized expression", fill = "Response") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = function(x) gsub("_", "\n", x))

# --- Dot plot ---
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
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        panel.grid.major = element_line(color = "#EFEFEF", linewidth = 0.2),
        panel.grid.minor = element_blank())

# --- Patient-level beeswarm ---
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

# --- Save all explicitly ---
ggsave("output/CD8T_subset_piechart.pdf",          plot = p_pie,        width = 12, height = 6)
ggsave("output/CD8T_subset_barplot.pdf",           plot = p_bar_comp,   width = 8,  height = 6)
ggsave("output/CD8T_gene_expression_boxplot.pdf",  plot = p_boxcluster, width = 14, height = 10)
ggsave("output/CD8T_marker_dotplot.pdf",           plot = p_dot,        width = 10, height = 6)
ggsave("output/CD8T_pretreatment_expression_6Astyle.pdf", plot = p_beeswarm, width = 8, height = 10)
ggsave("output/Macrophage_M1_M2_UMAP.pdf",        plot = p1 | p2 | p3 | p4, width = 22, height = 6)

genes <- c("IFNG", "TNF", "GZMK", "GZMB", "GZMH", "MKI67")

expr_df <- FetchData(scRNA.origin, 
                     vars = c("Celltype_Minor", "Celltype_Subset", 
                              "new_group", genes)) %>%
  filter(Celltype_Minor == "CD8+T cell",
         new_group %in% c("pre-pCR", "pre-RD"))

head(expr_df)

library(tidyr)

# Compute per-group dot plot data
dotplot_grouped <- expr_long %>%
  group_by(gene, cluster, new_group) %>%
  summarise(
    mean_expr   = mean(expression, na.rm = TRUE),
    pct_express = mean(expression > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(gene, new_group) %>%
  mutate(z_score = scale(mean_expr)[,1]) %>%
  ungroup() %>%
  mutate(new_group = factor(new_group, levels = c("pre-pCR", "pre-RD")))

ggplot(dotplot_grouped, aes(x = gene, y = cluster)) +
  geom_point(aes(size = pct_express, color = z_score)) +
  scale_color_gradient2(
    low      = "#5C88DAFF",
    mid      = "white",
    high     = "#CC0C00FF",
    midpoint = 0,
    name     = "Z-score\n(mean expr)"
  ) +
  scale_size_continuous(
    range  = c(0.5, 7),
    labels = scales::percent_format(),
    name   = "% Expressed"
  ) +
  facet_wrap(~ new_group, ncol = 2) +
  labs(
    title = "CD8+ T cell marker gene expression by subtype and response",
    x     = "Marker Gene",
    y     = "CD8+ Subtype"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y      = element_text(size = 9),
    axis.title       = element_text(size = 11),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
    strip.text       = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size = 9)
  )

ggsave("output/CD8T_dotplot_by_response.pdf", width = 12, height = 7)

### UMAP ###
table(scRNA.origin@meta.data$Celltype_Minor)
mac <- subset(scRNA.origin, Celltype_Minor == "Macrophage")
m1_genes <- c("TNF", "IL6", "IL1B", "CXCL10", "CXCL9", "CD86", "NOS2", "STAT1")
m2_genes <- c("MRC1", "CD163", "ARG1", "IL10", "TGFB1", "CCL18", "VEGFA", "PDCD1LG2")

# Filter to genes present in the object
m1_genes <- m1_genes[m1_genes %in% rownames(mac)]
m2_genes <- m2_genes[m2_genes %in% rownames(mac)]

# Add module scores
mac <- AddModuleScore(mac, features = list(m1_genes), name = "M1_score")
mac <- AddModuleScore(mac, features = list(m2_genes), name = "M2_score")

m1_p <- wilcox.test(M1_score1 ~ Celltype_Subset, 
                    data = mac@meta.data %>% dplyr::filter(grepl("c4|c11", Celltype_Subset)))$p.value
m2_p <- wilcox.test(M2_score1 ~ Celltype_Subset, 
                    data = mac@meta.data %>% dplyr::filter(grepl("c4|c11", Celltype_Subset)))$p.value

padj <- p.adjust(c(m1_p, m2_p), method = "BH")
names(padj) <- c("M1-like", "M2-like")

fdr_labels <- data.frame(
  signature = c("M1-like", "M2-like"),
  label = paste0("FDR = ", signif(padj, 2))
)

# Merge into score_df for annotating
score_df <- score_df %>%
  left_join(fdr_labels, by = "signature")

unique(score_df$cluster)


# take the macrophages out 
# Plot
p1 <- FeaturePlot(mac, features = "M1_score1", order = TRUE, pt.size = 0.5) +
  scale_color_gradientn(colors = c("lightgrey", "#CC0C00FF")) +
  labs(title = "M1-like signature") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- FeaturePlot(mac, features = "M2_score1", order = TRUE, pt.size = 0.5) +
  scale_color_gradientn(colors = c("lightgrey", "#5C88DAFF")) +
  labs(title = "M2-like signature") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Add annotation to UMAP plots
p1 <- FeaturePlot(mac, features = "M1_score1", order = TRUE, pt.size = 0.5) +
  scale_color_gradientn(colors = c("lightgrey", "#CC0C00FF")) +
  labs(title = "M1-like signature",
       caption = format_p(p_m1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(hjust = 0.5, size = 10, face = "italic"))

p2 <- FeaturePlot(mac, features = "M2_score1", order = TRUE, pt.size = 0.5) +
  scale_color_gradientn(colors = c("lightgrey", "#5C88DAFF")) +
  labs(title = "M2-like signature",
       caption = format_p(p_m2)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(hjust = 0.5, size = 10, face = "italic"))

p1 <- p1 + labs(caption = paste0("FDR = ", signif(padj["M1-like"], 2)))
p2 <- p2 + labs(caption = paste0("FDR = ", signif(padj["M2-like"], 2)))

p3 <- DimPlot(mac, group.by = "Celltype_Subset",
              label = TRUE, label.size = 3, repel = TRUE, pt.size = 0.5) +
  labs(title = "Macrophage clusters") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

y_pos <- score_df %>%
  group_by(signature) %>%
  summarise(y = max(score, na.rm = TRUE) * 1.08, .groups = "drop") %>%
  left_join(fdr_labels, by = "signature")


p4 <- ggplot(score_df, aes(x = signature, y = score, fill = cluster)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8, width = 0.6,
               position = position_dodge(0.75)) +
  geom_text(
    data = y_pos,
    aes(x = signature, y = y, label = label),
    inherit.aes = FALSE,
    size = 3, fontface = "italic"
  ) +
  scale_fill_manual(values = c("c4_Macro_CCL4" = "#CC0C00FF", 
                               "c11_Cycling_MKI67" = "#5C88DAFF"))
  labs(title = "M1 vs M2 scores: c4 & c11",
       x = NULL, y = "Module score", fill = "Cluster") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "right")

pdf("output/Macrophage_M1_M2_UMAP.pdf", width = 22, height = 6)
p1 | p2 | p3 | p4
dev.off()

## fdr for p values
data.frame(MP = mp_cols, p_raw = pvals, p_adj = p.adjust(pvals, method = "BH"))
mp_results <- data.frame(
  MP    = mp_cols,
  p_raw = pvals,
  p_adj = p.adjust(pvals, method = "BH")
) %>%
  mutate(significant = p_adj < 0.05)

write.csv(mp_results, "output/MP_wilcox_FDR.csv", row.names = FALSE)

# save plots
ggsave("output/Macrophage_M1_M2_UMAP.pdf", width = 18, height = 6)
ggsave("output/CD8T_pretreatment_expression_6Astyle.pdf", width = 8, height = 10)
ggsave("output/CD8T_subset_piechart.pdf", width = 12, height = 6)
ggsave("output/CD8T_gene_expression_boxplot.pdf", width = 14, height = 10)
ggsave("output/CD8T_marker_dotplot.pdf", width = 10, height = 6)