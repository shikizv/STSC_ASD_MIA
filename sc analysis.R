library(Seurat)
library(tidyverse)
library(dplyr)
library(harmony)
library(ggplot2)
library(ggsci)
library(dplyr)
library(reshape2)
library(DoubletFinder)
library(patchwork)
library(dittoSeq)
library(magrittr)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(pheatmap)
library(clustree)


col1 <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_lancet()(9))[-c(8,12,17,20,21)]#UMAP


####Data load####
set.seed(123)

dir = c('./data/SC/ASD/poly_IC_1/filtered_feature_bc_matrix/', 
        './data/SC/ASD/poly_IC_3/filtered_feature_bc_matrix/',
        './data/SC/ASD/poly_IC_4/filtered_feature_bc_matrix/',
        './data/SC/ASD/WT_1/filtered_feature_bc_matrix/',
        './data/SC/ASD/WT_2/filtered_feature_bc_matrix/',
        './data/SC/ASD/WT_3/filtered_feature_bc_matrix/')

names(dir) = c('poly_IC_1', 'poly_IC_3', 'poly_IC_4', 'WT_1', 'WT_2', 'WT_3')

ASD.sc <- list()

for (i in 1:length(dir)) {
  # ASD.sc[[i]] <- CreateSeuratObject(counts =  Read10X(data.dir = dir[i], gene.column = 2),
  #                                     project = names(dir)[i])
  ASD.sc[[i]]$orig.ident <- names(dir)[i]
}


# dir = c('./data/SC/ASD/poly_IC_1/raw_feature_bc_matrix/', 
#         './data/SC/ASD/poly_IC_3/raw_feature_bc_matrix/',
#         './data/SC/ASD/poly_IC_4/raw_feature_bc_matrix/',
#         './data/SC/ASD/WT_1/raw_feature_bc_matrix/',
#         './data/SC/ASD/WT_2/raw_feature_bc_matrix/',
#         './data/SC/ASD/WT_3/raw_feature_bc_matrix/')
# ASD.sc.raw <- list()
# for (i in 1:length(dir)) {
#   ASD.sc.raw[[i]] <- CreateSeuratObject(counts =  Read10X(data.dir = dir[i], gene.column = 2),
#                                       project = names(dir)[i])
#   # ASD.sc.raw[[i]]$orig.ident <- names(dir)[i]
# }

#raw lists of single samples
save(ASD.sc, file = './object/sc/ASD.sc.rawlist.RData')
rm(ASD.sc)


ASD.sce <- merge(x=ASD.sc[[1]], y=ASD.sc[-1]) 

DotPlot(ASD.sce, features = 'Igfbp7', group.by = 'orig.ident')
DotPlot(ASD.sce, features = 'B2m', group.by = 'orig.ident')

table(ASD.sce$orig.ident)

VlnPlot(ASD.sce, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(ASD.sce, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
ASD.sce[["percent.mt"]] <- PercentageFeatureSet(ASD.sce, pattern = "^mt-")
VlnPlot(ASD.sce, features = "percent.mt", pt.size = 0.1) + NoLegend()

# BiocManager::install("scDblFinder")
library(scDblFinder)
library(BiocParallel)

ASD.sce.df <- as.SingleCellExperiment(ASD.sce)
ASD.sce.df <- scDblFinder(ASD.sce.df, samples="orig.ident", BPPARAM=MulticoreParam(3))
table(ASD.sce.df$scDblFinder.class)

ASD.sce.qc <- as.Seurat(ASD.sce.df)
table(ASD.sce.qc$orig.ident)
table(ASD.sce.qc$scDblFinder.class)
ASD.sce.qc <- subset(ASD.sce.qc, scDblFinder.class == 'singlet')

ASD.sce.qc <- subset(ASD.sce.qc, subset = nFeature_RNA > 800 & nFeature_RNA < 6000 & percent.mt < 10 & 
                       nCount_RNA > 1000 & nCount_RNA < 20000)
colnames(ASD.sce.qc@meta.data)

p <- VlnPlot(ASD.sce.qc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 1,
        group.by = 'orig.ident') & NoLegend() & scale_fill_npg()

pdf('./figures/SC/main/Vlnplot_QC.pdf', width = 4, height = 7)
print(p)
dev.off()

# metadata <- ASD.sce.qc@meta.data
ASD.sce.qc@meta.data <- ASD.sce.qc@meta.data[,c(1:5)]
ASD.sce.qc <- SCTransform(ASD.sce.qc, assay = "RNA", vars.to.regress = c("percent.mt"), 
                          ncells = round(dim(ASD.sce.qc)[2]/10))
# ASD.sce.qc <- SCTransform(ASD.sce.qc, assay = "RNA", vars.to.regress = "percent.mt")
ASD.sce.qc <- RunPCA(ASD.sce.qc)
ASD.sce.qc <- RunHarmony(ASD.sce.qc, group.by.vars = "orig.ident")
ASD.sce.qc <- RunUMAP(ASD.sce.qc,reduction = "harmony",dims = 1:30, min.dist = 0.4)#
# ASD.sce.qc <- RunTSNE(ASD.sce.qc,reduction = "harmony",dims = 1:30)
ASD.sce.qc <- FindNeighbors(ASD.sce.qc, reduction = "harmony", dims = 1:30)
ASD.sce.qc <- FindClusters(ASD.sce.qc, 
                        resolution =0.4,
                        method ='igraph',
                        verbose = T)
#0.1 17 0.2 22 

ASD.sce.qc <- subset(ASD.sce.qc, seurat_clusters != '19')
colnames(ASD.sce.qc@meta.data)
# metadata <- ASD.sce.qc@meta.data
ASD.sce.qc@meta.data <- ASD.sce.qc@meta.data[,c(1:5)]
ASD.sce.qc <- SCTransform(ASD.sce.qc, assay = "RNA", vars.to.regress = c("percent.mt"), 
                          ncells = round(dim(ASD.sce.qc)[2]/10))
# ASD.sce.qc <- SCTransform(ASD.sce.qc, assay = "RNA", vars.to.regress = "percent.mt")
ASD.sce.qc <- RunPCA(ASD.sce.qc)
ASD.sce.qc <- RunHarmony(ASD.sce.qc, group.by.vars = "orig.ident")
ASD.sce.qc <- RunUMAP(ASD.sce.qc,reduction = "harmony",dims = 1:30, min.dist = 0.3)#
# ASD.sce.qc <- RunTSNE(ASD.sce.qc,reduction = "harmony",dims = 1:30)
ASD.sce.qc <- FindNeighbors(ASD.sce.qc, reduction = "harmony", dims = 1:30)
ASD.sce.qc <- FindClusters(ASD.sce.qc, 
                           resolution =0.4,
                           method ='igraph',
                           verbose = T)
#0.1 0.2 21 0.3 25 0.4 27
Idents(ASD.sce.qc) <- ASD.sce.qc$SCT_snn_res.4
DimPlot(ASD.sce.qc, cols = pal_igv()(50))#
DimPlot(ASD.sce.qc, cols = pal_igv()(50), group.by = 'seurat_clusters')#
# DimPlot(ASD.sce.qc, reduction = 'tsne', cols = pal_igv()(50))#


VlnPlot(ASD.sce.qc, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(ASD.sce.qc, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(ASD.sce.qc, features = "percent.mt", pt.size = 0.1) + NoLegend()

dittoBarPlot(ASD.sce.qc, 'seurat_clusters', group.by = 'orig.ident')


# save(ASD.sce.qc, file = './object/sc/ASD.sce.qc_df_new.RData')

Idents(ASD.sce.qc) <- ASD.sce.qc$seurat_clusters
Markers_all <- FindAllMarkers(ASD.sce.qc, assay = "SCT", logfc.threshold = 0.25, only.pos = T)
Markers_all_top10 <- Markers_all %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Markers_all_top10, file = './results/SC/Markers_all_top10.csv')
save(Markers_all, file = './object/sc/Markers_all.RData')

DotPlot(Brian.sc_ys, features = Markers_all_top10[which(Markers_all_top10$cluster == '11'),]$gene)

DotPlot(ASD.sce.qc, features = c("Hexb", "C1qc", "Ctss", "Tmem119", "Itgam","Cx3cr1",#Microgila 
                                 "Gja1", "Aldoc", "Plpp3", "Slc1a3", "Slc1a2", "Acsl6",  #Astrocyte
                                 "Cldn11", "Mog", "Mbp", "Mag", "Ermn", "Plp1", "Enpp2", #OLG
                                 "Pdgfra", "Tnr", "Epn2", #OPC
                                 "Ccdc153", "Tmem212", "Cfap54", "Cfap44" #EPC
), group.by = 'seurat_clusters')#,group.by = 'predicted_MCA'

DotPlot(ASD.sce.qc, features = c("Siglech", "Ctss",#Microgila 
                                 "Aqp4", "Mfge8", #Astrocyte
                                 "Ptgds", "Plp1", #OLG
                                 "Pdgfra", "Nnat", #OPC
                                 'Slc6a13', #VLMC
                                 'Slc47a1', #ABC
                                 'Ttr' #CPC
), group.by = 'seurat_clusters')

# DotPlot(ASD.sce.qc, features = c("Kcnj8", "Vtn", "Rgs5",#pericyte
#                                  "Acta2", "Slc6a13", "Dcn", #Vascular
#                                  "Nf2", "Cryab", "Matn2", "Aldoc", #Schwann cels
#                                  "Cldn5", "Igfbp7", "Ptprb", 'Flt1',  #endothelial
#                                  'Ttr', 'Kl', 'Clic6', 'Prlr','Drc7' #choroid plexus
# ), group.by = 'seurat_clusters')

DotPlot(ASD.sce.qc, features = c("Prph", "Disp2", "Tubb3", "Rbfox3", #Neuron
  # "Epha7", "Stmn2", "Foxd2", "Sox11", "Epha3", #ImmunN
  "Sst", "Gpr88", "Baiap3", "Tac1", "Penk", #Nendc,
  "Syt1", "Snap25", "Meg3", "Atp1b1", #matureNeuron
  'Gad1', 'Gad2', 'Slc6a1', 'Gabbr2', #GABA
  'Slc17a7', 'Slc17a6', 'Slc17a8', 'Slc1a1', #Glu                                                                         
  'Slc6a9', 'Slc32a1' #gly
), group.by = 'seurat_clusters')
  
DotPlot(ASD.sce.qc, features = c("Pcp2", "Calb1", "Grik2", "Rora"

), group.by = 'seurat_clusters')



####genesorteR####
library(genesorteR)
library(dplyr)
gs <- sortGenes(ASD.sce.qc@assays$SCT@data, ASD.sce.qc$seurat_clusters)
# gs1 <- getMarkers(gs, quant = 0.99)
gs_spec <- gs$specScore
genes <- rownames(gs_spec)
gs_top_genes <- data.frame()
for (cluster in colnames(gs_spec)) {
  cluster_scores <- data.frame(
    gene = genes,
    specScore = gs_spec[, cluster],
    cluster = cluster
  )
  top_genes <- cluster_scores %>%
    arrange(desc(specScore)) %>%
    head(10)
  gs_top_genes <- rbind(gs_top_genes, top_genes)
}

DotPlot(ASD.sce.qc, features = unique(gs_top_genes$gene), group.by = 'seurat_clusters')
save(gs_top_genes, file = './object/sc/gs_topgenes.RData')



# ('Neurons', 'Microglia', 'Astrocytes', 'Oligodendrocytes', 'OPCs', 'VLMC', 'CPC', 'Endothelial cells', 'ABC', 'EPC')
Annocelltypes <- c("Neurons", "Neurons", "OLGs", "Neurons", "Neurons", "Astrocytes", "Neurons", "Neurons", "Neurons",
                   "Astrocytes", "Neurons", "Neurons", "Endothelial", "Microgila", "OPCs", "Neurons", 'Neurons', 'VLMCs',
                   "Neurons", "Neurons", "ABCs", "Neurons", "Neurons", "Neurons", "EPCs", "CPCs", "Neurons")

Idents(ASD.sce.qc) <- ASD.sce.qc$seurat_clusters
names(Annocelltypes)<- levels(ASD.sce.qc)
ASD.sce.qc <- RenameIdents(ASD.sce.qc, Annocelltypes)
ASD.sce.qc@meta.data$celltype_main <- ASD.sce.qc@active.ident

ASD.sce.qc$celltype_main <- factor(ASD.sce.qc$celltype_main, levels = row.names(sort(table(ASD.sce.qc$celltype_main), 
                                                                           decreasing = T)))

pdf('./figures/SC/main/Dimplot_celltype_all.pdf', width = 6, height = 5)
p1<-DimPlot(ASD.sce.qc, group.by = 'celltype_main', cols = col1)+ggtitle('celltype')
print(p1)
dev.off()

pdf('./figures/SC/main/Dimplot_celltype_all_newcol.pdf', width = 6, height = 5)
p1<-DimPlot(ASD.sce.qc, group.by = 'celltype_main', cols = sccol)+ggtitle('celltype')
print(p1)
dev.off()

pdf('./figures/SC/main/Dimplot_seurat_clusters_all_newcol.pdf', width = 6, height = 5)
p<-DimPlot(ASD.sce.qc, group.by = 'seurat_clusters', cols = col1)+ggtitle('clusters')
print(p)
dev.off()



ASD.sce.qc$group <- ASD.sce.qc$ident
ASD.sce.qc$group <- gsub('poly', 'poly_IC', ASD.sce.qc$group)
pdf('./figures/SC/main/Dimplot_group_all.pdf', width = 11, height = 10)
p1<-DimPlot(ASD.sce.qc, group.by = 'group', cols = col1, pt.size = 0.01)+ggtitle('celltype')+theme(
  axis.text.x=element_text(size = 20),axis.text.y=element_text(size = 20),
  axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), plot.title=element_text(size = 30))
print(p1)
dev.off()

pdf('./figures/SC/main/DittoBarPlot_celltype_all_newcol.pdf', width = 4, height = 5)
p2<-dittoBarPlot(ASD.sce.qc, 'celltype_main', group.by = 'group', main = '', color.panel = sccol, xlab = NULL,
                 var.labels.reorder = c(7,2,8,4,6,9,10,1,5,3))+theme(
                   axis.text.x=element_text(angle = 0, hjust = 0.5, vjust=0.5, size = 12),
                   axis.text.y=element_text(size = 10)
                 )
print(p2)
dev.off()

pdf('./figures/SC/main/Dotplot_markers_cf.pdf', width = 10, height = 5)
p3<-DotPlot(ASD.sce.qc, features = c("Syt1", "Snap25", "Meg3",#neuron
                                 "Plpp3", "Slc1a3", "Slc1a2", #Astrocyte
                                 "Mbp", "Mag", "Plp1", #OLG
                                 "Cldn5", "Ptprb", 'Flt1',#Endo
                                 "Hexb", "Ctss", "Itgam",#Microgila 
                                 "Pdgfra", "Tnr", "Epn2", #OPC
                                 'Slc6a13', #VLMC
                                 'Slc47a1', #ABC
                                 "Tmem212", "Cfap54", "Cfap44", #EPC
                                 'Ttr'#CPC
), group.by = 'celltype_main')+
  scale_color_gradientn(colours = c("#330066","#336699","#66CC66","#FFCC33"))+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 18,color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 16,color = "black"),
        axis.title = element_text(size = 16,color = 'black'),
        legend.text = element_text(size = 16,color = "black"))+labs(x="",y="")
print(p3)
dev.off()

######FeaturePlot######
feature_gene <- c("Syt1",
                  "Plpp3",
                  "Mbp",
                  "Cldn5",
                  "Hexb",
                  "Pdgfra",
                  'Slc6a13',
                  'Slc47a1',
                  "Tmem212",
                  'Ttr')
FeaturePlot(ASD.sce.qc, features = feature_gene, ncol = 2)&scale_color_viridis()
pdf('./figures/SC/main/FeaturePlot_all_celltype.pdf', width = 7, height = 15)
p4 <- FeaturePlot(ASD.sce.qc, features = feature_gene, ncol = 2)&scale_color_viridis()
print(p4)
dev.off()

pdf('./figures/SC/main/FeaturePlot_all_celltype11.pdf', width = 15, height = 5)
p4 <- FeaturePlot(ASD.sce.qc, features = feature_gene, ncol = 5)&scale_color_viridis()
print(p4)
dev.off()

#####OR#####



# DimPlot(ASD.sce.qc, group.by = 'seurat_clusters', cols = pal_igv()(50))


Annocelltypes <- c("ExN", "ExN", "OLGs", "ExN", "ExN", "Astrocytes", "ExN", "InN", "InN",
                   "Astrocytes", "InN", "MSN", "Endothelial", "Microgila", "OPCs", "ExN", 'ExN', 'VLMCs',
                   "InN", "ExN", "ABCs", "ExN", "InN", "InN", "EPCs", "CPCs", "ExN")

Idents(ASD.sce.qc) <- ASD.sce.qc$seurat_clusters
names(Annocelltypes)<- levels(ASD.sce.qc)
ASD.sce.qc <- RenameIdents(ASD.sce.qc, Annocelltypes)
ASD.sce.qc@meta.data$celltype <- ASD.sce.qc@active.ident
DimPlot(ASD.sce.qc, cols = pal_igv()(50), group.by = 'celltype')
DotPlot(ASD.sce.qc, features = c('Slc17a7', 'Slc17a6', 'Slc17a8', 'Camk2a', 'Neurod6', 'Satb2', 'Tbr1',
                                  'Gad2', "Gad1", "Slc32a1", 'Dlx1', 'Dlx2', 'Pvalb', 'Sst', 'Vip', 'Lhx6'
), group.by = 'celltype')

ASD.sce.qc$celltype <- factor(ASD.sce.qc$celltype, levels = row.names(sort(table(ASD.sce.qc$celltype), 
                                                                             decreasing = T)))

DotPlot(ASD.sce.qc, features = 'Cdkn1a', group.by = 'celltype')
p<-DotPlot(subset(ASD.sce.qc, subset = celltype_main == 'Endothelial'), features = 'Igfbp7', group.by = 'group', scale = F)
pdata <- p$data
# pdata$sex <- sub(".*_", "", pdata$id)
# pdata$group <- sub("_.*", "", pdata$id)
# pdata <- pdata[,c("avg.exp","sex","group", "pct.exp")]

p<-ggplot(pdata, aes(x=id,y=features.plot))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  scale_size_continuous(range = c(3, 7))+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 0.5,vjust=0.5))+
  # scale_color_gradient(low="lightgrey",high="blue")+
  scale_color_viridis()+
  labs(x=NULL,y=NULL)+theme(axis.text.x = element_text(colour = 'black', size = 16),
                            axis.text.y = element_text(colour = 'black', size = 16))

pdf('./figures/SC/main/Igfbp7.pdf', width = 4, height = 1.5)
print(p)
dev.off()

p<-DotPlot(subset(ASD.sce.qc, subset = celltype_main == 'Microgila'), features = 'Cdkn1a', group.by = 'group', scale = F)
pdata <- p$data
pdata$celltype <- sub(".*_", "", pdata$id)
pdata$group <- sub("_.*", "", pdata$id)
pdata <- pdata[,c("avg.exp","celltype","group", "pct.exp")]

p<-ggplot(pdata, aes(x=group,y=celltype))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  scale_size_continuous(range = c(1, 9))+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+
  # scale_color_gradient(low="lightgrey",high="blue")+
  scale_color_viridis()+ ggtitle('Igfbp7')+
  labs(x=NULL,y=NULL)+theme(axis.text.x = element_text(colour = 'black', size = 16),
                            axis.text.y = element_text(colour = 'black', size = 16),
                            axis.title = element_text(colour = 'black', size = 16))

pdf('./figures/SC/main/Igfbp7_celltype1.pdf', width = 4, height = 5)
print(p)
dev.off()

ASD.sce.qc$ident_celltype <- paste(ASD.sce.qc$group, ASD.sce.qc$celltype_main, sep = '_')
p<-DotPlot(ASD.sce.qc, features = 'Igfbp7', group.by = 'ident_celltype', scale = F)+scale_color_viridis()
pdata <- p$data

p<-ggplot(pdata, aes(x=id,y=features.plot))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+coord_flip()+
  scale_size_continuous(range = c(3, 7))+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 0.5,vjust=0.5))+
  # scale_color_gradient(low="lightgrey",high="blue")+
  scale_color_viridis()+
  labs(x=NULL,y=NULL)+theme(axis.text.x = element_text(colour = 'black', size = 16),
                            axis.text.y = element_text(colour = 'black', size = 16))

pdf('sc_igfbp7.pdf', width = 4, height = 5)
print(p)
dev.off()

DotPlot(subset(ASD.sce.qc, subset = celltype_main == 'Microgila'), features = 'Igfbp7', group.by = 'group', scale = F)


Endo <- subset(ASD.sce.qc, subset = celltype_main == 'Endothelial')

VlnPlot(subset(ASD.sce.qc, subset = celltype_main == 'Microgila'), features = "Cdkn1a", pt.size = 0, group.by = 'group',
        cols = pal_npg()(2))+
  stat_compare_means(comparisons = list(c('poly_IC', 'WT'))# ,method = 't.test',
                     )+
  ylim(-0.1, 3)+NoLegend()
# 
VlnPlot(subset(ASD.sce.qc, subset = celltype_main == 'Astrocytes'), features = "Cd74", pt.size = 0, group.by = 'group',
        cols = pal_npg()(2))+
  stat_compare_means(comparisons = list(c('poly_IC', 'WT')))+
  ylim(-0.1, 3)+NoLegend()


ASD.sn.exp <- data.frame(FetchData(object = ASD.sce.qc,
                                   vars = c('Igfbp7', 'Cdkn1a', 'ident', 'celltype_main')))
ASD.sn.hm <- ASD.sn.exp %>% group_by(celltype_main, ident) %>% 
  summarise(mean_exp = mean(Igfbp7, na.rm = T))
ASD.sn.hm <- dcast(ASD.sn.hm,ident~celltype_main,value.var = 'mean_exp')
row.names(ASD.sn.hm) <- ASD.sn.hm$diagnosis
ASD.sn.hm <- ASD.sn.hm[,-1]
ASD.sn.hm <- t(ASD.sn.hm)
p<-pheatmap(ASD.sn.hm,
            show_rownames = T,
            show_colnames = T,
            cluster_cols = F,
            cluster_rows= F,
            fontsize_row = 10,
            fontsize_col = 12,
            cols <- viridis(8),
            border_color = 'grey',
            legend = T,
            angle_col = '90',
            clustering_method = 'single'
) 
pdf('./figures/sn_heatmap_IGFBP7exp.pdf', height = 5, width = 4)
print(p)
dev.off()





DotPlot(subset(ASD.sce.qc, subset = celltype_main == 'ABCs'), features = 'Cdkn1a', group.by = 'group')

# 1DotPlot(ASD.sce.qc, features = 'Cdkn1a', group.by = 'celltype_main')
# save(ASD.sce.qc, file = './object/sc/ASD.sce.qc.RData')

####SASP####
ASD.sce.qc <- AddModuleScore(ASD.sce.qc, features = geneset_RE_SASP, ctrl = 100, name = 'SASP')
p<-VlnPlot(ASD.sce.qc, features = "SASP1", pt.size = 0, group.by = 'group',
        cols = pal_npg()(2))+
  stat_compare_means(comparisons = list(c('poly', 'WT')))+
  ylim(-0.1, 0.22)+NoLegend()
pdf('./figures/SC/Vlnplot_SASP_sc.pdf', width = 6, height = 7.5)
print(p)
dev.off()

p<-VlnPlot(subset(ASD.sce.qc, subset = celltype_main == 'Endothelial'), features = "SASP1", pt.size = 0, group.by = 'group',
        cols = pal_npg()(2))+
  stat_compare_means(comparisons = list(c('poly', 'WT')))+
  ylim(-0.1, 0.22)+NoLegend()
pdf('./figures/SC/Vlnplot_SASP_Endo.pdf', width = 6, height = 7.5)
print(p)
dev.off()


# VlnPlot(subset(ASD.sce.qc, subset = celltype_main == 'Microgila'), features = "SASP1", pt.size = 0, group.by = 'group',
#         cols = pal_npg()(2))+
#   stat_compare_means(comparisons = list(c('poly', 'WT')))+
#   ylim(-0.1, 0.3)+NoLegend()
# 
# VlnPlot(subset(ASD.sce.qc, subset = celltype_main == 'Astrocytes'), features = "SASP1", pt.size = 0, group.by = 'group',
#         cols = pal_npg()(2))+
#   stat_compare_means(comparisons = list(c('poly_IC', 'WT')))+
#   ylim(-0.1, 0.3)+NoLegend()

####Neuron####
ASD.sce.neu <- subset(ASD.sce.qc, subset = celltype_main == 'Neurons')
ASD.sce.neu <- SCTransform(ASD.sce.neu, assay = "RNA", vars.to.regress = "percent.mt")
ASD.sce.neu <- RunPCA(ASD.sce.neu)
ASD.sce.neu <- RunHarmony(ASD.sce.neu, group.by.vars = "orig.ident")
ASD.sce.neu <- RunUMAP(ASD.sce.neu,reduction = "harmony",dims = 1:30, min.dist = 0.3, n.neighbors = 30)
ASD.sce.neu <- FindNeighbors(ASD.sce.neu, reduction = "harmony", dims = 1:30)
ASD.sce.neu <- FindClusters(ASD.sce.neu, 
                             resolution =0.4,
                             method ='igraph',
                             verbose = T)
#0.1 11 0.2 17 0.3 19
DimPlot(ASD.sce.neu, cols = pal_igv()(50), label = T, group.by = 'SCT_snn_res.0.1')

DotPlot(ASD.sce.neu, features = c('Neurod2', 'Slc17a7', 'Slc17a6', 'Slc17a8', 'Camk2a', 'Neurod6', 'Satb2', 'Tbr1',
                                 'Gad2', "Gad1", "Slc32a1", 'Gabrb1', 'Dlx2', 'Pvalb', 'Sst', 'Vip', 'Lhx6'
), group.by = 'seurat_clusters')

DotPlot(ASD.sce.neu, features = c('Camk2a', 'Neurod1', 'Ldb2', 'Slc18a3',
                                  "Gad1", 'Lhfpl3', 'Pcdh15', 'Gria1', 'Vip', 'Gal'
), group.by = 'seurat_clusters')

DotPlot(ASD.sce.neu, features = c('Cdkn1a', 'Igfbp7'), group.by = 'seurat_clusters')
DotPlot(ASD.sce.neu, features = c('Cdkn1a', 'Igfbp7'), group.by = 'group')

Annocelltypes <- c('unknown', 'unknown', 'unknown', 'ExN', 'unknown', 'InN', 'InN', 'ExN', 'InN', 'InN', 
                   'ExN', 'ExN', 'InN', 'InN', 'ExN', 'InN', 'ExN', 'unknown', 'ExN')
Idents(ASD.sce.neu) <- ASD.sce.neu$seurat_clusters
names(Annocelltypes)<- levels(ASD.sce.neu)
ASD.sce.neu <- RenameIdents(ASD.sce.neu, Annocelltypes)
ASD.sce.neu@meta.data$celltype_neu <- ASD.sce.neu@active.ident

Annocelltypes <- c('Neuron_C0', 'Neuron_C1', 'Neuron_C2', 'Neuron_C3', 'Neuron_C4', 'Neuron_C5', 'Neuron_C6', 'Neuron_C7', 'Neuron_C8',
                   'Neuron_C9', 'Neuron_C10', 'Neuron_C11', 'Neuron_C12', 'Neuron_C13', 'Neuron_C14', 
                   'Neuron_C15', 'Neuron_C16')
Idents(ASD.sce.neu) <- ASD.sce.neu$seurat_clusters
names(Annocelltypes)<- levels(ASD.sce.neu)
ASD.sce.neu <- RenameIdents(ASD.sce.neu, Annocelltypes)
ASD.sce.neu@meta.data$cluster <- ASD.sce.neu@active.ident

DimPlot(ASD.sce.neu, cols = pal_igv()(50), group.by = 'celltype_neu')
DimPlot(ASD.sce.neu, cols = pal_igv()(50), group.by = 'cluster')
DimPlot(ASD.sce.neu, cols = pal_igv()(50), group.by = 'celltype')

pdf('./figures/SC/neu/Dimplot_cluster_neuron.pdf', width = 6, height = 5)
p1<-DimPlot(ASD.sce.neu, group.by = 'cluster', cols = col1)+ggtitle('celltype')
print(p1)
dev.off()

#####NEURON_exp#####
ASD.sce.neu$celltype_group <- paste(ASD.sce.neu$cluster, ASD.sce.neu$group, sep = '&')

DotPlot(ASD.sce.neu, features = c('Cdkn1a', 'Igfbp7'), group.by = 'seurat_clusters')

p<-DotPlot(ASD.sce.neu, features = c('Igfbp7'), group.by = 'celltype_group')
pdata <- p$data
pdata$sex <- sub(".*&", "", pdata$id)
pdata$group <- sub("&.*", "", pdata$id)
pdata <- pdata[,c("avg.exp","sex","group", "pct.exp")]
p<-ggplot(pdata, aes(x=sex,y=group))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 0.5,vjust=0.5))+
  # scale_color_gradient(low="lightgrey",high="blue")+
  scale_color_viridis()+
  labs(x=NULL,y=NULL)
p

p<-DotPlot(ASD.sce.neu, features = c('Cdkn1a'), group.by = 'celltype_group')
pdata <- p$data
pdata$sex <- sub(".*&", "", pdata$id)
pdata$group <- sub("&.*", "", pdata$id)
pdata <- pdata[,c("avg.exp","sex","group", "pct.exp")]
p<-ggplot(pdata, aes(x=sex,y=group))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 0.5,vjust=0.5))+
  # scale_color_gradient(low="lightgrey",high="blue")+
  scale_color_viridis()+
  labs(x=NULL,y=NULL)
p

p<-DotPlot(ASD.sce.neu, features = c('Cd74'), group.by = 'celltype_group')
pdata <- p$data
pdata$sex <- sub(".*&", "", pdata$id)
pdata$group <- sub("&.*", "", pdata$id)
pdata <- pdata[,c("avg.exp","sex","group", "pct.exp")]
p<-ggplot(pdata, aes(x=sex,y=group))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 0.5,vjust=0.5))+
   scale_color_viridis()+
  labs(x=NULL,y=NULL)
p







Markers_all_neu <- FindAllMarkers(ASD.sce.neu, assay = "SCT", logfc.threshold = 0.25, only.pos = T)
Markers_all_neu_top10 <- Markers_all_neu %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Markers_all_neu_top10, file = './results/SC/Markers_all_neu_top10.csv')
save(Markers_all_neu, file = './object/sc/Markers_all_neu.RData')
save(Markers_neu, file = './object/sc/Markers_neu.RData')


Markers_all_neu_top10[which(Markers_all_neu_top10$cluster == '15'),]$gene
#up 12 
#down 10 15 EX
#celltype marker enrich

Markers_neu <- FindAllMarkers(ASD.sce.neu, assay = "SCT", logfc.threshold = 0, only.pos = F)

genelist <- Markers_all_neu[which(Markers_all_neu$cluster == '12'),]$gene
enrich.go.12 <- enrichGO(gene = genelist,
                      OrgDb = 'org.Mm.eg.db',
                      keyType = 'SYMBOL',
                      ont = 'ALL',
                      pAdjustMethod = 'fdr',
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
enrich.go.df.12 <- as.data.frame(enrich.go.12)




geneset_RE = msigdbr(species = "Mus musculus",#Homo sapiens
                     category = "C2",
                     subcategory = "CP:REACTOME"
) %>% dplyr::select(gs_name,gene_symbol)
geneset_RE <- geneset_RE %>% split(x = .$gene_symbol, f = .$gs_name)

geneset_GOBP <- msigdbr(species = "Mus musculus",
                        category = "C5", #Oncology
                        subcategory = "BP"#("BP" "CC" "MF")
) %>% dplyr::select(gs_name,gene_symbol)
geneset_GOBP <- geneset_GOBP %>% split(x = .$gene_symbol, f = .$gs_name)

library(fgsea)
library(msigdbr)
genelist <- Markers_neu[which(Markers_neu$cluster == 'Neuron_C12'),]
rank <- genelist$avg_log2FC[order(genelist$avg_log2FC, decreasing = T)]
names(rank) <- genelist$gene[order(genelist$avg_log2FC, decreasing = T)]
GSEA_res_RE_12 <- fgsea(pathways = geneset_RE, stats = rank, minSize=5, maxSize=500, nperm=10000)
GSEA_df_RE_12<-as.data.frame(GSEA_res_RE_12)
GSEA_df_RE_12$leadingEdge<-as.character(GSEA_df_RE_12$leadingEdge)
write.csv(GSEA_df_RE_12, file = './results/SC/neu/GSEA_res_RE_neu_12.csv')

genelist <- Markers_neu[which(Markers_neu$cluster == 'Neuron_C10'),]
rank <- genelist$avg_log2FC[order(genelist$avg_log2FC, decreasing = T)]
names(rank) <- genelist$gene[order(genelist$avg_log2FC, decreasing = T)]
GSEA_res_RE_10 <- fgsea(pathways = geneset_RE, stats = rank, minSize=5, maxSize=500, nperm=10000)
GSEA_df_RE_10<-as.data.frame(GSEA_res_RE_10)
GSEA_df_RE_10$leadingEdge<-as.character(GSEA_df_RE_10$leadingEdge)
write.csv(GSEA_df_RE_10, file = './results/SC/neu/GSEA_res_RE_neu_10.csv')

genelist <- Markers_neu[which(Markers_neu$cluster == 'Neuron_C15'),]
rank <- genelist$avg_log2FC[order(genelist$avg_log2FC, decreasing = T)]
names(rank) <- genelist$gene[order(genelist$avg_log2FC, decreasing = T)]
GSEA_res_RE_15 <- fgsea(pathways = geneset_RE, stats = rank, minSize=5, maxSize=500, nperm=10000)
GSEA_df_RE_15<-as.data.frame(GSEA_res_RE_15)
GSEA_df_RE_15$leadingEdge<-as.character(GSEA_df_RE_15$leadingEdge)
write.csv(GSEA_df_RE_15, file = './results/SC/neu/GSEA_res_RE_neu_15.csv')

result_df <- data.frame()
result_df <- rbind(
  result_df,
  data.frame(Value = GSEA_df_RE_10 %>% top_n(n = 3, wt = NES), Source = "Neuron10"),
  data.frame(Value = GSEA_df_RE_12 %>% top_n(n = 4, wt = NES), Source = "Neuron12"),
  data.frame(Value = GSEA_df_RE_15 %>% top_n(n = 3, wt = NES), Source = "Neuron15")
)
colnames(result_df) <- c(colnames(GSEA_df_RE_15), 'cluster')

result_df <- result_df[-7,]

result_df$labelx <- rep(0.01, nrow(result_df))
result_df$labely <- seq(nrow(result_df), 1)

result_df <- result_df[order(result_df$cluster),]
result_df$pathway <- factor(result_df$pathway, levels = rev(result_df$pathway[order(result_df$cluster)]))

p<-ggplot(data = result_df, aes(NES, pathway, fill = cluster)) +
  # geom_col(aes(x=NES,y=reorder(pathway,NES),fill=cluster),width = 0.8) +
  geom_bar(stat="identity",alpha=0.8,width = 0.8) + 
  geom_text(aes(x=labelx,y=labely,label = pathway),size=5, hjust =0)+
  theme_classic()+
  #scale_fill_distiller(palette="Oranges", direction = 1)+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 16),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.line.y = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 14),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 16))+
  xlab("Count")+ylab('GO Terms')+
  scale_x_continuous(expand = c(0.01,0))+
  scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF", "#00A087FF"))

pdf('./figures/SC/neu/enrich_RE_barplot.pdf', width = 8, height = 6)
print(p)
dev.off()




Markers_all_neu_top5 <- Markers_all_neu %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
neu_genes <- c(Markers_all_neu_top5[which(Markers_all_neu_top5$cluster == '10'),]$gene,
               Markers_all_neu_top5[which(Markers_all_neu_top5$cluster == '12'),]$gene,
               Markers_all_neu_top5[which(Markers_all_neu_top5$cluster == '15'),]$gene)
ASD.sce.neu_s <- subset(ASD.sce.neu, seurat_clusters %in% c(10,12,15))

pdf('./figures/SC/neu/Dotplot_markers_OR.pdf', width = 10, height = 5)
p3<-DotPlot(ASD.sce.neu_s, features = neu_genes, group.by = 'cluster', cols = c('lightyellow', 'red3'))+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 18,color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 18,color = "black"),
        axis.title = element_text(size = 16,color = 'black'),
        legend.text = element_text(size = 16,color = "black"))+labs(x="",y="")
print(p3)
dev.off()





ASD.sce.neu$group <- ASD.sce.neu$ident
dittoBarPlot(ASD.sce.neu, 'celltype_neu', group.by = 'group', main = '', color.panel = pal_npg()(10))
dittoBarPlot(ASD.sce.neu, 'seurat_clusters', group.by = 'group', main = '', color.panel = pal_igv()(50))

#####OR analysis##### 
library(sscVis)
library(data.table)
library(grid)
library(cowplot)
library(ggrepel)
library(readr)
library(plyr)
library(tidyr)
library(ggpubr)
meta <- ASD.sce.neu@meta.data
meta$loc <- meta$group
meta$meta.cluster <- meta$SCT_snn_res.0.4  
Neu.OR <- do.tissueDist(cellInfo.tb = meta,#这里的输入文件需要的是整理好的有分组和细胞类型的metadata文件
                         out.prefix = "./figures/SC/neu/ORplot_NEU_Clusters111",#设定导出的文件名，自己设置
                         pdf.width = 9,#热图的宽设置
                         pdf.height = 5,#热图的高度设置
                         verbose=1)
save(Neu.OR, file = './object/sc/Neu.OR.RData')

DEG_neu <- FindMarkers(object = ASD.sce.neu, ident.1 = 'poly', ident.2 = 'WT', min.pct = 0.1,
                          logfc.threshold = 0, group.by = 'group')
DEG_neu$gene <- row.names(DEG_neu)


# save(ASD.sce.neu, file = './object/sc/ASD.sce.neu.RData')

#####monocol_neuron_s#####
library(monocle)
data <- as(as.matrix(ASD.sce.neu@assays$RNA@counts), 'sparseMatrix') 
pd <- new('AnnotatedDataFrame', data = ASD.sce.neu@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#create CellDataSet
ASD_mono_neu <- newCellDataSet(data, phenoData = pd, featureData = fd,
                                 lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
rm(data, pd, fData, fd)

#Estimate size factors and dispersions
ASD_mono_neu <- estimateSizeFactors(ASD_mono_neu)
ASD_mono_neu <- estimateDispersions(ASD_mono_neu)
#QC
ASD_mono_neu <- detectGenes(ASD_mono_neu, min_expr = 0.1)
print(head(fData(ASD_mono_neu)))
print(head(pData(ASD_mono_neu)))
expressed_genes <- row.names(subset(fData(ASD_mono_neu), num_cells_expressed >= 10))
pData(ASD_mono_neu)$Total_mRNAs <- Matrix::colSums(exprs(ASD_mono_neu))
ASD_mono_neu <- ASD_mono_neu[,pData(ASD_mono_neu)$Total_mRNAs < 1e6]

# cds_DGT <- ASD_mono_neu
# diff_test_res <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~celltype_main")
# ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
# # ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:2000]#或者选择前多少的基因
# cds_DGT <- setOrderingFilter(cds_DGT, ordering_genes)
# plot_ordering_genes(cds_DGT)

#使用seurat计算
# Mdk.sce.tumor_s <- NormalizeData(Mdk.sce.tumor_s)
# seurat_variable_genes <- VariableFeatures(FindVariableFeatures(Mdk.sce.tumor_s, assay = "RNA"), assay = "RNA")
cds_seurat <- ASD_mono_neu
seurat_variable_genes <- ASD.sce.neu@assays$SCT@var.features
cds_seurat <- setOrderingFilter(cds_seurat, seurat_variable_genes)
plot_ordering_genes(cds_seurat)

#save(cds_seurat, file = 'cds_seurat.RData')
# 
# cds_monocle <- ASD_mono_neu
# disp_table <- dispersionTable(cds_monocle)
# disp_table <- subset(disp_table, mean_expression >= 0.1)
# head(disp_table)
# cds_monocle <- setOrderingFilter(cds_monocle, disp_table$gene_id)
# plot_ordering_genes(cds_monocle)

ASD_cds_neu <- cds_seurat
ASD_cds_neu <- reduceDimension(ASD_cds_neu, max_components = 2,reduction_method = 'DDRTree',
                                 residualModelFormulaStr = '~orig.ident')
ASD_cds_neu <- orderCells(ASD_cds_neu)
plot_cell_trajectory(ASD_cds_neu, color_by = "Pseudotime")




ASD.sce.non <- subset(ASD.sce.qc, subset = celltype_main != 'Neurons')
ASD.sce.non@meta.data$celltype_main <- gsub('Microglia', 'Microgila', ASD.sce.non@meta.data$celltype_main)
ASD.sce.non@meta.data$celltype <- gsub('Microglia', 'Microgila', ASD.sce.non@meta.data$celltype)

ASD.sce.non <- SCTransform(ASD.sce.non, assay = "RNA", vars.to.regress = "percent.mt")
ASD.sce.non <- RunPCA(ASD.sce.non)
ASD.sce.non <- RunHarmony(ASD.sce.non, group.by.vars = "orig.ident")
ASD.sce.non <- RunUMAP(ASD.sce.non,reduction = "harmony",dims = 1:30, min.dist = 0.4, n.neighbors = 30)
ASD.sce.non <- FindNeighbors(ASD.sce.non, reduction = "harmony", dims = 1:30)
ASD.sce.non <- FindClusters(ASD.sce.non,
                            resolution =0.1,
                            method ='igraph',
                            verbose = T)
#0.1 13 0.2 16 0.3 19 0.4 23
DimPlot(ASD.sce.non, cols = pal_igv()(50))

DimPlot(ASD.sce.non, cols = pal_igv()(50), group.by = 'celltype')

DotPlot(ASD.sce.non, features = c("Hexb", "C1qvb", "C1qc", "Ctss", "Tmem119", "Itgam","Cx3cr1",#Microgila 
                                 "Gja1", "Aldoc", "Plpp3", "Slc1a3", "Slc1a2", "Acsl6",  #Astrocyte
                                 "Cldn11", "Mog", "Mbp", "Mag", "Ermn", "Plp1", "Enpp2", #OLG
                                 "Pdgfra", "Tnr", "Epn2", #OPC
                                 "Ccdc153", "Tmem212", "Cfap54", "Cfap44", #OP
                                 'Slc6a13', #VLMC
                                 'Slc47a1', #ABC
                                 'Ttr', #CPC
                                 "Cldn5", "Igfbp7", "Ptprb", 'Flt1'  #endothelial
), group.by = 'seurat_clusters')#,group.by = 'predicted_MCA'

ASD.sce.non$group <- ASD.sce.non$ident
dittoBarPlot(ASD.sce.non, 'celltype', group.by = 'group', main = '', color.panel = pal_igv()(50))

save(ASD.sce.non, file = './object/sc/ASD.sce.non.RData')


Idents(ASD.sce.non) <- ASD.sce.non$celltype
markers_MG <- FindMarkers(object = ASD.sce.qc, ident.1 = 'poly_IC', ident.2 = 'WT', min.pct = 0.1,
                       logfc.threshold = 0.25, group.by = 'group', subset.ident = 'Microgila')
markers_AST <- FindMarkers(object = ASD.sce.qc, ident.1 = 'poly_IC', ident.2 = 'WT', min.pct = 0.1,
                          logfc.threshold = 0.25, group.by = 'group', subset.ident = 'Astrocytes')


intersect(row.names(markers_MG), wgcna_gene$gene)

####Endo####
ASD.sce.Endo <- subset(ASD.sce.qc, celltype_main == 'Endothelial')
igfbp7_counts <- GetAssayData(ASD.sce.Endo, assay = "RNA", slot = "counts")["Igfbp7", ]
ASD.sce.Endo$IGF_Status <- ifelse(igfbp7_counts > 0, "IGFpos", "IGFneg")
table(ASD.sce.Endo$IGF_Status)

Idents(ASD.sce.Endo) <- ASD.sce.Endo$IGF_Status
markers_IGF <- FindMarkers(object = ASD.sce.Endo, ident.1 = 'IGFpos', ident.2 = 'IGFneg', min.pct = 0.1,
                       logfc.threshold = 0.25)
markers_IGF <- markers_IGF[which(markers_IGF$p_val < 0.05),]
write.csv(markers_IGF, file = './results/sc/markers_IGF_ENDO.csv')


##
sc_marker_list <- list()
for (i in unique(ASD.sce.qc$celltype_main)) {
  markers <- FindMarkers(object = ASD.sce.qc, ident.1 = 'poly_IC', ident.2 = 'WT', min.pct = 0.1,
                         logfc.threshold = 0.1, group.by = 'group', subset.ident = i)
  markers <- subset(markers, p_val < 0.05)
  sc_marker_list[[i]] <- markers
  rm(markers)
}

library(msigdbr)
library(fgsea)
geneset_KEGG = msigdbr(species = "Mus musculus",#Homo sapiens
                       category = "C2", 
                       subcategory = "KEGG"
) %>% dplyr::select(gs_name,gene_symbol)
geneset_KEGG <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)

geneset_RE = msigdbr(species = "Mus musculus",#Homo sapiens
                     category = "C2",
                     subcategory = "CP:REACTOME"
) %>% dplyr::select(gs_name,gene_symbol)
geneset_RE <- geneset_RE %>% split(x = .$gene_symbol, f = .$gs_name)


Run_fGSEA <- function(DEG, geneset){
  rank <- DEG$avg_log2FC[order(DEG$avg_log2FC, decreasing = T)]
  DEG$gene <- row.names(DEG)
  names(rank) <- DEG$gene[order(DEG$avg_log2FC, decreasing = T)]
  GSEA_res <- fgsea(pathways = geneset, stats = rank, minSize=5, maxSize=500, nperm=10000)
  GSEA_df <- as.data.frame(GSEA_res)
  GSEA_df$leadingEdge <- as.character(GSEA_res$leadingEdge)
  return(GSEA_df)
}

GSEA_Neu_KEGG <- Run_fGSEA(sc_marker_list$Neurons, geneset_KEGG)
GSEA_Neu_RE <- Run_fGSEA(sc_marker_list$Neurons, geneset_RE)
p<-GSEA_BARPLOT(GSEA_Neu_RE)
write.csv(GSEA_Neu_RE, file = './results/DEG/enrich/GSEA_Neu_RE.csv')
pdf("./figures/SC/GSEA/GSEA_Neu_RE_barplot.pdf", width = 15, height = 3)
print(p)
dev.off()


GSEA_CPC_KEGG <- Run_fGSEA(sc_marker_list$CPCs, geneset_KEGG)
GSEA_CPC_RE <- Run_fGSEA(sc_marker_list$CPCs, geneset_RE)
p<-GSEA_BARPLOT(GSEA_CPC_RE)
write.csv(GSEA_CPC_RE, file = './results/DEG/enrich/GSEA_CPC_RE.csv')
pdf("./figures/SC/GSEA/GSEA_CPC_RE_barplot.pdf", width = 15, height = 3)
print(p)
dev.off()


GSEA_VLMC_KEGG <- Run_fGSEA(sc_marker_list$VLMCs, geneset_KEGG)
GSEA_VLMC_RE <- Run_fGSEA(sc_marker_list$VLMCs, geneset_RE)
p<-GSEA_BARPLOT(GSEA_VLMC_RE)
write.csv(GSEA_VLMC_RE, file = './results/DEG/enrich/GSEA_VLMC_RE.csv')
pdf("./figures/SC/GSEA/GSEA_VLMC_RE_barplot.pdf", width = 8, height = 3)
print(p)
dev.off()

sc_marker_num <- data.frame()
for (i in 1:length(sc_marker_list)) {
  sc_marker_num[i,1] <- dim(sc_marker_list[[i]])[1]
}
row.names(sc_marker_num) <- unique(ASD.sce.qc$celltype_main)
colnames(sc_marker_num) <- 'count'
sc_marker_num$celltype <- unique(ASD.sce.qc$celltype_main)

# unique(ASD.sce.qc$celltype_main)
p<-ggplot(sc_marker_num,aes(celltype,count,fill=celltype))+
  geom_bar(stat="summary",fun=mean,position="dodge", width = 0.6)+
  theme_bw()+scale_fill_manual(values = col1)+xlab('')+
  theme(legend.title=element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size = 13, color = 'black',angle=45,hjust = 1,vjust=1),
        axis.text.y = element_text(size = 14, color = 'black'))+
  geom_text(aes(label = count), vjust = - 0.1, size = 5)
pdf("./figures/SC/Barplot_SC_DEG_NUM.pdf", width = 10, height = 5)
print(p)
dev.off()

####cellchat####
library(CellChat)
library(patchwork)
library(future)
set.seed(123)
showDatabaseCategory(CellChatDB.mouse)
future::plan("multisession", workers = 24) # do parallel


ASD.sce.wt <- subset(ASD.sce.qc, ident == 'WT')
cellchat_wt <- createCellChat(object = ASD.sce.wt, meta = ASD.sce.wt@meta.data, group.by = "celltype_main")
cellchat_wt@DB <- CellChatDB.mouse
cellchat_wt <- subsetData(cellchat_wt)

future::plan("multisession", workers = 24)
cellchat_wt <- identifyOverExpressedGenes(cellchat_wt)
cellchat_wt <- identifyOverExpressedInteractions(cellchat_wt)
cellchat_wt <- projectData(cellchat_wt, PPI.mouse)

cellchat_wt <- computeCommunProb(cellchat_wt, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_wt <- filterCommunication(cellchat_wt, min.cells = 10)
cellchat_wt <- computeCommunProbPathway(cellchat_wt)
cellchat_wt <- aggregateNet(cellchat_wt)
save(cellchat_wt, file = './object/sc/cellchat_wt.RData')
# df.netp <- subsetCommunication(cellchat_wt)
# write.csv(df.netp, "./results/cellchat/ctl_net_pathway_new.csv", row.names = F)


ASD.sce.poly <- subset(ASD.sce.qc, ident == 'poly')
cellchat_poly <- createCellChat(object = ASD.sce.poly, meta = ASD.sce.poly@meta.data, group.by = "celltype_main")
cellchat_poly@DB <- CellChatDB.mouse
cellchat_poly <- subsetData(cellchat_poly)

future::plan("multisession", workers = 24)
cellchat_poly <- identifyOverExpressedGenes(cellchat_poly)
cellchat_poly <- identifyOverExpressedInteractions(cellchat_poly)
cellchat_poly <- projectData(cellchat_poly, PPI.mouse)

cellchat_poly <- computeCommunProb(cellchat_poly, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_poly <- filterCommunication(cellchat_poly, min.cells = 10)
cellchat_poly <- computeCommunProbPathway(cellchat_poly)
cellchat_poly <- aggregateNet(cellchat_poly)
save(cellchat_poly, file = './object/sc/cellchat_poly.RData')
df.netp <- subsetCommunication(cellchat_poly)

object.list <- list()
object.list[[1]] <- cellchat_wt 
object.list[[2]] <- cellchat_poly
names(object.list) <- c('wt', 'poly')
cellchat.All <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
save(cellchat.All, file = './object/sc/cellchat.All.RData')
p1<-compareInteractions(cellchat.All, show.legend = F, group = c(1,2),
                        color.use = c("#E64B35", "#4DBBD5"))
p2<-compareInteractions(cellchat.All, show.legend = F, group = c(1,2), measure = "weight",
                        color.use = c("#E64B35", "#4DBBD5"))


groupSize1 <- as.numeric(table(cellchat_wt@idents))
groupSize2 <- as.numeric(table(cellchat_poly@idents))
netVisual_circle(cellchat_wt@net$count, vertex.weight = groupSize1, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions in ctl")

netVisual_circle(cellchat_poly@net$count, vertex.weight = groupSize2, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions in shMDK")

cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt)
cellchat_poly <- netAnalysis_computeCentrality(cellchat_poly)
pathway.union <- union(cellchat_wt@netP$pathways, cellchat_poly@netP$pathways)

netAnalysis_signalingRole_heatmap(cellchat_wt, pattern = "outgoing", 
                                  signaling = pathway.union, title = 'wt', 
                                  width = 12, height = 20)
netAnalysis_signalingRole_heatmap(cellchat_poly, pattern = "outgoing", 
                                  signaling = pathway.union, title = 'ployIC', 
                                  width = 12, height = 20)

rankNet(cellchat.All, mode = "comparison", comparison = c(1,2), stacked = T, do.stat = TRUE,targets.use = 'Endothelial',
        color.use = c("#E64B35", "#4DBBD5"), do.flip = T, measure = 'count')




####datasets####
library(Seurat)
library(readr)
library(stringr)


#####MIA——FATAL#####
folder_path <- "./data/GEO/"
tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE)
tsv_names <- substring(tsv_files, 12)
tsv_names <- str_remove(tsv_names, "\\..*")

# 初始化一个列表来存储 Seurat 对象
seurat_objects <- list()

for (file in tsv_files) {
  counts <- read.delim(file, row.names = 1)  
  counts <- as.data.frame(t(counts))
  seurat_obj <- CreateSeuratObject(counts = counts)
  tsv_names <- substring(file, 12)
  tsv_names <- str_remove(tsv_names, "\\..*")
  seurat_obj$orig.ident <- tsv_names
  seurat_objects[[tsv_names]] <- seurat_obj
  print(paste(tsv_names, 'over', sep = ' '))
}
save(seurat_objects, file = './object/FATAL_seurat_objects.RData')

merged_object <- merge(
  x = seurat_objects[[1]],  # 第一个对象作为基础
  y = seurat_objects[-1],   # 剩余的对象列表
  add.cell.ids = tsv_names  # 添加细胞标识符
  # merge.data = TRUE,  # 合并数据槽
  # project = "SeuratProject"  # 项目名称
)
rm(seurat_objects)

merged_object[["percent.mt"]] <- PercentageFeatureSet(merged_object, pattern = "^mt-")
VlnPlot(merged_object, features = c('nCount_RNA', 'nFeature_RNA'), pt.size = 0.1) + NoLegend()
VlnPlot(merged_object, features = "percent.mt", pt.size = 0.1) + NoLegend()



####
library(scDblFinder)
library(BiocParallel)

merged_object.sce <- as.SingleCellExperiment(merged_object)
merged_object.sce <- scDblFinder(merged_object.sce, samples="orig.ident", BPPARAM=MulticoreParam(12))
table(merged_object.sce$scDblFinder.class)

merged_object <- as.Seurat(merged_object.sce)
table(merged_object$orig.ident)
table(merged_object$scDblFinder.class)
merged_object <- subset(merged_object, scDblFinder.class == 'singlet')


MIA_Fatal <- subset(merged_object, subset = nFeature_RNA > 800 & nFeature_RNA < 6000 & percent.mt < 5 & 
                      nCount_RNA > 1000 & nCount_RNA < 20000)
rm(merged_object)
# save(MIA_Fatal, file = './object/MIA_Fatal.RData')

MIA_Fatal <- SCTransform(MIA_Fatal, assay = "RNA", vars.to.regress = c("percent.mt"), 
                          ncells = round(dim(MIA_Fatal)[2]/10))
# MIA_Fatal <- SCTransform(MIA_Fatal, assay = "RNA", vars.to.regress = "percent.mt")
MIA_Fatal <- RunPCA(MIA_Fatal)
MIA_Fatal <- RunHarmony(MIA_Fatal, group.by.vars = "orig.ident")
MIA_Fatal <- RunUMAP(MIA_Fatal,reduction = "harmony",dims = 1:30, min.dist = 0.4)#
# MIA_Fatal <- RunTSNE(MIA_Fatal,reduction = "harmony",dims = 1:30)
MIA_Fatal <- FindNeighbors(MIA_Fatal, reduction = "harmony", dims = 1:30)
MIA_Fatal <- FindClusters(MIA_Fatal, 
                           resolution =0.5,
                           method ='igraph',
                           verbose = T)
# 0.1 14 0.2 21 0.3 25 0.4 26 0.5 27 
MIA_Fatal$seurat_clusters <- MIA_Fatal$SCT_snn_res.0.4
DimPlot(MIA_Fatal, group.by = 'SCT_snn_res.0.4', label = T)

DotPlot(MIA_Fatal, features = c('Cck', 'Ebf3', 'Tbr1', 'Neurod6', 'Emoes', 'Dlx2', 'Gad2', 'Isl1', 'Vim', 'Pax6', 'Ly86',
                                'Cldn5', 'Alas2', 'Tal1'
), group.by = 'seurat_clusters')

DotPlot(MIA_Fatal, features = c('Ebf3', 'Tbr1', 'Tiam2', 'Sema3c', 'Plxna4', 'Syt4', 'Fezf1', 'Calb2', 'Zic1', 'Cck', 'Pbx3',
                                'Eomes', 'Hes5', 'Unc5d', 'Neurog2', 'Tk1', 'Top2a', 'Cenpa'
), group.by = 'seurat_clusters')

DotPlot(MIA_Fatal, features = c('Top2a', 'Rrm2', 'Maf', 'Gria1', 'Six3', 'Syt6', 'Sst', 'Ndn', 'Ccnd2', 'Npy', 'Nr2f2',
                                'Nkx2-1', 'Zic4', 'Otp', 'Zfp503', 'Sybu'
), group.by = 'seurat_clusters')

DotPlot(MIA_Fatal, features = c('Vim', 'Pax6', 'Aldoc', 'Tnc', 'Mest', 'Ptprz1', 'Top2a', 'Tk1', 'Mcm6', 'Elavl4', 'Rspo3',
                                'Wnt8b'
), group.by = 'seurat_clusters')

markers_9 <- FindMarkers(MIA_Fatal, ident.1 = '9', logfc.threshold = 0.25, only.pos = T)
markers_22 <- FindMarkers(MIA_Fatal, ident.1 = '22', logfc.threshold = 0.25, only.pos = T)
markers_23 <- FindMarkers(MIA_Fatal, ident.1 = '23', logfc.threshold = 0.25, only.pos = T)

DotPlot(MIA_Fatal, features = c('Igfbp7', 'Cdkn1a', 'Cd74'
), group.by = 'seurat_clusters')



 library(genesorteR)
library(dplyr)
gs <- sortGenes(MIA_Fatal@assays$SCT@data, MIA_Fatal$seurat_clusters)
# gs1 <- getMarkers(gs, quant = 0.99)
gs_spec <- gs$specScore
genes <- rownames(gs_spec)
gs_top_genes <- data.frame()
for (cluster in colnames(gs_spec)) {
  cluster_scores <- data.frame(
    gene = genes,
    specScore = gs_spec[, cluster],
    cluster = cluster
  )
  top_genes <- cluster_scores %>%
    arrange(desc(specScore)) %>%
    head(10)
  gs_top_genes <- rbind(gs_top_genes, top_genes)
}

# DotPlot(ASD.sce.qc, features = unique(gs_top_genes$gene), group.by = 'seurat_clusters')
# save(gs_top_genes, file = './object/sc/gs_topgenes.RData')

MIA_Fatal_Endo <- subset(MIA_Fatal, seurat_clusters == '19')

MIA_Fatal$time <- str_sub(MIA_Fatal$orig.ident, -3)
MIA_Fatal$gender <- sub(".*(.{1}).{3}$", "\\1", MIA_Fatal$orig.ident)
MIA_Fatal$group <- substr(MIA_Fatal$orig.ident, 1, 3)
MIA_Fatal$group_gender <- paste(MIA_Fatal$group, MIA_Fatal$gender, sep = '_')

MIA_Fatal_Endo <- subset(MIA_Fatal, seurat_clusters == '19')
subset(MIA_Fatal_Endo, time == 'E14')
VlnPlot(subset(MIA_Fatal_Endo, time == 'E14'), features = "Igfbp7", pt.size = 0, group.by = 'group_gender',
        cols = pal_npg()(4))+
  stat_compare_means(comparisons = list(c('MIA_F', 'PBS_F'), c('MIA_M', 'PBS_M')))+
  ylim(-0.1, 10)+NoLegend()

VlnPlot(subset(MIA_Fatal_Endo, time == 'E18'), features = "Igfbp7", pt.size = 0, group.by = 'group_gender',
        cols = pal_npg()(4))+
  stat_compare_means(comparisons = list(c('MIA', 'PBS')))+
  ylim(-0.1, 10)+NoLegend()

VlnPlot(subset(MIA_Fatal_Endo, time == 'E18'), features = "Cdkn1a", pt.size = 0, group.by = 'group_gender',
        cols = pal_npg()(4))+
  stat_compare_means(comparisons = list(c('MIA_F', 'PBS_F'), c('MIA_M', 'PBS_M')))+
  ylim(-0.1, 10)+NoLegend()

MIA_Fatal_MG <- subset(MIA_Fatal, seurat_clusters == '25')

VlnPlot(subset(MIA_Fatal_MG, time == 'E14'), features = "Igfbp7", pt.size = 0, group.by = 'group',
        cols = pal_npg()(2))+
  stat_compare_means(comparisons = list(c('MIA', 'PBS')))+
  ylim(-0.1, 10)+NoLegend()

MIA_Fatal_23 <- subset(MIA_Fatal, seurat_clusters == '23')

VlnPlot(subset(MIA_Fatal_23, time == 'E14'), features = "Igfbp7", pt.size = 0, group.by = 'group_gender',
        cols = pal_npg()(4))+
  # stat_compare_means(comparisons = list(c('MIA', 'PBS')))+
  ylim(-0.1, 10)+NoLegend()


MIA_Fatal$batch <- substr(MIA_Fatal$orig.ident, 4, 4)
MIA_Fatal$batch <- gsub('F', '1', MIA_Fatal$batch)
MIA_Fatal$batch <- gsub('M', '1', MIA_Fatal$batch)


vln_group <- function(ti, clu, ge, ba){
  p=VlnPlot(subset(MIA_Fatal, time == ti & seurat_clusters == clu & batch == ba), features = ge, pt.size = 0,
            group.by = 'group_gender', cols = pal_npg()(4))+
    stat_compare_means(comparisons = list(c('MIA_F', 'PBS_F'), c('MIA_M', 'PBS_M')))+
    ylim(-0.1, 10)+NoLegend()
  return(p)
}

dittoBarPlot(MIA_Fatal, 'group_gender', group.by = 'batch')

vln_group('E14', '19', 'Igfbp7', 'C')
vln_group('E14', '24', 'Cdkn1a', 'C')


###
#MG 24
#AST 
#OLG 
#OPC 
#Endo 19

MIA_Fatal <- AddModuleScore(MIA_Fatal, features = SASP, ctrl = 100, name = 'SASP')
p<-VlnPlot(MIA_Fatal, features = "SASP1", pt.size = 0, group.by = 'group_gender',sort = T,
           cols = pal_npg()(4))+
  # stat_compare_means(comparisons = list(c('MIA', 'PBS')))+
  ylim(-0.1, 0.8)+NoLegend()

VlnPlot(subset(MIA_Fatal,  batch == 'C'), features = "SASP1", pt.size = 0,sort = T,
        group.by = 'group_gender', cols = pal_npg()(4))+
  stat_compare_means(comparisons = list(c('MIA_F', 'PBS_F'), c('MIA_M', 'PBS_M')))+
  ylim(-0.1, 10)+NoLegend()

pdf('./figures/SC/Vlnplot_SASP_sc.pdf', width = 6, height = 7.5)
print(p)
dev.off()













MIA_Fatal$time <- str_sub(MIA_Fatal$orig.ident, -3)
MIA_Fatal$gender <- sub(".*(.{1}).{3}$", "\\1", MIA_Fatal$orig.ident)
MIA_Fatal$group <- substr(MIA_Fatal$orig.ident, 1, 3)
MIA_Fatal$group_gender <- paste(MIA_Fatal$group, MIA_Fatal$gender, sep = '_')


DotPlot(MIA_Fatal, features = 'Cdkn1a', group.by = 'orig.ident')

save(MIA_Fatal, file = './object/MIA_Fatal.RData')

MIA_Fatal_C <- subset(MIA_Fatal, subset = orig.ident %in% c('MIACFE14', 'MIACME14', 'PBSCFE14', 'PBSCME14'))
p<-DotPlot(MIA_Fatal_C, features = 'Igfbp7', group.by = 'group_gender')
# DotPlot(MIA_Fatal_C, features = 'Cdkn1a', group.by = 'group_gender')
pdata <- p$data
pdata$sex <- sub(".*_", "", pdata$id)
pdata$group <- sub("_.*", "", pdata$id)
pdata <- pdata[,c("avg.exp","sex","group", "pct.exp")]

p<-ggplot(pdata, aes(x=sex,y=group))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 0.5,vjust=0.5))+
  # scale_color_gradient(low="lightgrey",high="blue")+
  scale_color_viridis()+
  labs(x=NULL,y=NULL)
pdf('./figures/SC/MIA/Igfbp7.pdf', width = 3, height = 1.5)
print(p)
dev.off()

pdf('./figures/SC/MIA/Igfbp7_legend.pdf', width = 5, height = 3.5)
print(p)
dev.off()


MIA_Fatal_2 <- subset(MIA_Fatal, subset = orig.ident %in% c('MIA2FE18', 'MIA2ME18', 'PBS2FE18', 'PBS2ME18'))
# DotPlot(MIA_Fatal_2, features = 'Igfbp7', group.by = 'group_gender')

p<-DotPlot(MIA_Fatal_2, features = 'Cdkn1a', group.by = 'group_gender')
# DotPlot(MIA_Fatal_C, features = 'Cdkn1a', group.by = 'group_gender')
pdata <- p$data
pdata$sex <- sub(".*_", "", pdata$id)
pdata$group <- sub("_.*", "", pdata$id)
pdata <- pdata[,c("avg.exp","sex","group", "pct.exp")]

p<-ggplot(pdata, aes(x=sex,y=group))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 0.5,vjust=0.5))+
  # scale_color_gradient(low="lightgrey",high="blue")+
  scale_color_viridis()+
  labs(x=NULL,y=NULL)
pdf('./figures/SC/MIA/cdkn1a.pdf', width = 3, height = 1.5)
print(p)
dev.off()
pdf('./figures/SC/MIA/cdkn1a_legend.pdf', width = 5, height = 3.5)
print(p)
dev.off()


ASD.sce.qc$celltype_group <- paste(ASD.sce.qc$celltype_main, ASD.sce.qc$group, sep = '&')

p<-DotPlot(ASD.sce.qc, features = c('App'), group.by = 'celltype_group')
pdata <- p$data
pdata$sex <- sub(".*&", "", pdata$id)
pdata$group <- sub("&.*", "", pdata$id)
pdata <- pdata[,c("avg.exp","sex","group", "pct.exp")]
p<-ggplot(pdata, aes(x=sex,y=group))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 0.5,vjust=0.5))+
  scale_color_viridis()+
  labs(x=NULL,y=NULL)
p


p<-DotPlot(ASD.sce.qc, features = c('Cdkn1a'), group.by = 'celltype_group')
pdata <- p$data
pdata$sex <- sub(".*&", "", pdata$id)
pdata$group <- sub("&.*", "", pdata$id)
pdata <- pdata[,c("avg.exp","sex","group", "pct.exp")]
p<-ggplot(pdata, aes(x=sex,y=group))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=0,hjust = 0.5,vjust=0.5))+
  scale_color_viridis()+
  labs(x=NULL,y=NULL)
p

p<-DotPlot(ASD.sce.qc, features = c('Cd74'), group.by = 'celltype_group')
pdata <- p$data
pdata$sex <- sub(".*&", "", pdata$id)
pdata$group <- sub("&.*", "", pdata$id)
pdata <- pdata[,c("avg.exp","sex","group", "pct.exp")]
p<-ggplot(pdata, aes(x=sex,y=group))+
  geom_point(aes(size=pct.exp,
                 color=avg.exp))+
  theme_classic()+
  scale_size_continuous(range = c(1, 9))+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+
  # scale_color_gradient(low="lightgrey",high="blue")+
  scale_color_viridis()+ ggtitle('Igfbp7')+ coord_flip()+
  labs(x=NULL,y=NULL)+theme(axis.text.x = element_text(colour = 'black', size = 16),
                            axis.text.y = element_text(colour = 'black', size = 16),
                            axis.title = element_text(colour = 'black', size = 16))
p
pdf('./figures/SC/main/Cd74_celltype.pdf', width = 5, height = 3)
print(p)
dev.off()



DotPlot(ASD.sce.qc, features = c('App'), group.by = 'group')







####MG_A####
MG_A <- as.list(as.data.frame(c('H2-Aa', 'H2-Ab1', 'H2-Eb1', 'H2-Ea', 'Cd74')))
names(MG_A) <- 'MG_A'
ASD.sce.qc <- AddModuleScore(ASD.sce.qc, features = MG_A, ctrl = 100, name = 'MG_A')

p<-VlnPlot(ASD.sce.qc, features = "MG_A1", pt.size = 0, group.by = 'group',sort = T,
           cols = pal_npg()(4))+
  stat_compare_means(comparisons = list(c('poly_IC', 'WT')))+
  ylim(-0.1, 1.2)+NoLegend()




####Additiional####
VlnPlot(subset(ASD.sce.qc, subset = celltype_main == 'Astrocytes'), features = "SASP1", pt.size = 0, group.by = 'group',
        cols = pal_npg()(2))+
  stat_compare_means(comparisons = list(c('poly', 'WT')))+
  ylim(-0.1, 0.2)+NoLegend()
 

geneset_RE = msigdbr(species = "Homo sapiens",#Homo sapiens
                     category = "C2",
                     subcategory = "CP:REACTOME"
) %>% dplyr::select(gs_name,gene_symbol)
geneset_RE <- geneset_RE %>% split(x = .$gene_symbol, f = .$gs_name)
RE_SASP <- geneset_RE[['REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP']]
RE_SASP <- as.list(as.data.frame(RE_SASP))

ASD_sn <- AddModuleScore(ASD_sn, features = RE_SASP, ctrl = 100, name = 'SASP')

#IGFBP7
#CDKN1A

#SASP1
#APP
# AST-FB           AST-PP   Microglia
VlnPlot(subset(ASD_sn, subset = cluster %in% c('AST-FB', 'AST-PP')), features = "SASP1", pt.size = 0, 
        group.by = 'diagnosis',
        cols = pal_npg()(2))+
  stat_compare_means(comparisons = list(c('ASD', 'Control')))+
  ylim(-0.2, 0.5)+NoLegend()



s.genes <- sapply(cc.genes$s.genes, str_to_title)
g2m.genes <- sapply(cc.genes$g2m.genes, str_to_title)
g2m.genes[30] <- 'Jpt1'
g2m.genes[16] <- 'Pimreg'
s.genes[15] <- 'Cenpu'
ASD.sce.qc <- CellCycleScoring(ASD.sce.qc, s.features = s.genes, 
                            g2m.features = g2m.genes, 
                            set.ident = TRUE)
p<-dittoBarPlot(subset(ASD.sce.qc, subset = celltype_main == 'Endothelial'), 'Phase', group.by = 'group',
             color.panel = pal_npg()(9))+
  theme(axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 10))+xlab('')
pdf('./figures/SC/main/DittoBarPlot_Endo.pdf', width = 4, height = 5)
print(p)
dev.off()

dittoBarPlot(subset(ASD.sce.qc, subset = celltype_main == 'Microgila'), 'Phase', group.by = 'group',
             color.panel = pal_npg()(9))+
  theme(axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 10))+xlab('')+ggtitle('Microgila')

dittoBarPlot(ASD.sce.qc, 'Phase', group.by = 'group',
             color.panel = pal_npg()(9))+
  theme(axis.text.x = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 10))+xlab('')+ggtitle('all')

# s.genes <- sapply(cc.genes$s.genes, capitalize_first_letter) 
# g2m.genes <- sapply(cc.genes$g2m.genes, capitalize_first_letter)
# g2m.genes[30] <- 'Jpt1'
# g2m.genes[16] <- 'Pimreg'
# s.genes[15] <- 'Cenpu'
# ASD.brain <- CellCycleScoring(ASD.brain, s.features = s.genes, 
#                                g2m.features = g2m.genes, 
#                                set.ident = TRUE)
# dittoBarPlot(ASD.brain, 'Phase', group.by = 'group')







s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ASD_sn <- CellCycleScoring(ASD_sn
                           , s.features = s.genes, 
                               g2m.features = g2m.genes, 
                               set.ident = TRUE)
dittoBarPlot(subset(ASD_sn, subset = cluster == 'Endothelial'), 'Phase', group.by = 'diagnosis')

dittoBarPlot(subset(ASD.sce.qc, subset = celltype_main == 'Endothelial'), 'Phase', group.by = 'group')


ASD.sce.qc <- AddModuleScore(ASD.sce.qc, features = list(geneset_RE$REACTOME_APOPTOSIS), ctrl = 100, name = 'APOPTOSIS')

ASD.sce.qc <- AddModuleScore(ASD.sce.qc, features = list(geneset_RE$REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS),
                            ctrl = 100, name = 'STRESS')

ASD.sce.qc <- AddModuleScore(ASD.sce.qc, features = list(geneset_RE$REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE),
                            ctrl = 100, name = 'OXIDATIVE_STRESS')

ASD.sce.qc <- AddModuleScore(ASD.sce.qc, features = list(geneset_RE$REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES),
                            ctrl = 100, name = 'FOXO')

ASD.sce.qc <- AddModuleScore(ASD.sce.qc, features = list(geneset_RE$REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS),
                            ctrl = 100, name = 'INTRINSIC')


ASD.sce.qc$group <- gsub('poly_IC', 'Positive', ASD.sce.qc$group)
ASD.sce.qc$group <- gsub('WT', 'Negative', ASD.sce.qc$group)

p<-VlnPlot(subset(ASD.sce.qc, celltype_main == 'Endothelial'), features = 'APOPTOSIS1',group.by = 'group', sort = T, pt.size = 0, 
        cols = pal_npg()(4))+
  ylim(-0.05, 0.2)+
  stat_compare_means(comparisons = list(c('Positive', 'Negative')),label = 'p.signif')+NoLegend()+xlab('')+ylab('')+
  theme(axis.text.x = element_text(size = 14,color = "black"),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.title.y = element_text(size = 16,color = "black"))+ggtitle('APOPTOSIS')
pdf('./figures/SC/VlnPlot_APO.pdf', width = 3.5, height = 5)
print(p)
dev.off()


p<-VlnPlot(subset(ASD.sce.qc, celltype_main == 'Endothelial'), features = 'STRESS1',group.by = 'group', sort = T, pt.size = 0, 
           cols = pal_npg()(4))+
  ylim(-0.05, 0.25)+
  stat_compare_means(comparisons = list(c('Positive', 'Negative')),label = 'p.signif')+NoLegend()+xlab('')+ylab('')+
  theme(axis.text.x = element_text(size = 14,color = "black"),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.title.y = element_text(size = 16,color = "black"))+ggtitle('STRESS')
pdf('./figures/SC/VlnPlot_STRESS.pdf', width = 3.5, height = 5)
print(p)
dev.off()


p<-VlnPlot(subset(ASD.sce.qc, celltype_main == 'Microgila'), features = 'APOPTOSIS1',group.by = 'group', sort = T, pt.size = 0, 
           cols = pal_npg()(4))+
  ylim(-0.05, 0.1)+
  stat_compare_means(comparisons = list(c('Positive', 'Negative')),label = 'p.signif')+NoLegend()+xlab('')+ylab('')+
  theme(axis.text.x = element_text(size = 14,color = "black"),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.title.y = element_text(size = 16,color = "black"))+ggtitle('APOPTOSIS')
pdf('./figures/SC/VlnPlot_APO_MG.pdf', width = 3.5, height = 5)
print(p)
dev.off()

p<-VlnPlot(subset(ASD.sce.qc, celltype_main == 'Microgila'), features = 'STRESS1',group.by = 'group', sort = T, pt.size = 0, 
           cols = pal_npg()(4))+
  ylim(-0.05, 0.25)+
  stat_compare_means(comparisons = list(c('Positive', 'Negative')),label = 'p.signif')+NoLegend()+xlab('')+ylab('')+
  theme(axis.text.x = element_text(size = 14,color = "black"),
        axis.text.y = element_text(size = 12,color = "black"),
        axis.title.y = element_text(size = 16,color = "black"))+ggtitle('STRESS')
pdf('./figures/SC/VlnPlot_STRESS_MG.pdf', width = 3.5, height = 5)
print(p)
dev.off()

