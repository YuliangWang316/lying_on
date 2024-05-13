library(patchwork)
library(Seurat)
library(dplyr)
library(future)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
options(warn=-1)
set.seed(1)
library(Hmisc)

GCB1.data<-Read10X("D:/GCB1_2.2.filtered_feature_bc_matrix/")
GCB2.data<-Read10X("D:/GCB2_2.2.filtered_feature_bc_matrix/")

GCB1.data<-as.data.frame(GCB1.data)
GCB2.data<-as.data.frame(GCB2.data)

# GCB1.metadata<-read.table("d:/GCB1metadata.txt",sep = "\t",header = TRUE,row.names = 1)
# GCB2.metadata<-read.table("d:/GCB2metadata.txt",sep = "\t",header = TRUE,row.names = 1)
# 
# GCB1.data<-GCB1.data[,rownames(GCB1.metadata)]
# GCB2.data<-GCB2.data[,rownames(GCB2.metadata)]

for (i in 1:length(colnames(GCB1.data))) {
  colnames(GCB1.data)[i] <- paste(colnames(GCB1.data)[i],"GCB1",i,sep = "-")  
}

for (i in 1:length(colnames(GCB2.data))) {
  colnames(GCB2.data)[i] <- paste(colnames(GCB2.data)[i],"GCB2",i,sep = "-")  
}

# rownames(GCB1.metadata)<-colnames(GCB1.data)
# rownames(GCB2.metadata)<-colnames(GCB2.data)
# 
# GCB1.metadata$group<-rep("GCB1",length(colnames(GCB1.data)))
# GCB2.metadata$group<-rep("GCB2",length(colnames(GCB2.data)))

# pbmc.metadata<-rbind(GCB1.metadata,GCB2.metadata)

pbmc.data<-cbind(GCB1.data,GCB2.data)
pbmc.data <- as.data.frame(pbmc.data)
pbmc.data_t<-as.data.frame(t(pbmc.data))
pbmc.data_t<-filter(pbmc.data_t,Ighd <=0)
pbmc.data<-as.data.frame(t(pbmc.data_t))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB",min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 30)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F,npcs = 100)
ElbowPlot(pbmc,ndims = 100)

library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ref <- readRDS("D:/CJJ immunity/E-MTAB-9005/SEURAT_OBJECTS/HumanTonsil_BCells_scRNA_SeuratObject.rds")
ref<-UpdateSeuratObject(ref)
ref.metadata<-ref@meta.data[,c("Sample","CellType")]
ref.metadata$mouse_ID<-rownames(ref.metadata)
colnames(ref.metadata)[2]<-"Clustername"
library(Hmisc)
ref.data<-ref@assays$RNA@counts
rownames(ref.data)<-capitalize(tolower(rownames(ref.data)))
ds_seurat<-SummarizedExperiment(assays=list(counts=ref.data),colData = ref.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc)<-pred.hesc@listData$labels
pbmc@meta.data$Cluster<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$Cluster)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:11)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:11)
DimPlot(pbmc,reduction = "ldaumap")
DimPlot(pbmc,reduction = "ldatsne")
VlnPlot(pbmc,features = c("Irf4","Gls","Kdm6b","Ly75"),sort = TRUE)
DotPlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"))
new.cluster.ids <- c("DZ_2", "Fraction3_2", "Fraction3_1", "Fraction3_3", "DZ_1", "DZ_3",
                     "Fraction2_3", "Fraction2_2", "prePB","Fraction2_1","prePB","prePB","prePB","prePB","prePB")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, label = TRUE, pt.size = 0.5) + NoLegend()
pbmc@meta.data$newcluster<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$newcluster)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:9)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:9)
DimPlot(pbmc,reduction = "ldaumap")
DimPlot(pbmc,reduction = "ldatsne")
VlnPlot(pbmc,features = c("Irf4","Gls","Kdm6b","Ly75"),sort = TRUE)
DotPlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"))

Idents(pbmc)<-pbmc$Cluster
new.cluster.ids <- c("DZ_2", "Fraction3_2", "Fraction3_1", "Fraction3_3", "DZ_1", "DZ_3",
                     "Fraction2_3", "Fraction2_2", "prePB","Fraction2_1","Fraction1_2","prePB","prePB","prePB","prePB")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, label = TRUE, pt.size = 0.5) + NoLegend()
pbmc@meta.data$newnewcluster<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$newnewcluster)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:9)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:9)
DimPlot(pbmc,reduction = "ldaumap")
DimPlot(pbmc,reduction = "ldatsne")
VlnPlot(pbmc,features = c("Irf4","Gls","Kdm6b","Ly75"),sort = TRUE)
DotPlot(pbmc,features =c("Irf4","Gls","Kdm6b","Ly75"))

# library(CelliD)
# pathway<-read.table("D:/GSE109732_GeneExpression_DZ_Fraction1-3_GC-PB.txt",sep = "\t",header = TRUE)
# pathway<-pathway[!duplicated(pathway$Gene),]
# rownames(pathway)<-pathway$Gene
# pathway<-pathway[,-1]
# pathway<-as.list(pathway)
# pbmc<-RunMCA(pbmc)
# HGT_organoid_gs <- RunCellHGT(pbmc, pathways = pathway, n.features = 200,dims = 1:50)
# HGT_organoid_gs_prediction <- rownames(HGT_organoid_gs)[apply(HGT_organoid_gs, 2, which.max)]
# HGT_organoid_signif <- ifelse(apply(HGT_organoid_gs, 2, max)>2, yes = HGT_organoid_gs_prediction, "unassigned")
# pbmc$HGT_organoid_prediction <- HGT_organoid_signif
