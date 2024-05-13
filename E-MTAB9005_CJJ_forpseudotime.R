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

library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ref <- readRDS("D:/CJJ immunity/E-MTAB-9005/SEURAT_OBJECTS/HumanTonsil_BCells_scRNA_SeuratObject.rds")
ref<-UpdateSeuratObject(ref)
pbmc2<-readRDS("e:/CJJGSE109732.rds")
pbmc2.data<-pbmc2@assays$RNA@counts
pbmc2.metadata<-pbmc2@meta.data[,c(5,22,23,25)]
remove(pbmc2)
gc()
pbmc<-ref
# pbmc.metadata<-ref.metadata
remove(ref)
gc()
Idents(pbmc)<-pbmc$CellType
pbmc_new<-subset(pbmc,idents = c("Cycling","GC","FCRL2/3high GC","DZ GC","LZ GC","prePB"))
pbmc2.metadata$mouse_ID<-rownames(pbmc2.metadata)
pbmc2.metadata$Clustername<-pbmc2.metadata$newcluster
rownames(pbmc2.data)<-toupper(rownames(pbmc2.data))
ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc2.data),colData = pbmc2.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
pbmc<-subset(pbmc,idents = c("Cycling","GC","FCRL2/3high GC","DZ GC","LZ GC","prePB","Plasmablast"))
pbmc@meta.data$newcluster<-as.character(pbmc@meta.data$CellType)
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(rownames(pbmc@meta.data)[j] ==  rownames(pred.hesc)[i]){
      pbmc@meta.data$newcluster[j]<-pred.hesc@listData$labels[i]
    }
  }
}

pbmc_new_new<-subset(pbmc,idents = c("GC","FCRL2/3high GC","DZ GC","LZ GC","prePB"))
pred.hesc <- SingleR(test = pbmc_new_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
pbmc@meta.data$newnewcluster<-as.character(pbmc@meta.data$CellType)
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(rownames(pbmc@meta.data)[j] ==  rownames(pred.hesc)[i]){
      pbmc@meta.data$newnewcluster[j]<-pred.hesc@listData$labels[i]
    }
  }
}
remove(i,j,pred.hesc)
gc()
