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
ref.data<-read.table("d:/CJJ immunity/GSE188617_counts_info.tsv/counts_info.tsv",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)

ref.metadata<-read.table("d:/CJJ immunity/GSE188617_counts_info.tsv/vdj_info_new.txt",sep = "\t",header = TRUE,row.names = 1)
ref.metadata$ID<-rownames(ref.metadata)

ref.data<-ref.data[,rownames(ref.metadata)]
# library(Hmisc)
# rownames(ref.data)<-capitalize(tolower(rownames(ref.data)))
ref <- CreateSeuratObject(counts = ref.data, project = "ref3k",meta.data = ref.metadata,min.cells = 3)
pbmc2<-readRDS("e:/CJJGSE109732.rds")
pbmc2.data<-pbmc2@assays$RNA@counts
pbmc2.metadata<-pbmc2@meta.data[,c(5,22,23,25)]
remove(pbmc2)
gc()
Idents(ref)<-ref$Cluster
ref_new<-subset(ref,idents = c("DZ-4","DZ-3","INT-4","LZ-2","INT-1","DZ-1","LZ-1","INT-2","DZ-2","PreM","Undetermined","INT-3","FCRL2/3","LZ-3"))
pbmc2.metadata$mouse_ID<-rownames(pbmc2.metadata)
pbmc2.metadata$Clustername<-pbmc2.metadata$newcluster
rownames(pbmc2.data)<-toupper(rownames(pbmc2.data)) ##big upper   
ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc2.data),colData = pbmc2.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = ref_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
ref@meta.data$newcluster<-rep("PC",length(rownames(ref@meta.data)))
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(ref@meta.data))) {
    if(rownames(ref@meta.data)[j] ==  rownames(pred.hesc)[i]){
      ref@meta.data$newcluster[j]<-pred.hesc@listData$labels[i]
    }
  }
}
remove(i,j,ref_new,pred.hesc,ds_seurat)
