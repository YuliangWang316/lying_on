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

GCB1.data<-read.table("D:/CJJ immunity/GSE139891_RAW/GSM4148370_GC1_umi_grch38.txt/RK10001_counts.txt",sep = "\t",header = TRUE)
GCB2.data<-read.table("D:/CJJ immunity/GSE139891_RAW/GSM4148371_GC2_umi_grch38.txt/RK10004_counts.txt",sep = "\t",header = TRUE)
GCB3.data<-read.table("D:/CJJ immunity/GSE139891_RAW/GSM4148372_DZ1_umi_grch38.txt/RK10002_counts.txt",sep = "\t",header = TRUE)
GCB4.data<-read.table("D:/CJJ immunity/GSE139891_RAW/GSM4148373_DZ2_umi_grch38.txt/RK10005_counts.txt",sep = "\t",header = TRUE)
GCB5.data<-read.table("D:/CJJ immunity/GSE139891_RAW/GSM4148374_LZ1_umi_grch38.txt/RK10003_counts.txt",sep = "\t",header = TRUE)
GCB6.data<-read.table("D:/CJJ immunity/GSE139891_RAW/GSM4148375_LZ2_umi_grch38.txt/RK10006_counts.txt",sep = "\t",header = TRUE)
GCB7.data<-read.table("D:/CJJ immunity/GSE139891_RAW/GSM4560814_GC3a_umi.txt/counts_no_bad.txt",sep = "\t",header = TRUE)
GCB8.data<-read.table("D:/CJJ immunity/GSE139891_RAW/GSM4560815_GC3b_umi.txt/counts_no_bad.txt",sep = "\t",header = TRUE)

GCB1.data<-GCB1.data[!duplicated(GCB1.data$X),]
rownames(GCB1.data)<-GCB1.data[,1]
GCB1.data<-GCB1.data[,-1]

GCB2.data<-GCB2.data[!duplicated(GCB2.data$X),]
rownames(GCB2.data)<-GCB2.data[,1]
GCB2.data<-GCB2.data[,-1]

GCB3.data<-GCB3.data[!duplicated(GCB3.data$X),]
rownames(GCB3.data)<-GCB3.data[,1]
GCB3.data<-GCB3.data[,-1]

GCB4.data<-GCB4.data[!duplicated(GCB4.data$X),]
rownames(GCB4.data)<-GCB4.data[,1]
GCB4.data<-GCB4.data[,-1]

GCB5.data<-GCB5.data[!duplicated(GCB5.data$X),]
rownames(GCB5.data)<-GCB5.data[,1]
GCB5.data<-GCB5.data[,-1]

GCB6.data<-GCB6.data[!duplicated(GCB6.data$X),]
rownames(GCB6.data)<-GCB6.data[,1]
GCB6.data<-GCB6.data[,-1]

GCB7.data<-GCB7.data[!duplicated(GCB7.data$X),]
rownames(GCB7.data)<-GCB7.data[,1]
GCB7.data<-GCB7.data[,-1]

GCB8.data<-GCB8.data[!duplicated(GCB8.data$X),]
rownames(GCB8.data)<-GCB8.data[,1]
GCB8.data<-GCB8.data[,-1]

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

for (i in 1:length(colnames(GCB3.data))) {
  colnames(GCB3.data)[i] <- paste(colnames(GCB3.data)[i],"GCB3",i,sep = "-")  
}

for (i in 1:length(colnames(GCB4.data))) {
  colnames(GCB4.data)[i] <- paste(colnames(GCB4.data)[i],"GCB4",i,sep = "-")  
}

for (i in 1:length(colnames(GCB5.data))) {
  colnames(GCB5.data)[i] <- paste(colnames(GCB5.data)[i],"GCB5",i,sep = "-")  
}

for (i in 1:length(colnames(GCB6.data))) {
  colnames(GCB6.data)[i] <- paste(colnames(GCB6.data)[i],"GCB6",i,sep = "-")  
}

for (i in 1:length(colnames(GCB7.data))) {
  colnames(GCB7.data)[i] <- paste(colnames(GCB7.data)[i],"GCB7",i,sep = "-")  
}

for (i in 1:length(colnames(GCB8.data))) {
  colnames(GCB8.data)[i] <- paste(colnames(GCB8.data)[i],"GCB8",i,sep = "-")  
}

GCB1.metadata<-data.frame(colnames(GCB1.data))
GCB2.metadata<-data.frame(colnames(GCB2.data))
GCB3.metadata<-data.frame(colnames(GCB3.data))
GCB4.metadata<-data.frame(colnames(GCB4.data))
GCB5.metadata<-data.frame(colnames(GCB5.data))
GCB6.metadata<-data.frame(colnames(GCB6.data))
GCB7.metadata<-data.frame(colnames(GCB7.data))
GCB8.metadata<-data.frame(colnames(GCB8.data))

rownames(GCB1.metadata)<-colnames(GCB1.data)
rownames(GCB2.metadata)<-colnames(GCB2.data)
rownames(GCB3.metadata)<-colnames(GCB3.data)
rownames(GCB4.metadata)<-colnames(GCB4.data)
rownames(GCB5.metadata)<-colnames(GCB5.data)
rownames(GCB6.metadata)<-colnames(GCB6.data)
rownames(GCB7.metadata)<-colnames(GCB7.data)
rownames(GCB8.metadata)<-colnames(GCB8.data)

GCB1.metadata$group<-rep("GCB1",length(colnames(GCB1.data)))
GCB2.metadata$group<-rep("GCB2",length(colnames(GCB2.data)))
GCB3.metadata$group<-rep("DZ1",length(colnames(GCB3.data)))
GCB4.metadata$group<-rep("DZ2",length(colnames(GCB4.data)))
GCB5.metadata$group<-rep("LZ1",length(colnames(GCB5.data)))
GCB6.metadata$group<-rep("LZ2",length(colnames(GCB6.data)))
GCB7.metadata$group<-rep("GCB3a",length(colnames(GCB7.data)))
GCB8.metadata$group<-rep("GCB3b",length(colnames(GCB8.data)))

GCB1.metadata$group2<-rep("GCB",length(colnames(GCB1.data)))
GCB2.metadata$group2<-rep("GCB",length(colnames(GCB2.data)))
GCB3.metadata$group2<-rep("DZ",length(colnames(GCB3.data)))
GCB4.metadata$group2<-rep("DZ",length(colnames(GCB4.data)))
GCB5.metadata$group2<-rep("LZ",length(colnames(GCB5.data)))
GCB6.metadata$group2<-rep("LZ",length(colnames(GCB6.data)))
GCB7.metadata$group2<-rep("GCB3",length(colnames(GCB7.data)))
GCB8.metadata$group2<-rep("GCB3",length(colnames(GCB8.data)))

colnames(GCB1.metadata)<-"ID"
colnames(GCB2.metadata)<-"ID"
colnames(GCB3.metadata)<-"ID"
colnames(GCB4.metadata)<-"ID"
colnames(GCB5.metadata)<-"ID"
colnames(GCB6.metadata)<-"ID"
colnames(GCB7.metadata)<-"ID"
colnames(GCB8.metadata)<-"ID"
pbmc.metadata<-rbind(GCB1.metadata,GCB2.metadata,GCB3.metadata,GCB4.metadata,GCB5.metadata,GCB6.metadata,GCB7.metadata,GCB8.metadata)
remove(GCB1.metadata,GCB2.metadata,GCB3.metadata,GCB4.metadata,GCB5.metadata,GCB6.metadata,GCB7.metadata,GCB8.metadata)

remove(i)

a<-intersect(rownames(GCB1.data),rownames(GCB2.data))
b<-intersect(a,rownames(GCB3.data))
c<-intersect(b,rownames(GCB4.data))
d<-intersect(c,rownames(GCB5.data))
e<-intersect(d,rownames(GCB6.data))
f<-intersect(e,rownames(GCB7.data))
g<-intersect(f,rownames(GCB8.data))

GCB1.data<-GCB1.data[g,]
GCB2.data<-GCB2.data[g,]
GCB3.data<-GCB3.data[g,]
GCB4.data<-GCB4.data[g,]
GCB5.data<-GCB5.data[g,]
GCB6.data<-GCB6.data[g,]
GCB7.data<-GCB7.data[g,]
GCB8.data<-GCB8.data[g,]

pbmc.data<-cbind(GCB1.data,GCB2.data,GCB3.data,GCB4.data,GCB5.data,GCB6.data,GCB7.data,GCB8.data)
remove(a,b,c,d,e,f,g)
remove(GCB1.data,GCB2.data,GCB3.data,GCB4.data,GCB5.data,GCB6.data,GCB7.data,GCB8.data)
gc()
colnames(pbmc.metadata)[2]<-"group"
colnames(pbmc.metadata)[3]<-"group2"
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "GCB",min.cells = 3, min.features = 200,meta.data = pbmc.metadata)
remove(pbmc.data,pbmc.metadata)
gc()
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 30)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose=F)
VariableFeaturePlot(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes, verbose=F)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose=F,npcs = 100)
ElbowPlot(pbmc,ndims = 100)

gc()
library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
ref.data<-read.table("d:/CJJ immunity/GSE139891_RAW/GSE139891.txt",sep = "\t",header = TRUE,check.names = FALSE)
ref.data<-ref.data[!duplicated(ref.data$Gene),]
ref.data<-na.omit(ref.data)
rownames(ref.data)<-ref.data$Gene
ref.data<-ref.data[,-1]
ref.metadata<-read.table("d:/CJJ immunity/GSE139891_RAW/GSE139891metadata.txt",sep = "\t",header = TRUE)
rownames(ref.metadata)<-ref.metadata$ID

ref.data<-ref.data[,rownames(ref.metadata)]
# library(Hmisc)
# rownames(ref.data)<-capitalize(tolower(rownames(ref.data)))
ref <- CreateSeuratObject(counts = ref.data, project = "ref3k",meta.data = ref.metadata,min.cells = 3)
colnames(ref.metadata)<-c("mouse_ID","Clustername")
ds_seurat<-SummarizedExperiment(assays=list(counts=ref@assays$RNA@counts),colData = ref.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)

Idents(pbmc)<-pred.hesc@listData$labels
pbmc@meta.data$Cluster<-Idents(pbmc)
pbmc<-RunLDA(pbmc,labels = pbmc$Cluster)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "ldaumap",dims = 1:12)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "ldatsne",dims = 1:12)
DimPlot(pbmc,reduction = "ldaumap")
DimPlot(pbmc,reduction = "ldatsne")
VlnPlot(pbmc,features = c("IRF4","GLS","KDM6B","LY75"),sort = TRUE)
DotPlot(pbmc,features =c("IRF4","GLS","KDM6B","LY75"))
remove(all.genes,ds_seurat,pred.hesc)
remove(ref,ref.data,ref.metadata)
gc()
pbmc2<-readRDS("e:/CJJGSE109732.rds")
pbmc2.data<-pbmc2@assays$RNA@counts
pbmc2.metadata<-pbmc2@meta.data[,c(5,22,23,25)]
remove(pbmc2)
gc()
rownames(pbmc2.data)<-toupper(rownames(pbmc2.data))

Idents(pbmc)<-pbmc$Cluster
pbmc_new<-subset(pbmc,idents = c("INT e (Intermediate) avg log2(tpm + 1)","DZ b (Dark Zone) avg log2(tpm + 1)","INT c (Intermediate) avg log2(tpm + 1)",
                                 "DZ a (Dark Zone) avg log2(tpm + 1)","INT b (Intermediate) avg log2(tpm + 1)","LZ b (Light Zone) avg log2(tpm + 1)",
                                 "INT a (Intermediate) avg log2(tpm + 1)","LZ a (Light Zone) avg log2(tpm + 1)","DZ c (Dark Zone) avg log2(tpm + 1)",
                                 "INT d (Intermediate) avg log2(tpm + 1)","PreM (Memory precursors) avg log2(tpm + 1)"))
pbmc2.metadata$mouse_ID<-rownames(pbmc2.metadata)
pbmc2.metadata$Clustername<-pbmc2.metadata$newcluster

ds_seurat<-SummarizedExperiment(assays=list(counts=pbmc2.data),colData = pbmc2.metadata[,c("mouse_ID","Clustername")])
ds_seurat<-logNormCounts(ds_seurat)
library(SingleR)
pred.hesc <- SingleR(test = pbmc_new@assays$RNA@counts, ref = ds_seurat, assay.type.test=1,
                     labels = ds_seurat$Clustername)
pbmc@meta.data$newcluster<-rep("PC",length(rownames(pbmc@meta.data)))
for (i in 1:length(rownames(pred.hesc))) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(rownames(pbmc@meta.data)[j] ==  rownames(pred.hesc)[i]){
      pbmc@meta.data$newcluster[j]<-pred.hesc@listData$labels[i]
    }
  }
}
remove(i,j,pbmc_new,pred.hesc,ds_seurat)
gc()
