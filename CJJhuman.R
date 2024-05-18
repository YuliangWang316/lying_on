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
library(cowplot)
library(tidyverse)
library(SummarizedExperiment)
library(scuttle)
library(SingleR)
PC.data <- Read10X(data.dir = "E:/P23042711/PC/20230630/Matrix/")
DP.data <- Read10X(data.dir = "E:/P23042711/DP/20230630/Matrix/")
DN.data <- Read10X(data.dir = "E:/P23042711/DN/20230630/Matrix/")
NB.data <- Read10X(data.dir = "E:/P23042711/NB/20230630/Matrix/")

PC.data <- as.data.frame(PC.data)
DP.data <- as.data.frame(DP.data)
DN.data <- as.data.frame(DN.data)
NB.data <- as.data.frame(NB.data)

for (i in 1:length(colnames(PC.data))) {
  colnames(PC.data)[i] <- paste(colnames(PC.data)[i],"PC",i,sep = "-")  
}
for (i in 1:length(colnames(DP.data))) {
  colnames(DP.data)[i] <- paste(colnames(DP.data)[i],"DP",i,sep = "-")  
}
for (i in 1:length(colnames(DN.data))) {
  colnames(DN.data)[i] <- paste(colnames(DN.data)[i],"DN",i,sep = "-")  
}
for (i in 1:length(colnames(NB.data))) {
  colnames(NB.data)[i] <- paste(colnames(NB.data)[i],"NB",i,sep = "-")  
}


PC.metadata<-data.frame(colnames(PC.data),rep("PC",length(colnames(PC.data))))
colnames(PC.metadata)<-c("barcode","group")
DP.metadata<-data.frame(colnames(DP.data),rep("DP",length(colnames(DP.data))))
colnames(DP.metadata)<-c("barcode","group")
DN.metadata<-data.frame(colnames(DN.data),rep("DN",length(colnames(DN.data))))
colnames(DN.metadata)<-c("barcode","group")
NB.metadata<-data.frame(colnames(NB.data),rep("NB",length(colnames(NB.data))))
colnames(NB.metadata)<-c("barcode","group")

pbmc.metadata<-rbind(PC.metadata,DP.metadata,DN.metadata,NB.metadata)
rownames(pbmc.metadata)<-pbmc.metadata[,1]
pbmc.data<-cbind(PC.data,DP.data,DN.data,NB.data)

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata,min.cells = 3, min.features = 200)
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500& percent.mt < 10)
remove(DN.data,DP.data,NB.data,PC.data)
remove(DN.metadata,DP.metadata,NB.metadata,PC.metadata)
remove(pbmc.data,pbmc.metadata,i)
Idents(pbmc)<-pbmc$group
pbmc<-NormalizeData(pbmc)
pbmc<- ScaleData(pbmc, features = rownames(pbmc))
NB<-subset(pbmc,idents = c("NB"))
DN<-subset(pbmc,idents = c("DN"))
DP<-subset(pbmc,idents = c("DP"))
PC<-subset(pbmc,idents = c("PC"))

NB_new<-NB[,sample(1:ncol(NB),500)]
DN_new<-DN[,sample(1:ncol(DN),500)]
DP_new<-DP[,sample(1:ncol(DP),500)]
PC_new<-PC[,sample(1:ncol(PC),500)]
pbmc_new<-pbmc.big <- merge(NB_new, y = c(DN_new, DP_new,PC_new), add.cell.ids = c("NB", "DN", "DP","PC"), project = "Total")
pbmc_new <- NormalizeData(pbmc_new)
pbmc_new <- FindVariableFeatures(pbmc_new, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc_new)
pbmc_new <- ScaleData(pbmc_new, features = all.genes)
pbmc_new <- RunPCA(pbmc_new, features = VariableFeatures(object = pbmc_new))
ElbowPlot(pbmc_new)
pbmc_new <- FindNeighbors(pbmc_new, dims = 1:20)
pbmc_new <- FindClusters(pbmc_new, resolution = 1.2)
pbmc_new <- RunUMAP(pbmc_new, dims = 1:20)
pbmc_new <- RunTSNE(pbmc_new, dims = 1:20)
DimPlot(pbmc_new, reduction = "umap")
DimPlot(pbmc_new, reduction = "umap",split.by = "group")
DimPlot(pbmc_new, reduction = "tsne")
DimPlot(pbmc_new, reduction = "tsne",split.by = "group")
DotPlot(pbmc_new,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"))
VlnPlot(pbmc_new,features = c("GLS","GOT1","GOT2","GLUD1","GPT2","MYC","IRF4","LY75","KDM6B"),pt.size = 0,sort = TRUE)
DimPlot(pbmc_new,reduction = "pca")
DimPlot(pbmc_new,reduction = "pca",split.by = "group")
pbmc_new <- RunLDA(pbmc_new, labels = pbmc_new$group)
pbmc_new <- RunUMAP(pbmc_new ,reduction = "lda",dims = 1:3,reduction.name = "lda_umap")
pbmc_new <- RunTSNE(pbmc_new, reduction = "lda",dims = 1:3,reduction.name = "lda_tsne" )
Idents(pbmc_new)<-pbmc_new$group
DimPlot(pbmc_new,reduction = "lda")
DimPlot(pbmc_new,reduction = "lda_umap")
DimPlot(pbmc_new,reduction = "lda_tsne")
pbmc_new_new<-subset(pbmc_new,idents =c("DN","DP"))

Glutamine_1<-read.table("d:/GSE60927.txt",sep = "\t")
Glutamine_1_new<-Glutamine_1[1:length(rownames(Glutamine_1)),]
Glutamine_1_new<-toupper(Glutamine_1_new)
Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmc_new_new))
Glutamine1<-FetchData(pbmc_new_new,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc_new_new,vars = rownames(pbmc_new_new))
Data<-data.frame()
remove(i)
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine1<-mean(Glutamine1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
remove(j)
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_Glutamine1)
Sortdatalo<-filter(Sortdata,S < average_Glutamine1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((Glutamine_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc_new_new,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc_new_new@meta.data$plasmacell<-Glutamine1$Average-control$Average
remove(p,control,Data,Data_T,Glutamine_1,Glutamine1,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total,average_Glutamine1,Glutamine_1_new,higene,logene,k)
DotPlot(pbmc_new_new,features = "plasmacell")
remove(all.genes,DN,DN_new,DP,DP_new,NB,NB_new,PC,PC_new)
remove(pbmc_new,pbmc.big)
gc()
saveRDS(pbmc_new_new,"HumanDNDP_sampling_enrichusingGSE60927.rds")
remove(pbmc_new_new)
gc()
pbmc_new_new<-subset(pbmc,idents = c("DN","DP"))

Glutamine_1<-read.table("d:/GSE60927.txt",sep = "\t")
Glutamine_1_new<-Glutamine_1[1:length(rownames(Glutamine_1)),]
Glutamine_1_new<-toupper(Glutamine_1_new)
Glutamine_1_new<-intersect(Glutamine_1_new,rownames(pbmc_new_new))
Glutamine1<-FetchData(pbmc_new_new,vars = Glutamine_1_new)
for (i in 1:length(rownames(Glutamine1))) {
  Glutamine1$Average[i]<-mean(as.numeric(Glutamine1[i,]))
}

Total<-FetchData(pbmc_new_new,vars = rownames(pbmc_new_new))
Data<-data.frame()
remove(i)
for (j in 1:length(colnames(Total))) {
  Data[1,j]<-mean(Total[,j])
}
average_Glutamine1<-mean(Glutamine1$Average)
colnames(Data)<-colnames(Total)
Data_T<-as.data.frame(t(Data))
remove(j)
rownamesDataT<-as.data.frame(rownames(Data_T))
Data_T<-cbind(Data_T,rownamesDataT)
colnames(Data_T)<-c("S","Gene")
Sortdata<-Data_T[order(Data_T$S,decreasing = TRUE),]

Sortdatahi<-filter(Sortdata,S > average_Glutamine1)
Sortdatalo<-filter(Sortdata,S < average_Glutamine1)
higene<-rownames(Sortdatahi[(length(rownames(Sortdatahi))-(10*length((Glutamine_1_new)))/2+1):length(rownames(Sortdatahi)),])
logene<-rownames(Sortdatalo[1:(10*length((Glutamine_1_new))/2),])
k<-c(higene,logene)
control<-FetchData(pbmc_new_new,vars = k )
for (p in 1:length(rownames(control))) {
  control$Average[p]<-mean(as.numeric(control[p,]))
}

pbmc_new_new@meta.data$plasmacell<-Glutamine1$Average-control$Average
remove(p,control,Data,Data_T,Glutamine_1,Glutamine1,rownamesDataT,Sortdata,Sortdatahi,Sortdatalo,Total,average_Glutamine1,Glutamine_1_new,higene,logene,k)
DotPlot(pbmc_new_new,features = "plasmacell")

saveRDS(pbmc_new_new,"d:/HumanDNDP_enrichusingGSE60927.rds")
remove(pbmc_new_new)
gc()
remove(pbmc)
gc()
pbmc1<-readRDS("h:/20240514_scHumanDNDP_GSE60927enrichment/HumanDNDP_enrichusingGSE60927.rds")
pbmc2<-readRDS("h:/20240514_scHumanDNDP_GSE60927enrichment/HumanDNDP_sampling_enrichusingGSE60927.rds")


pbmc1_DP<-subset(pbmc1,idents =c("DP"))
pbmc1_DN<-subset(pbmc1,idents =c("DN"))
pbmc1.data<-cbind(pbmc1_DP@assays$RNA@counts,pbmc1_DN@assays$RNA@counts)
pbmc1.data<-as.data.frame(pbmc1.data)
pbmc1.data<-cbind(rownames(pbmc1.data),pbmc1.data)
colnames(pbmc1.data)[1]<-"Gene_Symbol"
write.table(pbmc1.data,"HumanDNDP_enrichusingGSE60927.txt",sep = "\t",row.names = FALSE,quote = FALSE)
pbmc2_DP<-subset(pbmc2,idents =c("DP"))
pbmc2_DN<-subset(pbmc2,idents =c("DN"))
pbmc2.data<-cbind(pbmc2_DP@assays$RNA@counts,pbmc2_DN@assays$RNA@counts)
pbmc2.data<-as.data.frame(pbmc2.data)
pbmc2.data<-cbind(rownames(pbmc2.data),pbmc2.data)
colnames(pbmc2.data)[1]<-"Gene_Symbol"
write.table(pbmc2.data,"HumanDNDP_sampling_enrichusingGSE60927.txt",sep = "\t",row.names = FALSE,quote = FALSE)

pbmc1.data<-cbind(pbmc1.data[,c(1)],rep("NA",length(rownames(pbmc1.data))),pbmc1.data[,c(2:18835)])
colnames(pbmc1.data)[1:2]<-c("Gene_Symbol","Gene_Description")
line1<-c("#1.2",rep("",18835))
line2<-c("24250","18834",rep("",18834))
line3<-colnames(pbmc1.data)
# pbmc1.data$Gene_Symbol<-factor(pbmc1.data$Gene_Symbol,levels = pbmc1.data$Gene_Symbol)
pbmc1.data<-rbind(line1,line2,line3,pbmc1.data)
write.table(pbmc1.data,"HumanDNDP_enrichusingGSE60927_forGSEA.txt",sep = "\t",row.names = FALSE,quote = FALSE,col.names = FALSE)

pbmc2.data<-cbind(pbmc2.data[,c(1)],rep("NA",length(rownames(pbmc2.data))),pbmc2.data[,c(2:1001)])
colnames(pbmc2.data)[1:2]<-c("Gene_Symbol","Gene_Description")
line1<-c("#1.2",rep("",1001))
line2<-c("24250","1000",rep("",1000))
line3<-colnames(pbmc2.data)
# pbmc2.data$Gene_Symbol<-factor(pbmc2.data$Gene_Symbol,levels = pbmc2.data$Gene_Symbol)
pbmc2.data<-rbind(line1,line2,line3,pbmc2.data)
write.table(pbmc2.data,"HumanDNDP_sampling_enrichusingGSE60927_forGSEA.txt",sep = "\t",row.names = FALSE,quote = FALSE,col.names = FALSE)
line1<-c("18834","2","1",rep("",18831))
line2<-c("#","DP","DN",rep("",18831))
line3<-c(rep("DP",1668),rep("DN",17166))
pbmc1.metadata<-rbind(line1,line2,line3)
write.table(pbmc1.metadata,file = "HumanDNDP_enrichusingGSE60927_forGSEA_sample.txt",sep = "\t",row.names = FALSE,quote = FALSE,col.names = FALSE)
line1<-c("1000","2","1",rep("",997))
line2<-c("#","DP","DN",rep("",997))
line3<-c(rep("DP",500),rep("DN",500))
pbmc2.metadata<-rbind(line1,line2,line3)
write.table(pbmc2.metadata,file = "HumanDNDP_sample_enrichusingGSE60927_forGSEA_sample.txt",sep = "\t",row.names = FALSE,quote = FALSE,col.names = FALSE)


library(monocle)
trace('project2MST',edit = T,where = asNamespace("monocle"))
data<-as.sparse(pbmc_new_new@assays$RNA@counts)
pd <-pbmc_new_new@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,
                              expressionFamily = VGAM::negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 
# 
# 
pbmc.marker<-FindAllMarkers(pbmc_new_new,only.pos = TRUE,min.pct = 0)
diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~group",cores = 20)

pbmcmarkers_new<-pbmc.marker[which(pbmc.marker$p_val_adj < 0.05 ),]
pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
# ordering_genes <- diff_test_res$gene_short_name
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-90))
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2)

monocle_cds <-orderCells(monocle_cds )
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "group",cell_size = 0.75,)
plot_cell_trajectory(monocle_cds, color_by = "group",cell_size = 0.75,)+ facet_wrap(~group, nrow = 1)
plotdf=pData(monocle_cds)
library(ggridges)
mycolor<-c("#619CFF","#00BA38","#F8766D")

ggplot(plotdf, aes(x=Pseudotime,y=group,fill=group))+
  geom_density_ridges(scale=1) +
  
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )+scale_fill_manual(values = mycolor)
