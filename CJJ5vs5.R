library(tidyverse)
library(DESeq2)
#import data
setwd("J:/大宝贝/")
mycounts<-read.table("Mergedata_noreplication.txt",header = TRUE,row.names = 1,sep ="\t" )
mycounts_new<-mycounts
a<-colnames(mycounts)
setwd("e:/5vs5/")
for (j in 1:6) {
  for (k in (j+1):7) {
    
    for (l in (k+1):8) {
      for (m in (l+1):9) {
        for (n in (m+1):10) {
          for (p in 11:16) {
            for (q in (p+1):17) {
              
              for (r in (q+1):18) {
                for (s in (r+1):19) {
                  
                  for (t in (s+1):20) {
                    mycounts<-mycounts_new[,c(a[j],a[k],a[l],a[m],a[n],a[p],a[q],a[r],a[s],a[t])]
                    condition<-factor(c(rep("WT",5),rep("KO",5)),levels = c("KO","WT"))
                    colData<-data.frame(row.names = colnames(mycounts),condition)
                    
                    dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
                    dds <- DESeq(dds)
                    
                    res= as.data.frame(results(dds))
                    res<-na.omit(res)
                    Irf4<-res[which(rownames(res) == "Irf4"),]
                    if(length(rownames(Irf4))!=0){
                      if(Irf4$pvalue < 0.05 & Irf4$log2FoldChange >0){
                        write.csv(res,file=paste("4vs4",j,"_",k,"_",l,"_",m,"_",n,"_",p,"_",q,"_",r,"_",s,"_",t,"_","diffgene.csv",sep = ""))
                        diff_gene=as.data.frame(res)
                        gene_list=diff_gene[,c("log2FoldChange","pvalue")]
                        colnames(gene_list)=c("logFC","padj")
                        gene_list$threshold = as.factor(abs(gene_list$logFC) > 1 & gene_list$padj < 0.05)
                        colored_point<-gene_list[gene_list$threshold == "TRUE",]
                        Spgenes<-gene_list[rownames(gene_list) == "Irf4" ,]
                        gene_list$threshold<-as.character(gene_list$threshold)
                        gene_list$threshold[which(rownames(gene_list) == "Irf4" )]<-"Irf4"
                        gene_list$threshold[which(gene_list$logFC >1 & gene_list$threshold ==TRUE)] <- "UP"
                        colnames(gene_list)[3]<-"Significant"
                        gene_list$Significant[which(gene_list$Significant == "TRUE")]<-"Down"
                        gene_list$Significant[which(gene_list$Significant == "FALSE")]<-"Not Sig"
                        Mycolors<-c("Gray","Red","Gray","Gray")
                        library("ggplot2")
                        pdf(paste("4vs4",j,"_",k,"_",l,"_",m,"_",n,"_",p,"_",q,"_",r,"_",s,"_",t,"_","vocano.pdf",sep = ""))
                        up<-floor(max(gene_list$logFC))+1
                        down<-floor(min(gene_list$logFC))-1
                        sig<-floor(max(-log10(gene_list$padj)))+1
                        g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=Significant)) + geom_point(alpha=0.8, size=1.75)  + xlim(c(down, up)) + ylim(c(0, sig)) +xlab("log2 fold change") + ylab("-log10 p-value") + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) + scale_color_manual(values = Mycolors)
                        print(g)
                        dev.off()
                      }
                    }
                    
                  }
                  
                  
                  
                }
              }
            }
          }
        }
        
        
      }
    }
  }
}
