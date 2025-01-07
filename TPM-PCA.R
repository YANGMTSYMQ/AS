library(tidyverse) 
library(data.table)
library(ggplot2)
library(factoextra)
library(ggplot2)
library(ggh4x)
library(limma)
library(data.table)
library(clusterProfiler)
library(AnnotationHub)
library(biomaRt)
library(xlsx)
library("readxl")
A1<-read.csv("../Counts/counts.csv",row.names = 1)
A2<- as.data.frame(fread("../Counts/COUNTS.txt"))
rownames(A1)<-A1$Geneid
rownames(A2)<-A1$Geneid
counts <- A1[,7:ncol(A1)]
geneid_efflen <- subset(A2,select = c("Geneid","Length"))
library(dplyr)
colnames(geneid_efflen) <- c("geneid","efflen")
geneid_efflen_fc <- geneid_efflen
dim(geneid_efflen)
efflen <- geneid_efflen[match(rownames(counts),
                              geneid_efflen$geneid),
                        "efflen"]
counts2TPM <- function(count=count, efflength=efflen){
  RPK <- count/(efflength/1000)     
  PMSC_rpk <- sum(RPK)/1e6       
  RPK/PMSC_rpk                    
}  
tpm <- as.data.frame(apply(counts,2,counts2TPM))
colSums(tpm)
write.csv(tpm, file = '../Counts/tpm.csv',row.names = T)
library(ggpubr)
library(gmodels)
library(limma)
getwd()
mytpm <- read.csv(file = "../Counts/tpm.csv",row.names = 1)
mytpm <- log10(mytpm+1)
mytpm <- na.omit(mytpm)
pca.info <- fast.prcomp(mytpm)
length(colnames(mytpm))
col <- c(rep("AD",8),rep("CON",10))
head(pca.info)
pca.data <- data.frame(sample = rownames(pca.info$rotation), 
                       Type= col, pca.info$rotation)
percentage<-round(pca.info$sdev / sum(pca.info$sdev) * 100,2)
pca.pecent<- pca.data[,-c(1,2)]
percentage<-paste(colnames(pca.pecent),"(", paste(as.character(percentage), "%", ")", sep=""))
ggplot(pca.data,aes(x=PC1,y=PC2,label=sample,color=Type,shape=Type ))+ 
  geom_point(size=5)+ geom_text(size=5)+
  theme_bw() +xlab(percentage[1]) +
  ylab(percentage[2])+ 
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"))
ggplot(pca.data, aes(x = PC1, y = PC2, label = sample, color = Type, shape = Type)) +
  geom_point(size = 5) + 
  geom_text(size = 5) +
  stat_ellipse(aes(group = Type), level = 0.95, linetype = 2, size = 1, color = "black") +  
  theme_bw() +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
library(ggh4x)
mytheme <- theme(panel.grid = element_blank(),
                 panel.background = element_blank(),
                 legend.key = element_blank(),
                 legend.position = "top",
                 axis.line = element_line(colour = "grey30"),
                 axis.ticks.length = unit(1.8, "mm"),
                 ggh4x.axis.ticks.length.minor = rel(0.6))
p <- ggplot(pca.data,aes(x=PC1,y=PC2,fill=Type, label = T))+
  stat_centroid(aes(xend = PC1, yend = PC2, colour = Type),
                geom = "segment", crop_other = F,
                alpha=0.3,size = 3,show.legend = F)+
  geom_point(size=3,alpha=0.7,
             color="white",shape = 21,show.legend = T)+
  scale_color_manual(name="",
                     values = c("#FF9999","#c77cff","cyan","yellow","red"))+
  scale_fill_manual(name="",
                    values = c("#FF9999","#c77cff","cyan","yellow","red"))+
  guides(x = "axis_truncated",y = "axis_truncated")+
  theme_bw()
library(ggplot2)
library(ggrepel) 

p <- ggplot(pca.data, aes(x = PC1, y = PC2, fill = Type, label = sample)) + 
  stat_centroid(aes(xend = PC1, yend = PC2, colour = Type),
                geom = "segment", crop_other = FALSE,
                alpha = 0.3, size = 1, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.7,
             color = "white", shape = 21, show.legend = TRUE) +
  scale_color_manual(name = "",
                     values = c("#FF9999", "#c77cff", "cyan", "yellow", "red")) +
  scale_fill_manual(name = "",
                    values = c("#FF9999", "#c77cff", "cyan", "yellow", "red")) +
  guides(x = "axis_truncated", y = "axis_truncated") +
  theme_bw() +
  theme(
    text = element_text(size = 12),  
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),  
    legend.title = element_blank(),  
    legend.text = element_text(size = 12),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(title = "PCA Analysis of Samples") +
  geom_label_repel(show.legend = FALSE, size = 4, box.padding = 0.5)  
mart<-useMart("ensembl")
dataset_list=as.data.frame(listDatasets(mart))

mart_oas=useMart("ensembl", "hsapiens_gene_ensembl")
list_gene=getBM(attributes = c("ensembl_gene_id","entrezgene_id","external_gene_name","description","chromosome_name"),
                filters = "ensembl_gene_id",
                values = rownames(A1),
                mart = mart_oas)
listgene<-list_gene[!duplicated(list_gene$ensembl_gene_id),]
A1$geneid<-rownames(A1)
genes <- merge(listgene,A1,by.x = "ensembl_gene_id",by.y = "geneid")
getwd()
genes<-subset(genes,!duplicated(genes$external_gene_name))
write.csv(genes,"../Counts/genecounts.csv")

counts <-read.csv("../Counts/genecounts.csv",row.names = 4)
counts <-counts[,-(1:11)]

