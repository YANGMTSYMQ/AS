library(DESeq2)
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(forcats)
library(org.Hs.eg.db)
library(pheatmap)
library(ggrepel)
library(limma)
try({
  dir_path = paste("../counts_processed/")
  out_comes <- data.frame()
  file_list <- list.files(path = dir_path)
  file_list <- gsub("\\.[^.]+$", "", file_list)
  length(file_list)
  file_list1 <- file_list[-length(file_list)]
  length(file_list1)
  file_list2 <- file_list[2:length(file_list)]
  for (i in file_list1){
    
    for (j in file_list2) {
      
      
      file_path1 <- paste(dir_path,i,".csv",sep = "")
      file_path2 <- paste(dir_path,j,".csv",sep = "")
      data1<- read.csv(file_path1,head=T,row.names=1)
      num_columnsda1 <- as.numeric(ncol(data1))
      
      data2<- read.csv(file_path2,head=T,row.names=1)
      num_columnsda2 <- as.numeric(ncol(data2))
      if (identical(data1, data2)) {
        message(paste(i, "is same with", j))
      } else {
        DES<- merge(data1,data2,by= "row.names")
        rownames(DES)<-DES[,1]
        DES=DES[,-1]
        DES<- na.omit(DES)
        column_names <- colnames(DES)
        col<- data.frame( x = column_names,
                          condition = factor(c(rep("control", num_columnsda1)
                                               ,rep("model",num_columnsda2)
                          )))
        rownames(col)<-col[,1]
        col <- data.frame(col)
        dds <- DESeqDataSetFromMatrix(DES, col, design= ~ condition)
        dds <- DESeq(dds)
        des <- DES[rowMeans(DES)>1,] 
        des <-des[complete.cases(des),]
        dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
        res <- results(dds1,contrast = c("condition",'model', 'control'))
        res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
        res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
        res1_up<- res1[which( res1$pvalue < 0.05),]     
        res1_down<- res1[which(res1$pvalue < 0.05),]   
        res1_total <- rbind(res1_up,res1_down)
        df4<-merge(des, res1, by = "row.names", all = T)
        df4 <- na.omit(df4)
        df3 <- df4
        folder_path <- paste(dir_path,i,"vs",j)
        dir.create(folder_path)
        write.csv(df3,file= paste(folder_path,"/",i,"vs",j,".csv",sep = ""))
        sig_gene <- df3 %>%
          dplyr::filter(pvalue<0.05)
        DEG<-sig_gene%>%
          dplyr::filter(abs(log2FoldChange)>1)
        DEG_path <- paste(folder_path,"/DEG.csv",sep = "")
        nDEGs<-df3%>%
          dplyr::filter(abs(log2FoldChange)<1|pvalue>0.05)
        nDEG_path <- paste(folder_path,"/nDEG.csv",sep = "")
        write.csv(DEG,DEG_path)
        write.csv(nDEGs,nDEG_path)
        colnames(sig_gene)[1] <- "gene"
        num_rows <- length(sig_gene[[1]])
        num_rows<- as.numeric(num_rows)
        num_rows2<- as.numeric(num_rows-1)
        ids <- bitr(DEG$Row.names,'SYMBOL','ENTREZID','org.Hs.eg.db') 
        gene <- ids$ENTREZID
        ego_ALL <- enrichGO(gene = gene,
                              OrgDb=org.Hs.eg.db,
                              keyType = "ENTREZID",
                              ont = "ALL",
                              minGSSize = 5,
                              maxGSSize = 1000,
                              pvalueCutoff = 0.05,
                              readable = TRUE) 
        ids2 <- bitr(nDEGs$Row.names,'SYMBOL','ENTREZID','org.Hs.eg.db') 
        gene2 <- ids2$ENTREZID
        ego_ALL2 <- enrichGO(gene = gene2,
                            OrgDb=org.Hs.eg.db,
                            keyType = "ENTREZID",
                            ont = "ALL",
                            minGSSize = 5,
                            maxGSSize = 1000,
                            pvalueCutoff = 0.05,
                            readable = TRUE) 
          write.csv(ego_ALL,"GO_degs.csv",row.names = T)
          write.csv(ego_ALL2,"GO_ndegs.csv",row.names = T)
           
        }
        }
  }
}
)
  