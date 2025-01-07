library(MARVEL)
library(Seurat)
library(ShinyCell)
library(shiny)
library(foreign)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(wiggleplotr)
library(Seurat)
library(R.utils)
library(Matrix)
library(magrittr)
library(harmony)
dir.create("../work")
dir.create("../work/star")
dir.create("../work/gene")
dir.create("../work/gene/raw")
dir.create("../work/gene/norm")
file_list<-list.files("../work/star/")

####################################################
dir.create("../work/sj")
dir.create("../work/sj/raw")

# 整理SJ meta
samples<-list.files("../work/star/")
for(x in samples){
  sj_name = data.table::fread(paste0("../work/star/",x,"/Solo.out/SJ.out.tab"), data.table=F) %>% 
    dplyr::mutate(coord.intron=paste(V1, V2, V3, sep=":")) %>% 
    dplyr::select(coord.intron)
  write.table(sj_name, sep="\t", row.names=F, col.names=F, quote=F,
              file=paste0("../work/star/",x,"/Solo.out/SJ/raw/features.tsv"))
}

sce_sj.list = lapply(samples, function(x){
  print(x)
  fls = list.files(paste0("../work/star/",x,"/Solo.out/SJ/raw"), full.name=T)
  if(length(grep("gz",fls))!=3){
    for(fl in fls){ gzip(fl, remove=FALSE) }
  }
  count = Read10X(paste0("../work/star/",x,"/Solo.out/SJ/raw"),gene.column=1)
  sce_sj = CreateSeuratObject(count)
  sce_sj = RenameCells(sce_sj, add.cell.id=x)
  sce_sj = sce_sj[,colnames(sce.list[[x]])]
  return(sce_sj)
})

sce_sj = merge(sce_sj.list[[1]], sce_sj.list[-1])
sce_sj = JoinLayers(sce_sj)

sj_count_merge = sce_sj@assays$RNA@layers$counts
writeMM(sj_count_merge, file = '../work/sj/raw/matrix.mtx')

sj_count_merge_phe = sce_sj@meta.data %>% 
  tibble::rownames_to_column("cell.id") %>% 
  dplyr::select(cell.id)  # 仅需要1列
write.table(sj_count_merge_phe, sep="\t", row.names=F, 
            file="../work/sj/raw/phenoData.txt")

sj_count_merge_feat = data.frame(coord.intron=rownames(sce_sj)) %>% 
  dplyr::mutate(coord.intron=gsub("-","_", coord.intron))
write.table(sj_count_merge_feat, sep="\t", row.names=F, 
            file="../work/sj/raw/featureData.txt")
dir.create("../work/umap/test/")

sce = sce %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% 
  Seurat::RunPCA() %>% 
  Seurat::RunUMAP(dims = 1:30)

umap_embed = sce@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("cell.id") %>%
  dplyr::rename(x=umap_1, y=umap_2) 
write.table(umap_embed, row.names=F, sep="\t",
            file="../work/umap/UMAP_coordinate.txt")
df.gene.norm <- readMM("../work/gene/norm/matrix_normalised.mtx")
df.gene.norm.pheno <- read.table("../work/gene/norm/phenoData.txt", sep="\t", header=TRUE)
df.gene.norm.feature <- read.table("../work/gene/norm/featureData.txt", sep="\t", header=TRUE)
df.gene.count <- readMM("../work/gene/raw/matrix.mtx")
df.gene.count.pheno <- read.table("../work/gene/raw/phenoData.txt", sep="\t", header=TRUE)
df.gene.count.feature <- read.table("../work/gene/raw/featureData.txt", sep="\t", header=TRUE)
df.sj.count <- readMM("../work/sj/raw/matrix.mtx")
df.sj.count.pheno <- read.table("../work/sj/raw/phenoData.txt", sep="\t", header=TRUE)
df.sj.count.feature  <- read.table("../work/sj/raw/featureData.txt", sep="\t", header=TRUE)

rownames(df.sj.count) = df.sj.count.feature[,1]
colnames(df.sj.count) = df.sj.count.pheno[,1]

# (4) 细胞降维
df.coord <- read.table("../work/umap/UMAP_coordinate.txt", sep="\t", header=TRUE)

# (5) GTF
gtf <- data.table::fread("../../../annoation/mice/ensem/Mus_musculus.GRCm39.110.gtf", sep="\t", header=FALSE) %>% 
  as.data.frame() 
## 由于不同来源gtf版本差异而踩的坑，需要修改如下
gtf$V9 = gsub("gene_type","gene_biotype",gtf$V9)
gtf$V1 = gsub("chr","",gtf$V1)
gtf <- gtf %>%
  mutate(gene_short_name = str_extract(V9, '(?<=gene_name ")[^"]+'))

library(dplyr)
library(stringr)
library(MARVEL)
marvel <- CreateMarvelObject.10x(gene.norm.matrix=df.gene.norm,
                                 gene.norm.pheno=df.gene.norm.pheno,
                                 gene.norm.feature=df.gene.norm.feature,
                                 gene.count.matrix=df.gene.count,
                                 gene.count.pheno=df.gene.count.pheno,
                                 gene.count.feature=df.gene.count.feature,
                                 sj.count.matrix=df.sj.count,
                                 sj.count.pheno=df.sj.count.pheno,
                                 sj.count.feature=df.sj.count.feature,
                                 pca=df.coord,
                                 gtf=gtf
)
Aa<-marvel$sample.metadata
Aa$Group1<-sub("_.*", "", Aa$Group)
marvel$sample.metadata<-Aa
marvel <- AnnotateGenes.10x(MarvelObject=marvel)
marvel <- Annoation10X(MarvelObject=marvel)
marvelback<-marvel
marvel <- ValidateSJ.10x(MarvelObject=marvel)
marvel <- FilterGenes.10x(MarvelObject=marvel,
                          gene.type="protein_coding")
cell.group.list = split(marvel$sample.metadata$cell.id,marvel$sample.metadata$Group1)
## 基因在细胞群的表达率
#######################################保存中间结果
marvel <- CheckAlignment.10x(MarvelObject=marvel)
###########################################
saveRDS(marvel,"marvel.rds")
## 首先分析差异SJ
marvel <- comparevalues.sj(MarvelObject=marvel,
                           cell.group.g1=cell.group.list[[1]],
                           cell.group.g2=cell.group.list[[2]],
                           min.pct.cells.genes=10, #至少10%细胞表达
                           min.pct.cells.sj=10,    #至少10%细胞剪切
                           min.gene.norm=0.1,
                           seed=1,
                           n.iterations=100,
                           downsample=F,  #根据样本量少的一组下采样
                           show.progress=FALSE
)
## 然后进一步分析差异基因
marvel <- CompareValues.Genes.10x(MarvelObject=marvel,
                                  show.progress=FALSE
)

head(marvel$DE$SJ$Table)

marvel <- IsoSwitch.10x(MarvelObject=marvel,
                        pval.sj=0.05,
                        delta.sj=5,
                        min.gene.norm=0.1,
                        pval.adj.gene=0.05,
                        log2fc.gene=0.1
)



