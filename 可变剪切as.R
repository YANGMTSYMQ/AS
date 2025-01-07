library(IsoformSwitchAnalyzeR, quietly=TRUE)
library(data.table)
library(ggplot2)
library(ggrepel)
library(data.table)
library(DESeq2)
library(clusterProfiler)
library(tidyverse)
library(forcats)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
isoformExpr <- importIsoformExpression(parentDir="../AS/isoforms_results/", pattern="isoforms.results")
myDesign <- data.frame(
  sampleID = colnames(isoformExpr$abundance)[-1],
  condition = gsub('-.*', '', colnames(isoformExpr$abundance)[-1])
)
sar1 <- importRdata(
  isoformCountMatrix   = isoformExpr$counts,
  isoformRepExpression = isoformExpr$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "/生信文件/annoation/人/ensem/Homo_sapiens.GRCh38.110.gtf",
  isoformNtFasta       = "/生信文件/annoation/人/ensem/RSEM_index.transcripts.fa",
  ignoreAfterPeriod = T ,
  ignoreAfterBar = T ,
  ignoreAfterSpace = T,
  fixStringTieAnnotationProblem = TRUE
)
sar1$isoformFeatures$gene_expr<-NA
gene_log2fc<-sar1$isoformFeatures$gene_log2_fold_change
for (i in 1:length(gene_log2fc)){
  if(sar1$isoformFeatures$gene_log2_fold_change[i]>1){
    sar1$isoformFeatures$gene_expr[i]<-"UP"
  }else if(sar1$isoformFeatures$gene_log2_fold_change[i]< -1){
    sar1$isoformFeatures$gene_expr[i]<-"down"
  }else{
    sar1$isoformFeatures$gene_expr[i]<-"No-Sig"
  }}
table(sar1$isoformFeatures$gene_expr)
sar2 <- preFilter(sar1)#过滤低质量数据
sar3 <- isoformSwitchTestDEXSeq(sar2,dIFcutoff = 0.1,alpha=0.05)#异构剪切分析,计算q值
sar32 <-isoformSwitchTestSatuRn(sar2,dIFcutoff = 0.1,alpha = 0.5,diagplots = F )#样本量多用这个
a<-unique(sar3$isoformFeatures$gene_id)
##############################################
sar4 <- analyzeORF(
  sar3,
  orfMethod = "longest",  #有four different methods ，详见官方文档
  showProgress=FALSE
)#注释ORF
sar5 <- extractSequence(sar4,
                        pathToOutput = '../AS/AANT/',
                        writeToFile=TRUE)#生成转录本核苷酸与对应氨基酸序列
sar5.1 <- analyzeCPAT(
  switchAnalyzeRlist   = sar4,
  pathToCPATresultFile = "../AS/AANT/output.txt",
  codingCutoff         = 0.725, # the coding potential cutoff we suggested for human
  removeNoncodinORFs   = TRUE   
)#根据外界转录本核苷酸进行比对
sar5.2 <- analyzePFAM(
                      switchAnalyzeRlist= sar5.1,
                          pathToPFAMresultFile = "../AS/AANT/PFAM.txt",
                      showProgress=FALSE
)#根据转录本氨基酸序列进行比对（结构域）
sar6 <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = sar5.2,
  dIFcutoff = 0,
  quiet=TRUE
)#注释可变剪切类型
cIn <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')
########################################
sar6<-extractSequence(sar6)
sar7 <- analyzeSwitchConsequences(
  sar6,
  consequencesToAnalyze = cIn, 
  dIFcutoff = 0.1, # very high cutoff for fast runtimes
  showProgress=FALSE
)
##########################
##################################
extractSwitchSummary(
  sar7,
  filterForConsequences = TRUE
) 
sar.ge_df=extractTopSwitches( 
  sar7, 
  filterForConsequences = FALSE, 
  n = NA,                
  extractGenes = TRUE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)
sar.iso_df=extractTopSwitches( 
  sar7, 
  filterForConsequences = FALSE, 
  n = NA,                
  extractGenes = FALSE,    # when FALSE isoforms are returned
  sortByQvals = TRUE
)
###############################################################
P=switchPlot(sar7, gene ='GPM6B',plotTopology=FALSE,dIFcutoff = 0.1,IFcutoff = 0.3,localTheme = theme_bw(base_size = 11))#localTheme = theme_bw(base_size = 11))
BBBB<-sar7$isoformFeatures
BBBB<-subset(BBBB,BBBB$gene_id%in%numberok2$Var1)
sar7$isoformFeatures<-BBBB
switchPlotTopSwitches(
  sar7,
  alpha = 0.05,
  dIFcutoff = 0.1,
  onlySigIsoforms = FALSE,
  n=1000,
  sortByQvals=TRUE,
  filterForConsequences = FALSE,
  pathToOutput = "../AS/iso_png0604///",
  splitComparison=TRUE,
  splitFunctionalConsequences = TRUE,
  IFcutoff=0.1,
  fileType = "png",
  additionalArguments=list(),
  quiet=FALSE
)
#######################
targetgene<-as.data.frame(unique(nDEG_sar7$gene_id))
write.csv(targetgene,"../../targetgene1.csv")
geneexpression<-extractGeneExpression(sar7)
sar7_backup<-sar7
extractConsequenceSummary(
  sar7,
  consequencesToAnalyze='all',  #也可以指定感兴趣的consequences
  plotGenes = T,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = T,      # enables analysis of fraction of significant features
)

extractSplicingGenomeWide(
  sar7,
  featureToExtract = 'isoformUsage',                 # all isoforms stored in the switchAnalyzeRlist
  splicingToAnalyze = c('A3','A5','ES','IR'), # Splice types significantly enriched in COAD
  plot=T,
  returnResult=F  # Preventing the summary statistics to be returned as a data.frame
)
extractConsequenceEnrichment(
  sar7,
  consequencesToAnalyze=c('ORF_seq_similarity','domains_identified','domain_isotype','ORF_length'),
  minEventsForPlotting= F,
  analysisOppositeConsequence = TRUE,
  returnResult = T # if TRUE returns a data.frame with the results
)
extractSwitchSummary(
  sar7,
  filterForConsequences = TRUE
) 
extractSplicingEnrichment(
  sar7,
  splicingToAnalyze = c('A3','MES','ATSS','ATTS','A5','ES','IR'), 
  returnResult = TRUE # if TRUE returns a data.frame with the results
)


AAAA=sar7$isoformFeatures
AAAA<-subset(AAAA,AAAA$isoform_switch_q_value<0.05)
ggplot(AAAA, aes(x=-dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  scale_color_manual('Signficant\nIsoform Switch', values = c('grey','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #geom_label_repel(
    #data = subset(sar7$isoformFeatures, (sar7$isoformFeatures$isoform_switch_q_value) < 0.00001 & abs(sar7$isoformFeatures$dIF) >= 0.4),
    #aes(label = isoform_id),
    #size = 3,
   # box.padding = unit(0.5, "lines"),
    #point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE,
    #)
##############################iso的火山图
ggplot(data=sar7$isoformFeatures, aes(x=-iso_log2_fold_change, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(iso_log2_fold_change) > 1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-1, 1), linetype='dashed') + # default cutoff
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='log2FC', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
###############################################limma_iso火山图




