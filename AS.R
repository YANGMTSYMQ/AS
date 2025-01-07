library(IsoformSwitchAnalyzeR, quietly=TRUE)
library(data.table)
library(ggplot2)
library(ggrepel)
library(data.table)
library(DESeq2)
library(clusterProfiler)
library(tidyverse)
library(forcats)
library(org.Hs.eg.db)
isoformExpr <- importIsoformExpression(parentDir="../AS/isoforms_results/", pattern="isoforms.results")
myDesign <- data.frame(
  sampleID = colnames(isoformExpr$abundance)[-1],
  condition = gsub('_.*', '', colnames(isoformExpr$abundance)[-1])
)
STEP1 <- importRdata(
  isoformCountMatrix   = isoformExpr$counts,
  isoformRepExpression = isoformExpr$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "/Bioinfo//annoation/human//ensem/Homo_sapiens.GRCh38.110.gtf",
  isoformNtFasta       = "/Bioinfo//annoation/human//ensem/RSEM_index.transcripts.fa",
  ignoreAfterPeriod = T ,
  ignoreAfterBar = T ,
  ignoreAfterSpace = T,
  fixStringTieAnnotationProblem = TRUE
)
STEP2 <- preFilter(STEP1)
STEP3 <- isoformSwitchTestSatuRn(STEP2)
STEP4 <- analyzeORF(
  STEP3,
  orfMethod = "longest",  
  showProgress=FALSE
)
STEP5 <- extractSequence(STEP4,
                         pathToOutput = '../AS/AANT/',
                         writeToFile=TRUE)

cpat.py -g isoformSwitchAnalyzeR_isoform_nt.fasta \
-d /home/liusai/altersplice/software/CPAT/CPAT-3.0.0/dat/Human_logitModel.RData \
-x /home/liusai/altersplice/software/CPAT/CPAT-3.0.0/dat/Human_Hexamer.tsv \
-o output.txt

pfam_scan.pl -fasta isoformSwitchAnalyzeR_isoform_AA.fasta -dir /home/liusai/altersplice/software/PFAM/ -outfile PFAM.txt -as


STEP5.1 <- analyzeCPAT(
  switchAnalyzeRlist   = STEP4,
  pathToCPATresultFile = "../AS/AANT/output.txt",
  codingCutoff         = 0.725, 
  removeNoncodinORFs   = TRUE   
)
STEP5.2 <- analyzePFAM(
  switchAnalyzeRlist= STEP5.1,
  pathToPFAMresultFile = "../AS/AANT/PFAM.txt",
  showProgress=FALSE
)
STEP6 <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = STEP5.2,
  dIFcutoff = 0,
  quiet=TRUE
)
cIn <- c('intron_retention','coding_potential','NMD_status','domains_identified','ORF_seq_similarity')
STEP6<-extractSequence(STEP6)
STEP7 <- analyzeSwitchConsequences(
  STEP6,
  consequencesToAnalyze = cIn, 
  dIFcutoff = 0.1,
  showProgress=FALSE
)
extractSwitchSummary(
  STEP7,
  filterForConsequences = TRUE
) 
STEP.ge_df=extractTopSwitches( 
  STEP7, 
  filterForConsequences = FALSE, 
  n = NA,                
  extractGenes = TRUE,   
  sortByQvals = TRUE
)
STEP.iso_df=extractTopSwitches( 
  STEP7, 
  filterForConsequences = FALSE, 
  n = NA,                
  extractGenes = FALSE, 
  sortByQvals = TRUE
)
STEP8<-STEP7
INFO<-STEP8$isoformFeatures
INFO<-subset(INFO,INFO$gene_id%in%numberok2$Var1)
STEP8$isoformFeatures<-INFO
switchPlotTopSwitches(
  STEP8,
  alpha = 0.05,
  dIFcutoff = 0.1,
  onlySigIsoforms = FALSE,
  n=1000,
  sortByQvals=TRUE,
  filterForConsequences = FALSE,
  pathToOutput = "../AS/iso_png///",
  splitComparison=TRUE,
  splitFunctionalConsequences = TRUE,
  IFcutoff=0.1,
  fileType = "png",
  additionalArguments=list(),
  quiet=FALSE
)
geneexpression<-extractGeneExpression(STEP7)
extractConsequenceSummary(
  STEP7,
  consequencesToAnalyze='all',  
  plotGenes = T,           
  asFractionTotal = T,      
)

extractSplicingGenomeWide(
  STEP7,
  featureToExtract = 'isoformUsage',
  splicingToAnalyze = c("A3","A5","ES","IR","MEE"), 
  plot=T,
  returnResult=F  
)
extractConsequenceEnrichment(
  STEP7,
  consequencesToAnalyze=c('ORF_seq_similarity','domains_identified','domain_isotype','ORF_length'),
  minEventsForPlotting= F,
  analysisOppositeConsequence = TRUE,
  returnResult = T 
)
extractSwitchSummary(
  STEP7,
  filterForConsequences = TRUE
) 
extractSplicingEnrichment(
  STEP7,
  splicingToAnalyze = c('A3','MES','ATSS','ATTS','A5','ES','IR'), 
  returnResult = TRUE 
)




