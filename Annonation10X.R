library(dplyr)
library(tibble)
Annoation10X<-function (MarvelObject) 
{
  MarvelObject <- MarvelObject
  gtf <- MarvelObject$gtf
  df.sj.count <- MarvelObject$sj.count.matrix
  message("Creating splice junction metadata...")
  df <- data.frame(coord.intron = rownames(df.sj.count), stringsAsFactors = FALSE)
  . <- strsplit(rownames(df.sj.count), split = ":", fixed = TRUE)
  df$chr <- sapply(., function(x) {
    as.character(x[1])
  })
  df$start <- sapply(., function(x) {
    as.character(x[2])
  })
  df$end <- sapply(., function(x) {
    as.character(x[3])
  })
  message("Parsing GTF...")
  gtf <- gtf[which(gtf$V3 == "exon"), ]
  . <- strsplit(gtf$V9, split = ";")
  . <- sapply(., function(x) grep("gene_name", x, value = TRUE))
  . <- gsub("gene_name", "", .)
  . <- gsub(" ", "", .)
  . <- gsub("\"", "", .)
  gtf$gene_short_name <- .
  gtf$V4 <- gtf$V4 - 1
  gtf$V5 <- gtf$V5 + 1
  message("Matching gene names with SJ start coordinates in GTF...")
  gtf.small <- gtf[, c("V1", "V5", "gene_short_name")]
  gtf.small$chr.pos <- paste(gtf$V1, ":", gtf$V5, sep = "")
  #gtf.small$chr.pos <- paste("chr", gtf$V1, ":", gtf$V5, sep = "")
  gtf.small <- gtf.small[, c("chr.pos", "gene_short_name")]
  gtf.small <- unique(gtf.small)
  gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
  . <- by(gtf.small[, "gene_short_name"], gtf.small[, "chr.pos"], 
          function(x) {
            paste(x, collapse = "|")
          })
  gtf.small.collapsed <- data.frame(chr = as.character(names(.)), 
                                    gene_short_name = as.character(.), stringsAsFactors = FALSE)
  gtf.small.collapsed.start <- gtf.small.collapsed
  names(gtf.small.collapsed.start) <- paste(names(gtf.small.collapsed.start), 
                                            ".start", sep = "")
  message("Matching gene names with SJ end coordinates in GTF...")
  gtf.small <- gtf[, c("V1", "V4", "gene_short_name")]
  #gtf.small$chr.pos <- paste("chr", gtf$V1, ":", gtf$V4, sep = "")
  gtf.small$chr.pos <- paste(gtf$V1, ":", gtf$V4, sep = "")
  gtf.small <- gtf.small[, c("chr.pos", "gene_short_name")]
  gtf.small <- unique(gtf.small)
  gtf.small$chr.pos <- as.factor(gtf.small$chr.pos)
  . <- by(gtf.small[, "gene_short_name"], gtf.small[, "chr.pos"], 
          function(x) {
            paste(x, collapse = "|")
          })
  gtf.small.collapsed <- data.frame(chr = as.character(names(.)), 
                                    gene_short_name = as.character(.), stringsAsFactors = FALSE)
  gtf.small.collapsed.end <- gtf.small.collapsed
  names(gtf.small.collapsed.end) <- paste(names(gtf.small.collapsed.end), 
                                          ".end", sep = "")
  message("Annotating splice junctions...")
  df$chr.start <- paste(df$chr, ":", df$start, sep = "")
  df$chr.end <- paste(df$chr, ":", df$end, sep = "")
  df$id <- c(1:nrow(df))
  df <- left_join(df, gtf.small.collapsed.start, by = "chr.start")
  df <- left_join(df, gtf.small.collapsed.end, by = "chr.end")
  df$sj.type <- NA
  index.l <- !grepl("|", df$gene_short_name.start, fixed = TRUE) & 
    !grepl("|", df$gene_short_name.end, fixed = TRUE) & 
    !is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & 
    df$gene_short_name.start == df$gene_short_name.end
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_known.single.gene|end_known.single.gene|same"
  }
  index.l <- is.na(df$sj.type) & !grepl("|", df$gene_short_name.start, 
                                        fixed = TRUE) & !grepl("|", df$gene_short_name.end, 
                                                               fixed = TRUE) & !is.na(df$gene_short_name.start) & !is.na(df$gene_short_name.end) & 
    !df$gene_short_name.start == df$gene_short_name.end
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_known.single.gene|end_known.single.gene|different"
  }
  index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.start) & 
    is.na(df$gene_short_name.end)
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_unknown.gene|end_unknown.gene"
  }
  index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.start) & 
    !is.na(df$gene_short_name.end) & !grepl("|", df$gene_short_name.end, 
                                            fixed = TRUE)
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_unknown.gene|end_known.single.gene"
  }
  index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.start) & 
    !is.na(df$gene_short_name.end) & grepl("|", df$gene_short_name.end, 
                                           fixed = TRUE)
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_unknown.gene|end_known.multi.gene"
  }
  index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.end) & 
    !is.na(df$gene_short_name.start) & !grepl("|", df$gene_short_name.start, 
                                              fixed = TRUE)
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_known.single.gene|end_unknown.gene"
  }
  index.l <- is.na(df$sj.type) & is.na(df$gene_short_name.end) & 
    !is.na(df$gene_short_name.start) & grepl("|", df$gene_short_name.start, 
                                             fixed = TRUE)
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_known.multi.gene|end_unknown.gene"
  }
  index.l <- is.na(df$sj.type) & !is.na(df$gene_short_name.end) & 
    !is.na(df$gene_short_name.start) & grepl("|", df$gene_short_name.start, 
                                             fixed = TRUE) & grepl("|", df$gene_short_name.end, fixed = TRUE)
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_known.multi.gene|end_known.multi.gene"
  }
  index.l <- is.na(df$sj.type) & !is.na(df$gene_short_name.end) & 
    !is.na(df$gene_short_name.start) & grepl("|", df$gene_short_name.start, 
                                             fixed = TRUE) & !grepl("|", df$gene_short_name.end, 
                                                                    fixed = TRUE)
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_known.multi.gene|start_known.single.gene"
  }
  index.l <- is.na(df$sj.type) & !is.na(df$gene_short_name.end) & 
    !is.na(df$gene_short_name.start) & !grepl("|", df$gene_short_name.start, 
                                              fixed = TRUE) & grepl("|", df$gene_short_name.end, fixed = TRUE)
  index <- which(index.l == TRUE)
  if (length(index) >= 1) {
    df$sj.type[index] <- "start_known.single.gene|end_known.multi.gene"
  }
  if (sum(is.na(df$sj.type)) == 0) {
    message("All SJ successfully annotated")
  }
  else {
    message("Some SJ NOT successfully annotated")
  }
  df <- df[, c("coord.intron", "gene_short_name.start", "gene_short_name.end", 
               "sj.type")]
  MarvelObject$sj.metadata <- df
  return(MarvelObject)
}

