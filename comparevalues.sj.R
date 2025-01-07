comparevalues.sj<-function (MarvelObject, coord.introns = NULL, cell.group.g1, 
          cell.group.g2, min.pct.cells.genes = 10, min.pct.cells.sj = 10, 
          min.gene.norm = 1, seed = 1, n.iterations = 100, downsample = FALSE, 
          show.progress = TRUE) 
{
  MarvelObject <- MarvelObject
  sample.metadata <- MarvelObject$sample.metadata
  sj.metadata <- MarvelObject$sj.metadata
  df.gene.norm <- MarvelObject$gene.norm.matrix
  df.gene.count <- MarvelObject$gene.count.matrix
  df.sj.count <- MarvelObject$sj.count.matrix
  coord.introns.custom <- coord.introns
  cell.group.g1 <- cell.group.g1
  cell.group.g2 <- cell.group.g2
  min.pct.cells.genes <- min.pct.cells.genes
  min.pct.cells.sj <- min.pct.cells.sj
  seed <- seed
  n.iterations <- n.iterations
  downsample <- downsample
  show.progress <- show.progress
  min.gene.norm <- min.gene.norm
  if (!is.null(coord.introns.custom[1])) {
    overlap <- intersect(coord.introns.custom, row.names(df.sj.count))
    df.sj.count <- df.sj.count[overlap, ]
    sj.metadata <- sj.metadata[which(sj.metadata$coord.intron %in% 
                                       overlap), ]
    gene_short_names <- unique(sj.metadata$gene_short_name.start)
    df.gene.norm <- df.gene.norm[gene_short_names, ]
    df.gene.count <- df.gene.count[gene_short_names, ]
    message(paste(length(coord.introns.custom), " SJs specified by user", 
                  sep = ""))
    message(paste(length(overlap), " overlapping SJs found and subset-ed", 
                  sep = ""))
    message(paste(length(gene_short_names), " corresponding genes found and subset-ed", 
                  sep = ""))
  }
  if (downsample == TRUE) {
    set.seed(seed)
    n.cells.downsample <- min(length(cell.group.g1), length(cell.group.g2))
    cell.group.g1 <- sample(cell.group.g1, size = n.cells.downsample, 
                            replace = FALSE)
    cell.group.g2 <- sample(cell.group.g2, size = n.cells.downsample, 
                            replace = FALSE)
  }
  message(paste(length(cell.group.g1), " cells from Group 1 and ", 
                length(cell.group.g2), " cells from Group 2 included", 
                sep = ""))
  df.gene.norm <- df.gene.norm[, c(cell.group.g1, cell.group.g2)]
  df.gene.count <- df.gene.count[, c(cell.group.g1, cell.group.g2)]
  df.sj.count <- df.sj.count[, c(cell.group.g1, cell.group.g2)]
  table(colnames(df.gene.norm) == colnames(df.gene.count))
  table(colnames(df.gene.count) == colnames(df.sj.count))
  df.gene.norm.small <- df.gene.norm[, cell.group.g1]
  . <- apply(df.gene.norm.small, 1, function(x) {
    sum(x != 0)
  })
  . <- data.frame(cell.group = "cell.group.g1", gene_short_name = names(.), 
                  n.cells.total = length(cell.group.g1), n.cells.expr = as.numeric(.), 
                  pct.cells.expr = round(as.numeric(.)/length(cell.group.g1) * 
                                           100, digits = 2), stringsAsFactors = FALSE)
  results.g1 <- .
  df.gene.norm.small <- df.gene.norm[, cell.group.g2]
  . <- apply(df.gene.norm.small, 1, function(x) {
    sum(x != 0)
  })
  . <- data.frame(cell.group = "cell.group.g2", gene_short_name = names(.), 
                  n.cells.total = length(cell.group.g2), n.cells.expr = as.numeric(.), 
                  pct.cells.expr = round(as.numeric(.)/length(cell.group.g2) * 
                                           100, digits = 2), stringsAsFactors = FALSE)
  results.g2 <- .
  results <- rbind.data.frame(results.g1, results.g2)
  results$cell.group <- factor(results$cell.group, levels = c("cell.group.g1", 
                                                              "cell.group.g2"))
  results <- results[which(results$pct.cells.expr > min.pct.cells.genes), 
  ]
  gene_short.names.1 <- results[which(results$cell.group == 
                                        "cell.group.g1"), "gene_short_name"]
  gene_short.names.2 <- results[which(results$cell.group == 
                                        "cell.group.g2"), "gene_short_name"]
  gene_short_names <- intersect(gene_short.names.1, gene_short.names.2)
  message(paste(length(gene_short.names.1), " genes expressed in cell group 1", 
                sep = ""))
  message(paste(length(gene_short.names.2), " genes expressed in cell group 2", 
                sep = ""))
  message(paste(length(gene_short_names), " genes expressed in BOTH cell group and retained", 
                sep = ""))
  df.gene.norm.small <- df.gene.norm[gene_short_names, ]
  . <- apply(df.gene.norm.small, 1, function(x) {
    mean(log2(x + 1))
  })
  mean.combined.df <- data.frame(gene_short_name = names(.), 
                                 mean.expr.gene.norm.g1.g2 = as.numeric(.), stringsAsFactors = FALSE)
  index <- which(mean.combined.df$mean.expr.gene.norm.g1.g2 > 
                   min.gene.norm)
  mean.combined.df <- mean.combined.df[index, ]
  gene_short_names <- mean.combined.df$gene_short_name
  message(paste(length(gene_short_names), " genes with mean log2(expression + 1) > ", 
                min.gene.norm, " retained", sep = ""))
  df.sj.count.small <- df.sj.count[, cell.group.g1]
  coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% 
                                       gene_short_names), "coord.intron"]
  length(coord.introns)
  df.sj.count.small <- df.sj.count.small[coord.introns, ]
  . <- apply(df.sj.count.small, 1, function(x) {
    sum(x != 0)
  })
  . <- data.frame(cell.group = "cell.group.g1", coord.intron = names(.), 
                  n.cells.total = length(cell.group.g1), n.cells.expr = as.numeric(.), 
                  pct.cells.expr = round(as.numeric(.)/length(cell.group.g1) * 
                                           100, digits = 2), stringsAsFactors = FALSE)
  results.g1 <- .
  df.sj.count.small <- df.sj.count[, cell.group.g2]
  coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% 
                                       gene_short_names), "coord.intron"]
  length(coord.introns)
  df.sj.count.small <- df.sj.count.small[coord.introns, ]
  . <- apply(df.sj.count.small, 1, function(x) {
    sum(x != 0)
  })
  . <- data.frame(cell.group = "cell.group.g2", coord.intron = names(.), 
                  n.cells.total = length(cell.group.g2), n.cells.expr = as.numeric(.), 
                  pct.cells.expr = round(as.numeric(.)/length(cell.group.g2) * 
                                           100, digits = 2), stringsAsFactors = FALSE)
  results.g2 <- .
  results <- rbind.data.frame(results.g1, results.g2)
  results$cell.group <- factor(results$cell.group, levels = c("cell.group.g1", 
                                                              "cell.group.g2"))
  results <- results[which(results$pct.cells.expr > min.pct.cells.sj), 
  ]
  coord.introns.1 <- results[which(results$cell.group == "cell.group.g1"), 
                             "coord.intron"]
  coord.introns.2 <- results[which(results$cell.group == "cell.group.g2"), 
                             "coord.intron"]
  coord.introns <- unique(c(coord.introns.1, coord.introns.2))
  message(paste(length(coord.introns.1), " SJ expressed in cell group 1", 
                sep = ""))
  message(paste(length(coord.introns.2), " SJ expressed in cell group 2", 
                sep = ""))
  message(paste(length(coord.introns), " SJ expressed in EITHER cell groups and retained", 
                sep = ""))
  n.sj <- length(coord.introns)
  n.genes <- length(unique(sj.metadata[which(sj.metadata$coord.intron %in% 
                                               coord.introns), "gene_short_name.start"]))
  message(paste("Total of ", n.sj, " SJ from ", n.genes, " genes included for DE analysis", 
                sep = ""))
  sj.metadata <- sj.metadata[which(sj.metadata$coord.intron %in% 
                                     coord.introns), ]
  df.sj.count.g1 <- df.sj.count[sj.metadata$coord.intron, 
                                cell.group.g1]
  df.gene.count.g1 <- df.gene.count[unique(sj.metadata$gene_short_name.start), 
                                    cell.group.g1]
  df.sj.count.g2 <- df.sj.count[sj.metadata$coord.intron, 
                                cell.group.g2]
  df.gene.count.g2 <- df.gene.count[unique(sj.metadata$gene_short_name.start), 
                                    cell.group.g2]
  table(row.names(df.sj.count.g1) == row.names(df.sj.count.g2))
  table(row.names(df.gene.count.g1) == row.names(df.gene.count.g2))
  message("Computing PSI for cell group 1...")
  n.cells.total <- ncol(df.sj.count.g1)
  n.cells.expr.sj <- apply(df.sj.count.g1, 1, function(x) {
    sum(x != 0)
  })
  pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 
                               100, digits = 2)
  sj.count.total <- apply(df.sj.count.g1, 1, function(x) {
    sum(x)
  })
  results.sj <- data.frame(coord.intron = names(sj.count.total), 
                           n.cells.total = n.cells.total, n.cells.expr.sj = n.cells.expr.sj, 
                           pct.cells.expr.sj = pct.cells.expr.sj, sj.count.total = sj.count.total, 
                           stringsAsFactors = FALSE)
  row.names(results.sj) <- NULL
  results.sj <- left_join(results.sj, sj.metadata[, c("coord.intron", 
                                                 "gene_short_name.start")], by = "coord.intron")
  n.cells.expr.gene <- apply(df.gene.count.g1, 1, function(x) {
    sum(x != 0)
  })
  pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 
                                 100, digits = 2)
  gene.count.total <- apply(df.gene.count.g1, 1, function(x) {
    sum(x)
  })
  results.gene <- data.frame(gene_short_name.start = names(gene.count.total), 
                             n.cells.expr.gene = n.cells.expr.gene, pct.cells.expr.gene = pct.cells.expr.gene, 
                             gene.count.total = gene.count.total, stringsAsFactors = FALSE)
  results <- left_join(results.sj, results.gene, by = "gene_short_name.start")
  results$psi <- round(results$sj.count.total/results$gene.count.total * 
                         100, digits = 2)
  cols <- c("coord.intron", "gene_short_name.start", "n.cells.total", 
            "n.cells.expr.sj", "pct.cells.expr.sj", "n.cells.expr.gene", 
            "pct.cells.expr.gene", "sj.count.total", "gene.count.total", 
            "psi")
  results <- results[, cols]
  names(results)[which(names(results) == "gene_short_name.start")] <- "gene_short_name"
  names(results)[-which(names(results) %in% c("coord.intron", 
                                              "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% 
                                                                                                    c("coord.intron", "gene_short_name"))], ".g1", sep = "")
  results.g1 <- results
  message("Computing PSI for cell group 2...")
  n.cells.total <- ncol(df.sj.count.g2)
  n.cells.expr.sj <- apply(df.sj.count.g2, 1, function(x) {
    sum(x != 0)
  })
  pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 
                               100, digits = 2)
  sj.count.total <- apply(df.sj.count.g2, 1, function(x) {
    sum(x)
  })
  results.sj <- data.frame(coord.intron = names(sj.count.total), 
                           n.cells.total = n.cells.total, n.cells.expr.sj = n.cells.expr.sj, 
                           pct.cells.expr.sj = pct.cells.expr.sj, sj.count.total = sj.count.total, 
                           stringsAsFactors = FALSE)
  row.names(results.sj) <- NULL
  results.sj <- left_join(results.sj, sj.metadata[, c("coord.intron", 
                                                 "gene_short_name.start")], by = "coord.intron")
  n.cells.expr.gene <- apply(df.gene.count.g2, 1, function(x) {
    sum(x != 0)
  })
  pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 
                                 100, digits = 2)
  gene.count.total <- apply(df.gene.count.g2, 1, function(x) {
    sum(x)
  })
  results.gene <- data.frame(gene_short_name.start = names(gene.count.total), 
                             n.cells.expr.gene = n.cells.expr.gene, pct.cells.expr.gene = pct.cells.expr.gene, 
                             gene.count.total = gene.count.total, stringsAsFactors = FALSE)
  results <- left_join(results.sj, results.gene, by = "gene_short_name.start")
  results$psi <- round(results$sj.count.total/results$gene.count.total * 
                         100, digits = 2)
  cols <- c("coord.intron", "gene_short_name.start", "n.cells.total", 
            "n.cells.expr.sj", "pct.cells.expr.sj", "n.cells.expr.gene", 
            "pct.cells.expr.gene", "sj.count.total", "gene.count.total", 
            "psi")
  results <- results[, cols]
  names(results)[which(names(results) == "gene_short_name.start")] <- "gene_short_name"
  names(results)[-which(names(results) %in% c("coord.intron", 
                                              "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% 
                                                                                                    c("coord.intron", "gene_short_name"))], ".g2", sep = "")
  results.g2 <- results
  table(results.g1$coord.intron == results.g2$coord.intron)
  table(results.g1$gene_short_name == results.g2$gene_short_name)
  index.l <- table(results.g1$coord.intron == results.g2$coord.intron)
  index.true <- length(which(names(index.l) == TRUE))
  index.false <- length(which(names(index.l) == FALSE))
  if (index.true == 1 & index.false == 0) {
    results.g2$coord.intron <- NULL
    results.g2$gene_short_name <- NULL
    results <- cbind.data.frame(results.g1, results.g2)
  }
  else {
    message("Error in merging tables from Group 1 and 2")
  }
  results$log2fc <- log2((results$psi.g2 + 1)/(results$psi.g1 + 
                                                 1))
  results$delta <- results$psi.g2 - results$psi.g1
  results.obs <- results
  message("Creating null distributions...")
  set.seed(seed)
  if (show.progress == TRUE) {
    pb <- txtProgressBar(1, n.iterations, style = 3)
  }
  .list.results.perm <- list()
  for (i in 1:n.iterations) {
    cell.ids.shuffled <- sample(colnames(df.sj.count), size = ncol(df.sj.count), 
                                replace = FALSE)
    df.sj.count.shuffled <- df.sj.count
    colnames(df.sj.count.shuffled) <- cell.ids.shuffled
    df.gene.count.shuffled <- df.gene.count
    colnames(df.gene.count.shuffled) <- cell.ids.shuffled
    table(colnames(df.sj.count.shuffled) == colnames(df.gene.count.shuffled))
    index.l <- table(colnames(df.sj.count.shuffled) == colnames(df.gene.count.shuffled))
    index.true <- length(which(names(index.l) == TRUE))
    index.false <- length(which(names(index.l) == FALSE))
    if (index.true == 1 & index.false == 0) {
    }
    else {
      return(message(paste("Error in iteration ", 1, " ...", 
                           sep = "")))
    }
    df.sj.count.g1 <- df.sj.count.shuffled[results.obs$coord.intron, 
                                           cell.group.g1]
    df.gene.count.g1 <- df.gene.count.shuffled[unique(results.obs$gene_short_name), 
                                               cell.group.g1]
    df.sj.count.g2 <- df.sj.count.shuffled[results.obs$coord.intron, 
                                           cell.group.g2]
    df.gene.count.g2 <- df.gene.count.shuffled[unique(results.obs$gene_short_name), 
                                               cell.group.g2]
    table(row.names(df.sj.count.g1) == row.names(df.sj.count.g2))
    table(row.names(df.gene.count.g1) == row.names(df.gene.count.g2))
    sj.count.total <- apply(df.sj.count.g1, 1, function(x) {
      sum(x)
    })
    results.sj <- data.frame(coord.intron = names(sj.count.total), 
                             sj.count.total = sj.count.total, stringsAsFactors = FALSE)
    row.names(results.sj) <- NULL
    results.sj <- left_join(results.sj, sj.metadata[, c("coord.intron", 
                                                   "gene_short_name.start")], by = "coord.intron")
    gene.count.total <- apply(df.gene.count.g1, 1, function(x) {
      sum(x)
    })
    results.gene <- data.frame(gene_short_name.start = names(gene.count.total), 
                               gene.count.total = gene.count.total, stringsAsFactors = FALSE)
    results <- left_join(results.sj, results.gene, by = "gene_short_name.start")
    results$psi <- round(results$sj.count.total/results$gene.count.total * 
                           100, digits = 2)
    names(results)[which(names(results) == "gene_short_name.start")] <- "gene_short_name"
    names(results)[-which(names(results) %in% c("coord.intron", 
                                                "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% 
                                                                                                      c("coord.intron", "gene_short_name"))], ".g1", sep = "")
    results.g1 <- results
    sj.count.total <- apply(df.sj.count.g2, 1, function(x) {
      sum(x)
    })
    results.sj <- data.frame(coord.intron = names(sj.count.total), 
                             sj.count.total = sj.count.total, stringsAsFactors = FALSE)
    row.names(results.sj) <- NULL
    results.sj <- left_join(results.sj, sj.metadata[, c("coord.intron", 
                                                   "gene_short_name.start")], by = "coord.intron")
    gene.count.total <- apply(df.gene.count.g2, 1, function(x) {
      sum(x)
    })
    results.gene <- data.frame(gene_short_name.start = names(gene.count.total), 
                               gene.count.total = gene.count.total, stringsAsFactors = FALSE)
    results <- left_join(results.sj, results.gene, by = "gene_short_name.start")
    results$psi <- round(results$sj.count.total/results$gene.count.total * 
                           100, digits = 2)
    names(results)[which(names(results) == "gene_short_name.start")] <- "gene_short_name"
    names(results)[-which(names(results) %in% c("coord.intron", 
                                                "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% 
                                                                                                      c("coord.intron", "gene_short_name"))], ".g2", sep = "")
    results.g2 <- results
    results <- left_join(results.g1, results.g2, by = "coord.intron")
    results$delta.perm <- results$psi.g2 - results$psi.g1
    row.names(results) <- results$coord.intron
    results <- results[, "delta.perm", drop = FALSE]
    .list.results.perm[[i]] <- results
    if (show.progress == TRUE) {
      setTxtProgressBar(pb, i)
    }
  }
  results.perm <- do.call(cbind.data.frame, .list.results.perm)
  message("Computing P values...")
  index.l <- table(results.obs$coord.intron == row.names(results.perm))
  index.true <- length(which(names(index.l) == TRUE))
  index.false <- length(which(names(index.l) == FALSE))
  if (index.true == 1 & index.false == 0) {
    results.perm <- cbind.data.frame(results.perm, results.obs[, 
                                                               "delta", drop = FALSE])
  }
  else {
    return(message("Error in consolidating observed and permutated results"))
  }
  index.delta.perm <- grep("perm", names(results.perm))
  index.delta.obs <- which(names(results.perm) == "delta")
  pval <- apply(results.perm[, ], 1, function(x) {
    delta.obs <- x[index.delta.obs]
    delta.perm <- x[index.delta.perm]
    pval <- sum(abs(delta.perm) > abs(delta.obs))/n.iterations
    return(pval)
  })
  results.obs$pval <- pval
  results.obs <- results.obs[order(results.obs$pval), ]
  results.obs <- left_join(results.obs, mean.combined.df, by = "gene_short_name")
  MarvelObject$DE$SJ$Table <- results.obs
  MarvelObject$DE$SJ$cell.group.g1 <- cell.group.g1
  MarvelObject$DE$SJ$cell.group.g2 <- cell.group.g2
  return(MarvelObject)
}
