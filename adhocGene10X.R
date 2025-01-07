Gene10X<-function (MarvelObject, cell.group.list, gene_short_name, log2.transform = TRUE, 
          min.pct.cells = 10, downsample = FALSE, seed = 1) 
{
  MarvelObject <- MarvelObject
  df.gene.norm <- MarvelObject$gene.norm.matrix
  cell.group.list <- cell.group.list
  gene_short_name <- gene_short_name
  min.pct.cells <- min.pct.cells
  downsample <- downsample
  seed <- seed
  log2.transform <- log2.transform
  .list <- list()
  set.seed(seed)
  if (downsample == TRUE) {
    n.cells <- min(sapply(cell.group.list, length))
    message(paste("Downsampling to ", n.cells, " cells per group", 
                  sep = ""))
    for (i in 1:length(cell.group.list)) {
      cell.ids <- cell.group.list[[i]]
      cell.ids <- sample(cell.ids, size = n.cells, replace = FALSE)
      .list[[i]] <- cell.ids
    }
    names(.list) <- names(cell.group.list)
    cell.group.list <- .list
  }
  df.gene.norm <- df.gene.norm[gene_short_name, ]
  df.gene.norm <- data.frame(cell.id = names(df.gene.norm), 
                             exp.gene.norm = as.numeric(df.gene.norm), stringsAsFactors = FALSE)
  row.names(df.gene.norm) <- NULL
  if (log2.transform == TRUE) {
    df.gene.norm$exp.gene.norm <- log2(df.gene.norm$exp.gene.norm + 
                                         1)
  }
  .list <- list()
  for (i in 1:length(cell.group.list)) {
    . <- data.frame(cell.id = cell.group.list[[i]], group = names(cell.group.list)[i])
    .list[[i]] <- .
  }
  ref <- do.call(rbind.data.frame, .list)
  df.gene.norm <- left_join(df.gene.norm, ref, by = "cell.id")
  df.gene.norm <- df.gene.norm[!is.na(df.gene.norm$group), 
  ]
  df.gene.norm$group <- factor(df.gene.norm$group, levels = names(cell.group.list))
  df.gene.norm <- df.gene.norm[order(df.gene.norm$group), 
  ]
  groups <- levels(df.gene.norm$group)
  mean.expr <- NULL
  pct.cells.expr <- NULL
  for (i in 1:length(groups)) {
    group <- groups[i]
    df.gene.norm.small <- df.gene.norm[which(df.gene.norm$group == 
                                               group), ]
    mean.expr[i] <- mean(df.gene.norm.small$exp.gene.norm)
    pct.cells.expr[i] <- sum(df.gene.norm.small$exp.gene.norm != 
                               0)/nrow(df.gene.norm.small) * 100
  }
  results <- data.frame(group = groups, mean.expr = mean.expr, 
                        pct.cells.expr = pct.cells.expr, stringsAsFactors = FALSE)
  results$pct.cells.expr[which(results$pct.cells.expr < min.pct.cells)] <- NA
  results$group <- factor(results$group, levels = names(cell.group.list))
  results <- results[order(results$group), ]
  data <- results
  x <- as.character(rep(1, times = nrow(data)))
  y <- rep(c(nrow(data):1))
  z1 <- data$pct.cells.expr
  z2 <- data$mean.expr
  maintitle <- gene_short_name
  xtitle <- ""
  ytitle <- ""
  legendtitle.size <- "% cells"
  legendtitle.color <- "mean[log2(expr)]"
  labels.y <- data$group
  plot <- ggplot() + geom_point(data, mapping = aes(x = x, 
                                                    y = y, size = z1, color = z2)) + scale_color_gradientn(colors = c("gray", 
                                                                                                                      "cyan", "green", "yellow", "red")) + scale_y_continuous(breaks = y, 
                                                                                                                                                                              labels = labels.y) + labs(title = maintitle, x = xtitle, 
                                                                                                                                                                                                        y = ytitle, size = legendtitle.size, color = legendtitle.color) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_blank(), plot.title = element_text(size = 12, 
                                                                        hjust = 0.5), axis.line.x = element_blank(), 
          axis.line.y = element_line(colour = "black"), axis.ticks.x = element_blank(), 
          axis.title = element_text(size = 12), axis.text.x = element_blank(), 
          axis.text.y = element_text(size = 10, colour = "black"), 
          legend.title = element_text(size = 10), legend.text = element_text(size = 8), 
          legend.key = element_blank()) #+ scale_size_area(breaks = c(25, 
                                         #                            50, 75, 100), limits = c(1, 100))
  MarvelObject$adhocGene$Expression$Gene$Table <- results
  MarvelObject$adhocGene$Expression$Gene$Plot <- plot
  MarvelObject$adhocGene$cell.group.list <- cell.group.list
  MarvelObject$adhocGene$gene_short_name <- gene_short_name
  return(MarvelObject)
}
