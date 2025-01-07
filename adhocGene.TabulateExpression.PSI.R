PSI10X<-function (MarvelObject, min.pct.cells = 10) 
{
  MarvelObject <- MarvelObject
  sj.metadata <- MarvelObject$sj.metadata
  df.gene.count <- MarvelObject$gene.count.matrix
  df.sj.count <- MarvelObject$sj.count.matrix
  cell.group.list <- MarvelObject$adhocGene$cell.group.list
  gene_short_name <- MarvelObject$adhocGene$gene_short_name
  min.pct.cells <- min.pct.cells
  coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start == 
                                       gene_short_name), "coord.intron"]
  df.sj.count <- df.sj.count[coord.introns, ]
  .list <- list()
  for (i in 1:length(cell.group.list)) {
    df.sj.count.small <- df.sj.count[, cell.group.list[[i]]]
    sj.count.total <- apply(df.sj.count.small, 1, function(x) {
      sum(x)
    })
    n.cells.total <- ncol(df.sj.count.small)
    n.cells.expr.sj <- apply(df.sj.count.small, 1, function(x) {
      sum(x != 0)
    })
    pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 
                                 100, digits = 2)
    results <- data.frame(group = names(cell.group.list)[i], 
                          coord.intron = row.names(df.sj.count.small), n.cells.total = n.cells.total, 
                          n.cells.expr.sj = n.cells.expr.sj, pct.cells.expr.sj = pct.cells.expr.sj, 
                          sj.count.total = sj.count.total, stringsAsFactors = FALSE)
    row.names(results) <- NULL
    .list[[i]] <- results
  }
  results.sj.count <- do.call(rbind.data.frame, .list)
  df.gene.count <- df.gene.count[gene_short_name, , drop = FALSE]
  .list <- list()
  for (i in 1:length(cell.group.list)) {
    df.gene.count.small <- df.gene.count[, cell.group.list[[i]], 
                                         drop = FALSE]
    gene.count.total <- apply(df.gene.count.small, 1, function(x) {
      sum(x)
    })
    n.cells.total <- ncol(df.gene.count.small)
    n.cells.expr.gene <- apply(df.gene.count.small, 1, function(x) {
      sum(x != 0)
    })
    pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 
                                   100, digits = 2)
    results <- data.frame(group = names(cell.group.list)[i], 
                          gene_short_name = row.names(df.gene.count.small), 
                          n.cells.total = n.cells.total, n.cells.expr.gene = n.cells.expr.gene, 
                          pct.cells.expr.gene = pct.cells.expr.gene, gene.count.total = gene.count.total, 
                          stringsAsFactors = FALSE)
    row.names(results) <- NULL
    .list[[i]] <- results
  }
  results.gene.count <- do.call(rbind.data.frame, .list)
  results <- left_join(results.sj.count, results.gene.count[, c("group", 
                                                           "gene.count.total")], by = "group")
  results$psi <- round(results$sj.count.total/results$gene.count.total * 
                         100, digits = 2)
  results$psi[which(results$pct.cells.expr.sj < min.pct.cells)] <- NA
  results <- results[!is.na(results$psi), ]
  . <- aggregate(results$n.cells.expr.sj, list(results$coord.intron), 
                 function(x) {
                   sum(x)
                 })
  . <- .[order(.[, 2], decreasing = TRUE), ]
  results$coord.intron <- factor(results$coord.intron, levels = .[, 
                                                                  1])
  results$group <- factor(results$group, levels = names(cell.group.list))
  results <- results[order(results$group, results$coord.intron), 
  ]
  coord.introns <- levels(results$coord.intron)
  . <- data.frame(coord.intron = coord.introns, figure.column = paste("SJ-", 
                                                                      c(1:length(coord.introns)), sep = ""), stringsAsFactors = FALSE)
  results <- left_join(results, ., by = "coord.intron")
  cols.1 <- c("group", "figure.column", "coord.intron")
  cols.2 <- setdiff(names(results), cols.1)
  results <- results[, c(cols.1, cols.2)]
  data <- results
  x <- data$coord.intron
  y <- data$group
  z1 <- data$pct.cells.expr.sj
  z2 <- data$psi
  maintitle <- gene_short_name
  xtitle <- "location"
  ytitle <- "sample"
  legendtitle.size <- "% cells"
  legendtitle.color <- "PSI"
  labels.y <- data$group
  plot <- ggplot() + geom_point(data, mapping = aes(x = x, 
                                                    y = y, size = z1, color = z2)) + 
    scale_color_gradientn(colors = c("gray", "cyan", "green", "yellow", "red")) +
    scale_y_discrete(limits = rev(levels(y))) + 
    labs(title = maintitle, x = xtitle, y = ytitle, size = legendtitle.size, color = legendtitle.color) + 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                            plot.title = element_text(size = 12, hjust = 0.5), axis.line.x = element_blank(), 
                                          axis.line.y = element_line(colour = "black"), axis.ticks.x = element_blank(), 
                                           axis.title = element_text(size = 12), axis.text.x = element_blank(), 
                                            axis.text.y = element_text(size = 10, colour = "black"), 
                                          legend.title = element_text(size = 10), legend.text = element_text(size = 8), 
                                           legend.key = element_blank())# +
    #scale_size_area(breaks = c(0,100), limits = c(1, 100))
  plot
  MarvelObject$adhocGene$Expression$PSI$Table <- results
  MarvelObject$adhocGene$Expression$PSI$Plot <- plot
  return(MarvelObject)
}

