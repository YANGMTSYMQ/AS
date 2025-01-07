adhocGene.DE.PSI.10X_re<-function (MarvelObject) 
{
  MarvelObject <- MarvelObject
  df <- MarvelObject$adhocGene$Expression$PSI$Table
  groups <- factor(levels(df$group), levels = levels(df$group))
  pairs <- gtools::combinations(n = length(groups), r = 2, 
                                v = as.numeric(groups), repeats.allowed = FALSE)
  pairs <- as.data.frame(pairs)
  levels <- levels(groups)
  for (i in 1:length(levels)) {
    pairs[pairs == i] <- levels[i]
  }
  coord.introns <- levels(as.factor(df$coord.intron))
  .list <- list()
  for (j in 1:length(coord.introns)) {
    df.small <- df[which(df$coord.intron == coord.introns[j]), 
    ]
    if (nrow(df.small) >= 2) {
      pairs.small <- pairs[which(pairs$V2 %in% df.small$group), 
      ]
      pairs.small <- pairs.small[which(pairs.small$V1 %in% 
                                         df.small$group), ]
      delta <- NULL
      for (i in 1:nrow(pairs.small)) {
        g1 <- pairs.small[i, 1]
        g2 <- pairs.small[i, 2]
        mean.g1 <- df.small[which(df.small$group == 
                                    g1), "psi"]
        mean.g2 <- df.small[which(df.small$group == 
                                    g2), "psi"]
        delta[i] <- mean.g2 - mean.g1
      }
      results <- data.frame(coord.intron = coord.introns[j], 
                            group.pair = paste(pairs.small[, 2], " vs ", 
                                               pairs.small[, 1], sep = ""), delta = delta, 
                            stringsAsFactors = FALSE)
      .list[[j]] <- results
    }
  }
  results <- do.call(rbind.data.frame, .list)
  . <- unique(df[, c("coord.intron", "figure.column")])
  results <- left_join(results, ., by = "coord.intron")
  cols <- c("group.pair", "figure.column", "coord.intron", 
            "delta")
  results <- results[, cols]
  MarvelObject$adhocGene$DE$PSI$Data <- results
  return(MarvelObject)
}
