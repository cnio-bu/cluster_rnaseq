# Modified DESeq2 plotPCA to plot PC3 and sample labels
plotPCA.3 <- function(object, intgroup = "condition", ntop = 500, 
                      returnData = FALSE, pc = c(1, 2)) {
  # Calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # Select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  # Perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select, ]))
  
  # The contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  
  # Add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # Assembly the data for the plot (modified)
  pc <- setNames(pc, paste0("PC", pc))
  d <- data.frame(PCx = pca$x[,pc[1]], PCy = pca$x[,pc[2]], group = group, 
                  intgroup.df, name = colnames(object))
  colnames(d)[1:2] <- names(pc)
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pc]
    return(d)
  }
  xlab <- paste0(names(pc)[1], ": ", round(percentVar[pc[1]] * 100), 
                 "% variance")
  ylab <- paste0(names(pc)[2], ": ", round(percentVar[pc[2]] * 100), 
                 "% variance")
  p <- ggplot(data = d, aes_string(x = names(pc)[1], y = names(pc)[2], 
                                   color = "group", label = "name")) + 
    geom_point(size = 3) + xlab(xlab) + ylab(ylab) +
    geom_text_repel(size = 3, force = 1.5) + 
    guides(color = guide_legend(title = intgroup)) + theme_minimal()
  return(p)
}