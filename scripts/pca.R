log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("DESeq2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("patchwork"))
source("scripts/plotPCA.3.R")

## SNAKEMAKE I/O ##
vsd <- snakemake@input[["vst"]]

## SNAKEMAKE PARAMS ##
#condition <- snakemake@params[["condition"]]
#levels <- snakemake@params[["levels"]]

## CODE ##
# Get vst
vsd <- readRDS(vsd)

# PCA plots for the first 3 principal components
pca12 <- plotPCA.3(vsd)
pca13 <- plotPCA.3(vsd, pc = c(1, 3))
pca23 <- plotPCA.3(vsd, pc = c(2, 3))

# Make patchwork with the three plots
pca <- pca12 + pca13 + pca23 + 
  plot_layout(guides = "collect", widths = rep(1, 3)) & 
  theme(legend.position = "bottom")

# Save PCA plot
ggsave(filename = snakemake@output[["pdf"]], plot = pca, width = 20)
ggsave(filename = snakemake@output[["png"]], plot = pca, dpi = 600, width = 20)