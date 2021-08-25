log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("DESeq2"))
suppressMessages(library("scales"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("patchwork"))
source("scripts/plotPCA.3.R")

## SNAKEMAKE I/O ##
vsd <- snakemake@input[["vst"]]

## SNAKEMAKE PARAMS ##
levels <- snakemake@params[["levels"]]

## CODE ##
# Get vst
vsd <- readRDS(vsd)

# Assign a color to each level
all_levels <- levels(colData(vsd)$condition)
colors <- setNames(hue_pal()(length(all_levels)), all_levels)

# Subset the colors according to the specified levels
colors <- colors[levels]

# Subset the sample names according to the specified levels
samples <- rownames(colData(vsd))[colData(vsd)$condition %in% levels]

# Subset samples
vsd <- vsd[, samples]

# Relevel levels to keep original order

# PCA plots for the first 3 principal components
pca12 <- plotPCA.3(vsd) + 
  scale_color_manual(values = colors, breaks = all_levels)
pca13 <- plotPCA.3(vsd, pc = c(1, 3)) + 
  scale_color_manual(values = colors, breaks = all_levels)
pca23 <- plotPCA.3(vsd, pc = c(2, 3)) + 
  scale_color_manual(values = colors, breaks = all_levels)

# Make patchwork with the three plots
pca <- pca12 + pca13 + pca23 + 
  plot_layout(guides = "collect", widths = rep(1, 3)) & 
  theme(legend.position = "bottom")

# Save PCA plot
ggsave(filename = snakemake@output[["pdf"]], plot = pca, width = 20)
ggsave(filename = snakemake@output[["png"]], plot = pca, dpi = 600, width = 20)