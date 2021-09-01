log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("DESeq2"))
suppressMessages(library("dplyr"))
suppressMessages(library("scales"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ComplexHeatmap"))

## SNAKEMAKE I/O ##
vsd <- snakemake@input[["vst"]]

## SNAKEMAKE PARAMS ##
levels <- snakemake@params[["levels"]]

## CODE ##
# Get vst
vsd <- readRDS(vsd)

# colData
coldata <- as.data.frame(colData(vsd)) %>% select(-sizeFactor)

# All levels
all_levels <- levels_colors <- lapply(asplit(coldata, 2), 
                                      function(x) unname(unique(x)))

# Assign a color to each level of each variable
all_levels_colors <- hue_pal()(length(unlist(all_levels)))
i <- 1
for (n in names(levels_colors)) {
  x <- levels_colors[[n]]
  levels_colors[[n]] <- setNames(all_levels_colors[i:(i - 1 + length(x))], x)
  i <- length(x) + 1
}

# Subset the sample names according to the specified levels
coldata <- coldata %>% filter(condition %in% levels)
samples <- rownames(coldata)

levels_colors <- setNames(lapply(names(levels_colors), function(y) {
  levels_colors[[y]][names(levels_colors[[y]]) %in% unique(coldata[[y]])]
}), names(levels_colors))

# Euclidean sample distances
sampleDist <- as.matrix(dist(t(assay(vsd[, samples])), method = "euclidean"))

# Colors
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Heatmaps
pdf(snakemake@output[["pdf"]], width = 9, height = 7)
pheatmap(sampleDist, name = "Euclidean Distance", color = colors, 
         annotation_col = coldata, annotation_colors = levels_colors)
dev.off()

png(snakemake@output[["png"]], width = 9, height = 7, units = "in", res = 600)
pheatmap(sampleDist, name = "Euclidean Distance", color = colors,
         annotation_col = coldata, annotation_colors = levels_colors)
dev.off()