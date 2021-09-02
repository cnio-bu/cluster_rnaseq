log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("DESeq2"))
suppressMessages(library("dplyr"))
suppressMessages(library("scales"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ComplexHeatmap"))
source("scripts/levelFunctions.R")

## SNAKEMAKE I/O ##
vsd <- snakemake@input[["vst"]]

## SNAKEMAKE PARAMS ##
levels <- snakemake@params[["levels"]]
ref <- snakemake@params[["ref_levels"]]

## CODE ##
# Get vst
vsd <- readRDS(vsd)

# colData
coldata <- as.data.frame(colData(vsd)) %>% select(-sizeFactor)

# Set names to reference levels for each column in coldata
ref <- setNames(ref, colnames(coldata))

# Convert all coldata columns to factors and relevel
coldata <- coldata %>% 
  mutate(across(everything(), factor_relevel, reference = ref))

# Assign a color to each level of each variable
levels_colors <- color_levels(coldata)

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