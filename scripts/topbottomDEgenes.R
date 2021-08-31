log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("dplyr"))
suppressMessages(library("plyr"))
suppressMessages(library("tibble"))
suppressMessages(library("scales"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ComplexHeatmap"))

## SNAKEMAKE I/O ##
norm_counts <- snakemake@input[["normalized_counts"]]
diffexp <- snakemake@input[["diffexp"]]

## SNAKEMAKE PARAMS ##
levels <- snakemake@params[["levels"]]
designmatrix <- snakemake@params[["designmatrix"]]

## CODE ##
# Get normalized counts
norm_counts <- read.table(norm_counts)

# Get differential expression results
diffexp <- read.table(diffexp, header = TRUE, row.names = 1)

# Get design matrix
designmatrix <- read.table(designmatrix, header = TRUE, row.names = 1)

# Keep only significantly expressed genes
diffexp <- diffexp %>% rownames_to_column %>% filter(padj < 0.05)


# Format design matrix condition
designmatrix <- designmatrix %>% 
  mutate(condition = gsub("^\\*", "", condition)) %>% select(condition)

# All condition's levels
all_levels <- designmatrix %>% select(condition) %>% unique %>% pull

# Assign a color to each level
levels_colors <- setNames(hue_pal()(length(all_levels)), all_levels)
levels_colors <- list(condition = levels_colors)

# Subset the sample names according to the specified levels
designmatrix <- designmatrix %>% filter(condition %in% levels)
samples <- rownames(designmatrix)

# Get top 25 and bottom 25 DE genes
top <- diffexp %>% arrange(desc(log2FoldChange)) %>% 
  top_n(n = 25, wt = log2FoldChange) %>% select(rowname) %>% pull

bottom <-diffexp %>% arrange(log2FoldChange) %>% 
  top_n(n = -25, wt = log2FoldChange) %>% select(rowname) %>% pull

# Subset the normalized count matrix and scale by row
norm_counts <- t(apply(norm_counts[c(top, bottom), samples], 1, scale))
colnames(norm_counts) <- samples

# Specify the colors and breaks of the plot
colors <- colorRampPalette(c("blue", "white", "red"))(100)
break_limit <- round_any(max(max(norm_counts), abs(min(norm_counts))), 0.5,
                         ceiling)
breaks <- seq(-break_limit, break_limit, length.out = 101)

# Heatmaps
pdf(snakemake@output[["pdf"]], width = 9, height = 7)
pheatmap(norm_counts, name = "Gene Expression", color = colors, breaks = breaks,
         annotation_col = designmatrix, annotation_colors = levels_colors)
dev.off()

png(snakemake@output[["png"]], width = 9, height = 7, units = "in", res = 600)
pheatmap(norm_counts, name = "Gene Expression", color = colors, breaks = breaks,
         annotation_col = designmatrix, annotation_colors = levels_colors)
dev.off()