log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("plyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("tibble"))
suppressMessages(library("scales"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ComplexHeatmap"))
source("scripts/levelFunctions.R")

## SNAKEMAKE I/O ##
norm_counts <- snakemake@input[["normalized_counts"]]
diffexp <- snakemake@input[["diffexp"]]

## SNAKEMAKE PARAMS ##
levels <- snakemake@params[["levels"]]
designmatrix <- snakemake@params[["designmatrix"]]
ref <- snakemake@params[["ref_levels"]]
  
## CODE ##
# Get normalized counts
norm_counts <- read.table(norm_counts)

# Get differential expression results
diffexp <- read.table(diffexp, , header = TRUE, row.names = 1, sep = '\t', 
                      blank.lines.skip = FALSE, quote = "")
diffexp <- diffexp[c("baseMean", "log2FoldChange", "lfcSE", 
                   "stat", "pvalue", "padj")]
# Get design matrix
designmatrix <- read.table(designmatrix, header = TRUE, row.names = 1)

# Keep only significantly expressed genes
diffexp <- diffexp %>% rownames_to_column %>% filter(padj < 0.05)

# Subset the design matrix keeping only samples to be tested
designmatrix <- designmatrix[colnames(norm_counts), , drop = FALSE]

# Set names to reference levels for each column in design matrix
ref <- setNames(ref, colnames(designmatrix))

# Remove '*' prefix from design matrix cells, convert all columns to factors
# and relevel
designmatrix <- designmatrix %>% 
  mutate(across(everything(), factor_relevel, reference = ref))

# Assign a color to each level of each variable
levels_colors <- color_levels(designmatrix)

# Subset the sample names according to the specified levels
designmatrix <- designmatrix %>% filter(condition %in% levels) %>% 
  select(condition)
samples <- rownames(designmatrix)

levels_colors <- list(condition = levels_colors[["condition"]]
                      [names(levels_colors[["condition"]]) %in% levels])

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