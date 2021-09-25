log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("DESeq2"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("scales"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggrepel"))
suppressMessages(library("patchwork"))
source("scripts/plotPCA.3.R")
source("scripts/levelFunctions.R")

## SNAKEMAKE I/O ##
vsd <- snakemake@input[["vst"]]

## SNAKEMAKE PARAMS ##
levels <- snakemake@params[["levels"]]
design <- snakemake@params[["design"]]
ref <- snakemake@params[["ref_levels"]]

## CODE ##
# Get vst
vsd <- readRDS(vsd)

# Get variables to plot by from design
variables <- unlist(str_split(design, pattern = "\\*|:|\\+"))
variables <- rev(str_trim(str_remove(variables, "~")))

# colData
coldata <- as.data.frame(colData(vsd)) %>% select(-sizeFactor)

# Set names to reference levels for each column in coldata
ref <- setNames(ref, colnames(coldata))

# Convert all coldata columns to factors and relevel
coldata <- coldata %>% 
  mutate(across(everything(), factor_relevel, reference = ref))

# Get variables' all combined levels
all_levels <- expand.grid(rev(lapply(coldata, levels))) %>% 
  unite("combined", all_of(variables), sep = ":") %>% select(combined) %>% pull

# Assign a color to each level
colors <- setNames(hue_pal()(length(all_levels)), all_levels)

# Subset the sample names according to the specified levels
coldata <- coldata %>% filter(get(variables[1]) %in% levels)
samples <- rownames(coldata)

# Keep only specified levels
levels <- coldata %>% unite("combined", all_of(variables), sep = ":") %>% 
  select(combined) %>% pull %>% unique

# Subset the colors according to the specified levels
colors <- colors[levels]

# Subset samples
vsd <- vsd[, samples]

# PCA plots for the first 3 principal components
pca12 <- plotPCA.3(vsd, intgroup = variables) + 
  scale_color_manual(values = colors, breaks = all_levels)
pca13 <- plotPCA.3(vsd, intgroup = variables, pc = c(1, 3)) + 
  scale_color_manual(values = colors, breaks = all_levels)
pca23 <- plotPCA.3(vsd, intgroup = variables, pc = c(2, 3)) + 
  scale_color_manual(values = colors, breaks = all_levels)

# Make patchwork with the three plots
pca <- pca12 + pca13 + pca23 + 
  plot_layout(guides = "collect", widths = rep(1, 3)) & 
  theme(legend.position = "bottom")

# Save PCA plot
ggsave(filename = snakemake@output[["pdf"]], plot = pca, width = 20)
ggsave(filename = snakemake@output[["png"]], plot = pca, dpi = 600, width = 20)