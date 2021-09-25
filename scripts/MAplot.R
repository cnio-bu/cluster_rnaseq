log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("DESeq2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("patchwork"))

## PARALLELIZATION ##
parallel <- FALSE
if (snakemake@threads > 1) {
  suppressMessages(library("BiocParallel"))
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}

## SNAKEMAKE I/O ##
dds <- snakemake@input[["dds"]]

## SNAKEMAKE PARAMS ##
condition <- snakemake@params[["condition"]]
levels <- snakemake@params[["levels"]]

## CODE ##
# Get dds
dds <- readRDS(dds)

# Get results
res <- results(dds, contrast = c(condition, levels), alpha=0.05, 
               parallel = parallel)

# Log Fold Change shrinkage
coef <- paste0(c(condition, levels[1], "vs", levels[2]), collapse = "_")
res_apeglm <- lfcShrink(dds, coef = coef, type = "apeglm")
res_normal <- lfcShrink(dds, coef = coef, type = "normal")
res_ashr <- lfcShrink(dds, coef = coef, type = "ashr")

# MA plots
MA_apeglm <- ggmaplot(res_apeglm, fdr = 0.05, fc = 1.5, size = 0.7, top = 10,
                      main = "apeglm", legend = "top", font.main = "bold", 
                      font.legend = "bold", font.label = c("bold", 11), 
                      label.rectangle = TRUE, ggtheme = theme_minimal()) + 
  theme(plot.title = element_text(hjust = 0.5))
MA_normal <- ggmaplot(res_normal, fdr = 0.05, fc = 1.5, size = 0.7, top = 10,
                      main = "normal", legend = "top", font.main = "bold", 
                      font.legend = "bold", font.label = c("bold", 11), 
                      label.rectangle = TRUE, ggtheme = theme_minimal()) + 
  theme(plot.title = element_text(hjust = 0.5))
MA_ashr <- ggmaplot(res_ashr, fdr = 0.05, fc = 1.5, size = 0.7, top = 10,
                    main = "ashr", legend = "top", font.main = "bold", 
                    font.legend = "bold", font.label = c("bold", 11), 
                    label.rectangle = TRUE, ggtheme = theme_minimal()) + 
  theme(plot.title = element_text(hjust = 0.5))

# Make patchwork with the three plots
MA <- MA_apeglm + MA_normal + MA_ashr +
  plot_layout(widths = rep(5, 3))

# Save PCA plot
ggsave(filename = snakemake@output[["pdf"]], plot = MA, width = 20)
ggsave(filename = snakemake@output[["png"]], plot = MA, dpi = 600, width = 20)