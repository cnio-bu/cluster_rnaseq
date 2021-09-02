log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("DESeq2"))
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("limma"))

## PARALLELIZATION ##
parallel <- FALSE
if (snakemake@threads > 1) {
    suppressMessages(library("BiocParallel"))
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

## SNAKEMAKE I/O ##
counts <- snakemake@input[["counts"]]

## SNAKEMAKE PARAMS ##
samples <- snakemake@params[["samples"]]
designmatrix <- snakemake@params[["designmatrix"]]
ref <- snakemake@params[["ref_levels"]]
design <- snakemake@params[["design"]]
batch_variables <- snakemake@params[["batch"]]

## FUNCTION ##
factor_relevel <- function(x, reference) {
  relev <- relevel(as.factor(gsub("^\\*", "", as.character(x))),
                   ref = reference[cur_column()])
  return(relev)
}

## CODE ##
# Get counts
counts <- read.table(counts, header = TRUE, row.names = 1)

# Get samples for DEA
samples <- read.table(samples, header = TRUE, row.names = 1)
DEAsamples <- rownames(samples)[samples$diffexp == "True"]

# Get design matrix
designmatrix <- read.table(designmatrix, header = TRUE, row.names = 1)

# Subset the count matrix keeping only samples to be tested
counts <- counts %>% select(all_of(DEAsamples))

# Subset the design matrix keeping only samples to be tested
designmatrix <- designmatrix[colnames(counts), , drop = FALSE]

# Set names to reference levels for each column in design matrix
ref <- setNames(ref, colnames(designmatrix))

# Remove '*' prefix from design matrix cells, convert all columns to factors
# and relevel using ref
designmatrix <- designmatrix %>% 
  mutate(across(everything(), factor_relevel, reference = ref))

# DESeq2 from htseqCount output
dds <- DESeqDataSetFromMatrix(counts, colData = designmatrix, 
                              design = as.formula(design))

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# DEA
dds <- DESeq(dds, parallel = parallel)

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Data transformation
vsd <- vst(dds, blind = FALSE)

# Save objects
saveRDS(dds, file = snakemake@output[["dds"]])
write.table(norm_counts, file = snakemake@output[["normalized_counts"]], 
            sep = "\t", quote = FALSE, col.names = NA)
saveRDS(vsd, file = snakemake@output[["vst"]][1])

# If there is a batch, perform batch correction
if (!is.null(batch_variables)) {
  ### colData
  coldata <- as.data.frame(colData(vsd))
  
  # Batch correct the vst
  batch <- coldata %>% unite("combined", all_of(batch_variables), sep = ":") %>%
    select(combined) %>% pull
  assay(vsd) <- removeBatchEffect(assay(vsd), batch = batch)
  saveRDS(vsd, file = snakemake@output[["vst"]][2])
}