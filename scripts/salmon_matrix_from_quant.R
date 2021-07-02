library('DESeq2')
library('readr')
library('tximeta')


## SNAKEMAKE I/O ##
metadata_cache      <- snakemake@output[['metadata_cache']]
tximeta_transcripts <- snakemake@output[['transcript_estimates']]
gene_level_matrix   <- snakemake@output[['gene_level_matrix']]

# NOTE: We don't really fetch the expanded list of quant files from the
# snakemake object since we don't really need it. It's there so that
# this rule does not get executed until all samples have been quantified.

## SNAKEMAKE PARAMS ##
salmon_quant_directory <- snakemake@params[['salmon_quant_directory']]
samples_files          <- snakemake@params[['samples']]

## SNAKEMAKE LOG ##
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


# SET metadata cache folder
setTximetaBFC(metadata_cache)

samples           <- read.csv(samples_files, sep="\t", header=TRUE)
rownames(samples) <- samples$sample
quant_files   <- file.path(salmon_quant_directory, samples$sample, "quant.sf")

# tximeta looks for two columns: names for samples and files for paths
samples$files <- quant_files
samples$names <- samples$sample

se  <- tximeta(coldata = samples, type = 'salmon')
gse <- summarizeToGene(se)

## This may seem dumb but it let us use the DESeqDataSet wrapper
## to transform gene-level estimates to integer counts taking into 
## account the average transcript lengths from tximeta
## An alternative is to just round them up, as with tximport.

dds        <- DESeqDataSet(se=gse, design=~1) #blind design
raw_counts <- counts(dds, normalized=FALSE)


saveRDS(se, tximeta_transcripts)
write.table(raw_counts,file = gene_level_matrix, sep = "\t", 
            quote = FALSE, col.names = NA, row.names = TRUE)
