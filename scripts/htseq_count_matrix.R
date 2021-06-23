log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

## SNAKEMAKE I/O ##
htseqcount <- snakemake@input[["quant"]]

## CODE ##
# Get counts
counts <- lapply(htseqcount, read.table, header = FALSE, row.names = 1)

# Merge counts
counts <- do.call("cbind", counts)

# Remove last rows
counts <- counts[!startsWith(rownames(counts), "__"), ]

# Add colnames
samples <- sapply(htseqcount, function(x) {
  sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(x))
})
colnames(counts) <- samples

# Save object
write.table(counts, file = snakemake@output[["counts"]], sep = "\t", 
            quote = FALSE, col.names = NA, row.names = TRUE)