log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("DESeq2"))
suppressMessages(library("openxlsx"))

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
res <- results(dds, contrast = c(condition, levels), alpha=0.05, parallel = parallel)

# Sort by adjusted p-value
res <- res[order(res$padj, decreasing = FALSE), ]

# Save tsv
write.table(res, file = snakemake@output[["tsv"]], sep = "\t", quote = FALSE, 
            col.names = NA)

# Add Name column
res$Name <- rownames(res)
res <- res[c("Name", colnames(res)[1:(ncol(res)-1)])]

# Green and red styles for formatting excel
redStyle <- createStyle(fontColour = "#FF1F00", bgFill = "#F6F600")
redPlainStyle <- createStyle(fontColour = "#FF1F00")
greenStyle <- createStyle(fontColour = "#008000", bgFill = "#F6F600")
greenPlainStyle <- createStyle(fontColour = "#008000")
boldStyle <- createStyle(textDecoration = c("BOLD"))

# Excel file
wb <- createWorkbook()
sheet <- "Sheet1"
addWorksheet(wb, sheet)

# Legend 
legend <- t(data.frame(paste("Upregulated in", levels[1]),
                       paste("Downregulated in", levels[1]), "FDR=0.05"))
writeData(wb, sheet, legend[, 1])
addStyle(wb, sheet, cols = 1, rows = 1, 
         style = createStyle(fontColour = "#FF1F00", fgFill = "#F6F600"))
addStyle(wb, sheet, cols = 1, rows = 2, 
         style = createStyle(fontColour = "#008000", fgFill = "#F6F600"))
addStyle(wb, sheet, cols = 1, rows = 3, style = boldStyle)
invisible(sapply(1:3, function(i) mergeCells(wb, sheet, cols = 1:3, rows = i)))

# Reorder genes according to adjusted p-value
writeData(wb, sheet, res, startRow = 6)
addStyle(wb, sheet, cols = 1:ncol(res), rows = 6, style = boldStyle, 
         gridExpand = TRUE)
conditionalFormatting(wb, sheet, cols = 1:ncol(res), rows = 7:(nrow(res)+6),
                      rule = "AND($C7>0, $G7<0.05)", style = redStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res), rows = 7:(nrow(res)+6),
                      rule = "AND($C7>0, $G7>0.05)", style = redPlainStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res), rows = 7:(nrow(res)+6),
                      rule = "AND($C7<0, $G7<0.05)", style = greenStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res), rows = 7:(nrow(res)+6),
                      rule = "AND($C7<0, $G7>0.05)", style = greenPlainStyle)
setColWidths(wb, sheet, 1:ncol(res), widths = 13)

# Save excel
saveWorkbook(wb, file = snakemake@output[["xlsx"]], overwrite = TRUE)
