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
res_sort <- res[order(res$padj, decreasing = FALSE), ]

# Save tsv
write.table(res_sort, file = snakemake@output[["tsv"]], sep = "\t", quote = FALSE, 
            col.names = NA)

# Add Name column
res_sort$Name <- rownames(res_sort)
res_sort <- res_sort[c("Name", colnames(res_sort)[1:(ncol(res_sort)-1)])]

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
writeData(wb, sheet, res_sort, startRow = 6)
addStyle(wb, sheet, cols = 1:ncol(res_sort), rows = 6, style = boldStyle, 
         gridExpand = TRUE)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_sort), rows = 7:(nrow(res_sort)+6),
                      rule = "AND($C7>0, $G7<0.05, NOT(ISBLANK($G7)))", style = redStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_sort), rows = 7:(nrow(res_sort)+6),
                      rule = "AND($C7>0, OR($G7>0.05, ISBLANK($G7)))", style = redPlainStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_sort), rows = 7:(nrow(res_sort)+6),
                      rule = "AND($C7<0, $G7<0.05, NOT(ISBLANK($G7)))", style = greenStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_sort), rows = 7:(nrow(res_sort)+6),
                      rule = "AND($C7<0, OR($G7>0.05, ISBLANK($G7)))", style = greenPlainStyle)
setColWidths(wb, sheet, 1:ncol(res_sort), widths = 13)

# Save excel
saveWorkbook(wb, file = snakemake@output[["xlsx"]], overwrite = TRUE)



# GET SHRUNKEN LOG FOLD CHANGES.
coef <- paste0(c(condition, levels[1], "vs", levels[2]), collapse = "_")
res_shrink <- lfcShrink(dds, coef=coef, res = res, type="apeglm")

# Sort by adjusted p-value
res_shrink <- res_shrink[order(res_shrink$padj, decreasing = FALSE), ]

# Save tsv
write.table(res_shrink, file = snakemake@output[["tsv_lfcShrink"]], sep = "\t", quote = FALSE, 
            col.names = NA)

# Add Name column
res_shrink$Name <- rownames(res_shrink)
res_shrink <- res_shrink[c("Name", colnames(res_shrink)[1:(ncol(res_shrink)-1)])]

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
writeData(wb, sheet, res_shrink, startRow = 6)
addStyle(wb, sheet, cols = 1:ncol(res_shrink), rows = 6, style = boldStyle, 
         gridExpand = TRUE)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_shrink), rows = 7:(nrow(res_shrink)+6),
                      rule = "AND($C7>0, $F7<0.05, NOT(ISBLANK($F7)))", style = redStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_shrink), rows = 7:(nrow(res_shrink)+6),
                      rule = "AND($C7>0, OR($F7>0.05, ISBLANK($F7)))", style = redPlainStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_shrink), rows = 7:(nrow(res_shrink)+6),
                      rule = "AND($C7<0, $F7<0.05, NOT(ISBLANK($F7)))", style = greenStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_shrink), rows = 7:(nrow(res_shrink)+6),
                      rule = "AND($C7<0, OR($F7>0.05, ISBLANK($F7)))", style = greenPlainStyle)
setColWidths(wb, sheet, 1:ncol(res_shrink), widths = 13)

# Save excel
saveWorkbook(wb, file = snakemake@output[["xlsx_lfcShrink"]], overwrite = TRUE)