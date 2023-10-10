log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages(library("DESeq2"))
suppressMessages(library("openxlsx"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("AnnotationDbi"))

## PARALLELIZATION ##
parallel <- FALSE
if (snakemake@threads > 1) {
    suppressMessages(library("BiocParallel"))
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

## SNAKEMAKE I/O ##
dds <- snakemake@input[["dds"]]
geneID_to_geneSYMBOL <- snakemake@input[["geneID_to_geneSYMBOL"]]

## SNAKEMAKE PARAMS ##
condition <- snakemake@params[["condition"]]
levels <- snakemake@params[["levels"]]
specie <- snakemake@params[["specie"]]


# Call the proper library according to the specie.
if (specie == 'human'){
    database <- 'org.Hs.eg.db'
} else{
    database <- 'org.Mm.eg.db'
}

suppressMessages(library(database, character.only = TRUE))
object_db <- eval(parse(text=database))

## CODE ##
# Get dds
dds <- readRDS(dds)
geneID_to_geneSYMBOL <- read.table(geneID_to_geneSYMBOL, sep = "\t")

print(class(geneID_to_geneSYMBOL))
print(geneID_to_geneSYMBOL)

# Get results
res <- results(dds, contrast = c(condition, levels), alpha=0.05, parallel = parallel)

# Annotate the GeneSymbol and complete GeneName from the ENSEMBL Gene ID.
ensg_res <- rownames(res)
ensemblGene_DEA <- gsub("\\.[0-9]*$", "", ensg_res)

ensg_symbol <- as.data.frame(mapIds(org.Hs.eg.db, keys = ensemblGene_DEA,
                                    column = "SYMBOL", keytype = "ENSEMBL"))
colnames(ensg_symbol) <- "GeneSymbol"
ensg_genename <- as.data.frame(mapIds(org.Hs.eg.db, keys = ensemblGene_DEA,
                                      column = "GENENAME", keytype = "ENSEMBL"))
colnames(ensg_genename) <- "GeneName"
annot <- merge(ensg_symbol, ensg_genename, by = "row.names", all = TRUE)
rownames(res) <- ensemblGene_DEA
res <- merge(as.data.frame(res), annot, by.x = "row.names", by.y = 1, all = TRUE)
ENSG_res_list <- as.list(ensg_res)
names(ENSG_res_list) <- ensemblGene_DEA
rownames(res) <- ENSG_res_list[res$Row.names]
col_order <- c("GeneSymbol", "GeneName", "baseMean", "log2FoldChange",
               "lfcSE", "stat", "pvalue", "padj")
res <- res[, col_order]

# Sort by adjusted p-value
res <- res[order(res$padj, decreasing = FALSE), ]

# Save tsv
write.table(res, file = snakemake@output[["tsv"]], sep = "\t", quote = FALSE, 
            col.names = NA)

# Add Name column
res$EnsemblGeneID <- rownames(res)
res <- res[c("EnsemblGeneID", colnames(res)[1:(ncol(res)-1)])]

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
                      rule = "AND($E7>0, $I7<0.05, NOT(ISBLANK($I7)))", style = redStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res), rows = 7:(nrow(res)+6),
                      rule = "AND($E7>0, OR($I7>0.05, ISBLANK($I7)))", style = redPlainStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res), rows = 7:(nrow(res)+6),
                      rule = "AND($E7<0, $I7<0.05, NOT(ISBLANK($I7)))", style = greenStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res), rows = 7:(nrow(res)+6),
                      rule = "AND($E7<0, OR($I7>0.05, ISBLANK($I7)))", style = greenPlainStyle)
setColWidths(wb, sheet, 1:ncol(res), widths = 13)

# Save excel
saveWorkbook(wb, file = snakemake@output[["xlsx"]], overwrite = TRUE)



# GET SHRUNKEN LOG FOLD CHANGES.
coef <- paste0(c(condition, levels[1], "vs", levels[2]), collapse = "_")
res_shrink <- lfcShrink(dds, coef=coef, type="apeglm")

# Annotate the GeneSymbol and complete GeneName from the ENSEMBL Gene ID.
ensg_res <- rownames(res_shrink)
ensemblGene_DEA <- gsub("\\.[0-9]*$", "", ensg_res)

ensg_symbol <- as.data.frame(mapIds(org.Hs.eg.db, keys = ensemblGene_DEA,
                                    column = "SYMBOL", keytype = "ENSEMBL"))
colnames(ensg_symbol) <- "GeneSymbol"
ensg_genename <- as.data.frame(mapIds(org.Hs.eg.db, keys = ensemblGene_DEA,
                                      column = "GENENAME", keytype = "ENSEMBL"))
colnames(ensg_genename) <- "GeneName"
annot <- merge(ensg_symbol, ensg_genename, by = "row.names", all = TRUE)
rownames(res_shrink) <- ensemblGene_DEA
res_shrink <- merge(as.data.frame(res_shrink), annot, by.x = "row.names", by.y = 1, all = TRUE)
ENSG_res_list <- as.list(ensg_res)
names(ENSG_res_list) <- ensemblGene_DEA
rownames(res_shrink) <- ENSG_res_list[res_shrink$Row.names]
col_order <- c("GeneSymbol", "GeneName", "baseMean", "log2FoldChange",
               "lfcSE", "pvalue", "padj")
res_shrink <- res_shrink[, col_order]

# Sort by adjusted p-value
res_shrink <- res_shrink[order(res_shrink$padj, decreasing = FALSE), ]

# Save tsv
write.table(res_shrink, file = snakemake@output[["tsv_lfcShrink"]], sep = "\t", quote = FALSE, 
            col.names = NA)

# Add Name column
res_shrink$EnsemblGeneID <- rownames(res_shrink)
res_shrink <- res_shrink[c("EnsemblGeneID", colnames(res_shrink)[1:(ncol(res_shrink)-1)])]

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
                      rule = "AND($E7>0, $H7<0.05, NOT(ISBLANK($H7)))", style = redStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_shrink), rows = 7:(nrow(res_shrink)+6),
                      rule = "AND($E7>0, OR($H7>0.05, ISBLANK($H7)))", style = redPlainStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_shrink), rows = 7:(nrow(res_shrink)+6),
                      rule = "AND($E7<0, $H7<0.05, NOT(ISBLANK($H7)))", style = greenStyle)
conditionalFormatting(wb, sheet, cols = 1:ncol(res_shrink), rows = 7:(nrow(res_shrink)+6),
                      rule = "AND($E7<0, OR($H7>0.05, ISBLANK($H7)))", style = greenPlainStyle)
setColWidths(wb, sheet, 1:ncol(res_shrink), widths = 13)

# Save excel
saveWorkbook(wb, file = snakemake@output[["xlsx_lfcShrink"]], overwrite = TRUE)