library(Biobase)
data("sultan")
conds <- pData(sultan)$cell.line

########################################################################
## Matrix or data.frame
########################################################################

filter <- HTSFilter(exprs(sultan), conds, s.len=25, plot=FALSE)

########################################################################
## DGEExact
########################################################################

library(edgeR)
dge <- DGEList(counts=exprs(sultan), group=conds)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
et <- exactTest(dge)
et <- HTSFilter(et, DGEList=dge, s.len=25, plot=FALSE)$filteredData
## topTags(et)


########################################################################
## DESeq2
########################################################################

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exprs(sultan),
                              colData = data.frame(cell.line = conds),
                              design = ~ cell.line)
## Not run:
##
## dds <- DESeq(dds)
## filter <- HTSFilter(dds, s.len=25, plot=FALSE)$filteredData
## class(filter)
## res <- results(filter, independentFiltering=FALSE)
