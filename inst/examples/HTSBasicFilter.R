library(Biobase)
data("sultan")
conds <- pData(sultan)$cell.line
 
########################################################################
## Matrix or data.frame
########################################################################

## Filter genes with total (sum) normalized gene counts < 10
filter <- HTSBasicFilter(exprs(sultan), method="sum", cutoff.type="value", 
                        cutoff = 10)
                        
                        
########################################################################
## DGEExact
########################################################################

library(edgeR)
## Filter genes with CPM values less than 100 in more than 2 samples
dge <- DGEList(counts=exprs(sultan), group=conds)
dge <- calcNormFactors(dge)
filter <- HTSBasicFilter(dge, method="cpm", cutoff.type=2, cutoff=100)

########################################################################
## DESeq2
########################################################################

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = exprs(sultan),
                             colData = data.frame(cell.line = conds),
                               design = ~ cell.line)
                             
                             
## Not run: Filter genes with mean normalized gene counts < 40% quantile
## dds <- DESeq(dds)
## filter <- HTSBasicFilter(dds, method="mean", cutoff.type="quantile", 
##	cutoff = 0.4)
## res <- results(filter, independentFiltering=FALSE)
