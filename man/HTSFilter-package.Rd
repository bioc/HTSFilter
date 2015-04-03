\name{HTSFilter-package}
\alias{HTSFilter-package}
\docType{package}
\title{
Filter replicated high-throughput transcriptome sequencing data
}
\description{
This package implements a filtering procedure for replicated
transcriptome sequencing data based on a global Jaccard similarity index in 
order to identify genes with low, constant levels of expression across one or 
more experimental conditions.
}
\details{
\tabular{ll}{
Package: \tab HTSFilter\cr
Type: \tab Package\cr
Version: \tab 1.7.1\cr
Date: \tab 2015-03-04\cr
License: \tab Artistic-2.0\cr
}
}
\author{
Andrea Rau, Melina Gallopin, Gilles Celeux, and Florence Jaffrezic

Maintainer: Andrea Rau <andrea.rau@jouy.inra.fr>
}
\references{

R. Bourgon, R. Gentleman, and W. Huber. (2010) Independent filtering increases detection power for high-
throughput experiments. \emph{PNAS} \bold{107}(21):9546-9551.

P. Jaccard (1901). Etude comparative de la distribution orale dans une portion des Alpes et des Jura.
\emph{Bulletin de la Societe Vaudoise des Sciences Naturelles}, \bold{37}:547-549.

A. Rau, M. Gallopin, G. Celeux, F. Jaffrezic (2013). Data-based filtering 
for replicated high-throughput transcriptome sequencing experiments. \emph{Bioinformatics}, 
doi: 10.1093/bioinformatics/btt350.

}
\keyword{ package }
\examples{

data("sultan")
conds <- pData(sultan)$cell.line

########################################################################
## Matrix or data.frame
########################################################################

filter <- HTSFilter(exprs(sultan), conds, s.len=25, plot=FALSE)


########################################################################
## CountDataSet
########################################################################

library(DESeq)
cds <- newCountDataSet(exprs(sultan), conds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- HTSFilter(cds, s.len=25, plot=FALSE)$filteredData
class(cds)
## res <- nbinomTest(cds, levels(conds)[1], levels(conds)[2])


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



}
