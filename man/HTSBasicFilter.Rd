\name{HTSBasicFilter}
\docType{methods}
\alias{HTSBasicFilter}
\alias{HTSBasicFilter-methods}
\alias{HTSBasicFilter,matrix-method}
\alias{HTSBasicFilter,data.frame-method}
\alias{HTSBasicFilter,DGEList-method}
\alias{HTSBasicFilter,DGEExact-method}
\alias{HTSBasicFilter,DGEGLM-method}
\alias{HTSBasicFilter,DGELRT-method}
\alias{HTSBasicFilter,CountDataSet-method}
\title{
Implement basic filters for transcriptome sequencing data.
}
\description{
This function implements a variety of basic filters for transcriptome sequencing data.
}

\usage{
\S4method{HTSBasicFilter}{matrix}(x, method, cutoff.type="value", cutoff=10, 
	length=NA, normalization=c("TMM", "DESeq", "none"))

\S4method{HTSBasicFilter}{data.frame}(x, method, cutoff.type="value", cutoff=10, 
	length=NA, normalization=c("TMM", "DESeq", "none"))

\S4method{HTSBasicFilter}{CountDataSet}(x, method, cutoff.type="value", cutoff=10, 
	length=NA, normalization=c("DESeq", "TMM", "none"))

\S4method{HTSBasicFilter}{DGEList}(x, method, cutoff.type="value", cutoff=10, 
	length=NA, normalization=c("TMM", "DESeq", "pseudo.counts", "none"))

\S4method{HTSBasicFilter}{DGEExact}(x, method, cutoff.type="value", cutoff=10, 
	length=NA, normalization=c("TMM", "DESeq", "pseudo.counts", "none"))

\S4method{HTSBasicFilter}{DGEGLM}(x, method, cutoff.type="value", cutoff=10, 
	length=NA, normalization=c("TMM", "DESeq", "none"))

\S4method{HTSBasicFilter}{DGELRT}(x, method, cutoff.type="value", cutoff=10, 
	length=NA, normalization=c("TMM", "DESeq", "none"))
}

\arguments{
  \item{x}{
A numeric matrix or data.frame representing the counts of dimension (\emph{g} x \emph{n}), for \emph{g} 
genes in \emph{n} samples, a \code{CountDataSet} object, a \code{DGEList} object, a \code{DGEExact} object, 
a \code{DGEGLM} object, or a \code{DGELRT} object.
}
  \item{method}{
Basic filtering method to be used: \dQuote{mean}, \dQuote{sum}, \dQuote{rpkm}, \dQuote{variance}, 
\dQuote{cpm}, \dQuote{max}, \dQuote{cpm.mean}, \dQuote{cpm.sum}, \dQuote{cpm.variance}, \dQuote{cpm.max}, 
\dQuote{rpkm.mean}, \dQuote{rpkm.sum}, \dQuote{rpkm.variance}, or \dQuote{rpkm.max}
}
  \item{cutoff.type}{
Type of cutoff to be used: a numeric value indicating the number of samples to be used for filtering (when 
\code{method} = \dQuote{cpm} or \dQuote{rpkm}), or one of \dQuote{value}, \dQuote{number}, or \dQuote{quantile}
}
  \item{cutoff}{
Cutoff to be used for chosen filter
}
  \item{length}{
Optional vector of length \emph{n} containing the lengths of each gene in \code{x}; optional except in the 
case of \code{method} = \dQuote{rpkm}
}
  \item{normalization}{
Normalization method to be used to correct for differences in library sizes, with choices
\dQuote{TMM} (Trimmed Mean of M-values), \dQuote{DESeq} (normalization method proposed in the
DESeq package), \dQuote{pseudo.counts} (pseudo-counts obtained via quantile-quantile normalization in the
edgeR package, only available for objects of class \code{DGEList} and \code{DGEExact}), and \dQuote{none} 
(to be used only if user is certain no normalization is required, or if data have already been pre-normalized 
by an alternative method)
}
}


\details{
This function implements a basic filter for high-throughput sequencing data for a variety of filter types:
mean, sum, RPKM, variance, CPM, maximum, mean CPM values, the sum of CPM values, the variance of CPM values,
maximum CPM value, mean RPKM values, the sum of RPKM values, the variance of RPKM values, or the 
maximum RPKM value. The filtering criteria used may be for a given cutoff value, a number of genes, or a given
quantile value.
}

\value{
\item{filteredData }{An object of the same class as \code{x} containing the data that passed the filter}
%
\item{on }{A binary vector of length \emph{g}, where 1 indicates a gene with normalized expression
greater than the optimal filtering threshold \code{s.optimal} in at least one sample (irrespective of 
condition labels), and 0 indicates a gene with normalized expression less than or equal to the optimal 
filtering threshold in all samples}
%
\item{normFactor }{A vector of length \emph{n} giving the estimated library sizes estimated by the
normalization method specified in \code{normalization}}
%
\item{removedData }{A matrix containing the filtered data}
%
\item{filterCrit }{A vector or matrix containing the criteria used to perform filtering}
}

\references{

R. Bourgon, R. Gentleman, and W. Huber. (2010) Independent filtering increases detection power for high-
throughput experiments. \emph{PNAS} \bold{107}(21):9546-9551.

A. Rau, M. Gallopin, G. Celeux, F. Jaffrezic (2013). Independent data-based filtering 
for replicated high-throughput transcriptome sequencing experiments (submitted).

}
\author{
Andrea Rau, Melina Gallopin, Gilles Celeux, and Florence Jaffrezic
}
\examples{

data("sultan")
conds <- pData(sultan)$cell.line

########################################################################
## Matrix or data.frame
########################################################################

## Filter genes with total (sum) normalized gene counts < 10
filter <- HTSBasicFilter(exprs(sultan), method="sum", cutoff.type="value", 
	cutoff = 10)


########################################################################
## CountDataSet
########################################################################

## Filter genes with mean normalized gene counts less than the 40% quantile
cds <- newCountDataSet(exprs(sultan), conds)
filter <- HTSBasicFilter(cds, method="mean", cutoff.type="quantile", 
	cutoff = 0.4)


########################################################################
## DGEExact
########################################################################

## Filter genes with CPM values less than 100 in more than 2 samples
dge <- DGEList(counts=exprs(sultan), group=conds)
dge <- calcNormFactors(dge)
filter <- HTSBasicFilter(dge, method="cpm", cutoff.type=2, cutoff=100)


}
\keyword{ methods }

