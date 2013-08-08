\name{normalizeData}
\alias{normalizeData}
\title{
Normalize transcriptome sequencing data.
}
\description{
Normalize count-based measures of transcriptome sequencing data using the
Trimmed Means of M-values (TMM) or DESeq approach.
}

\usage{
normalizeData(data, normalization) 
}

\arguments{
  \item{data}{
A numeric matrix representing the counts of dimension (\emph{g} x \emph{n}), for \emph{g} genes in \emph{n} samples.
}
  \item{normalization}{
Normalization method to be used to correct for differences in library sizes, with choices
\dQuote{TMM} (Trimmed Mean of M-values), \dQuote{DESeq} (normalization method proposed in the
DESeq package), and \dQuote{none} 
}
}

\value{

\item{data.norm }{A numeric matrix representing the normalized counts of dimension (\emph{g} x \emph{n}), 
for \emph{g} genes in \emph{n} samples.}
%
\item{norm.factor }{A vector of length \emph{n} giving the estimated library sizes estimated by the
normalization method specified in \code{normalization}}
}

\references{
S. Anders and W. Huber (2010). Differential expression analysis for sequence count data.
\emph{Genome Biology}, 11(R106):1-28.

A. Rau, M. Gallopin, G. Celeux, F. Jaffrezic (2013). Data-based filtering 
for replicated high-throughput transcriptome sequencing experiments. \emph{Bioinformatics},
doi: 10.1093/bioinformatics/btt350.

M. D. Robinson and A. Oshlack (2010). A scaling normalization method for differential expression
analysis of RNA-seq data. \emph{Genome Biology}, 11(R25).
}
\author{
Andrea Rau, Melina Gallopin, Gilles Celeux, and Florence Jaffrezic
}
\examples{

data("sultan")
normData <- normalizeData(exprs(sultan), norm="DESeq") 

}
\keyword{ methods }

