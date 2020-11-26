#' Filter replicated high-throughput transcriptome sequencing data
#'
#' This package implements a filtering procedure for replicated
#' transcriptome sequencing data based on a global Jaccard similarity index in
#' order to identify genes with low, constant levels of expression across one or
#' more experimental conditions.
#'
#' \tabular{ll}{ Package: \tab HTSFilter\cr Type: \tab Package\cr Version:
#' \tab 1.31.1\cr Date: \tab 2020-11-26\cr License: \tab Artistic-2.0 \cr LazyLoad:
#' \tab yes\cr }
#'
#' @name HTSFilter-package
#' @aliases HTSFilter-package
#' @docType package
#' @author Andrea Rau, Melina Gallopin, Gilles Celeux, and Florence Jaffrezic
#'
#' Maintainer: Andrea Rau <\url{andrea.rau@inrae.fr}>
#' @references
#' R. Bourgon, R. Gentleman, and W. Huber. (2010) Independent filtering increases detection power for high-
#' throughput experiments. \emph{PNAS} \bold{107}(21):9546-9551.
#'
#' P. Jaccard (1901). Etude comparative de la distribution
#' orale dans une portion des Alpes et des Jura.
#' \emph{Bulletin de la Societe Vaudoise des Sciences Naturelles}, \bold{37}:547-549.
#'
#' A. Rau, M. Gallopin, G. Celeux, F. Jaffrezic (2013). Data-based filtering
#' for replicated high-throughput transcriptome sequencing experiments. \emph{Bioinformatics},
#' doi: 10.1093/bioinformatics/btt350.
#'
#' @keywords package
#' @example /inst/examples/HTSFilter-package.R
#' @import Biobase
NULL




#' RNA-seq data from humans in Sultan et al. (2008)
#'
#' This dataset represents RNA-seq data from humans in two conditions (Ramos B cell line and HEK293T), with two biological replicates per condition.
#' The ExpressionSet was downloaded from the ReCount online resource.
#'
#' @name sultan
#' @docType data
#' @references \url{data_blah.com}
#' @usage data(sultan)
#' @keywords datasets
#' @format An ExpressionSet named \code{sultan.eset} containing the phenotype data and
#' expression data for the Sultan et al. (2008) experiment. Phenotype data may be
#' accessed using the \code{pData} function, and expression data may be accessed
#' using the \code{exprs} function.
#' @return Object of class \sQuote{ExpressionSet}. Matrix of counts can be accessed after
#' loading the \sQuote{Biobase} package and calling \code{exprs(sultan))}.
#' @source ReCount online resource (http://bowtie-bio.sourceforge.net/recount).
#' @references
#' A. C. Frazee, B. Langmead, and J. T. Leek. ReCount: a multi-experiment resource
#' of analysis-ready RNA-seq gene count datasets. BMC Bioinformatics, 12(449), 2011.
#'
#' M. Sultan, M. H. Schulz, H. Richard, A. Magen, A. Klingenhoff, M. Scherf, M. Seifert,
#' T. Borodina, A. Soldatov, D. Parkhomchuk, D. Schmidt, S. O'Keefe, S. Haas,
#' M. Vingron, H. Lehrach, and M. L. Yaspo. A global view of gene activity and alternative
#' splicing by deep sequencing of the human transcriptome. Science, 15(5891):956-60, 2008.
NULL