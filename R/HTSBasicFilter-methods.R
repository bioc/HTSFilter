#' Implement basic filters for transcriptome sequencing data.
#' 
#' Implement a variety of basic filters for transcriptome sequencing data.
#' 
#' This function implements a basic filter for high-throughput sequencing data for a variety of filter types: 
#' mean, sum, RPKM, variance, CPM, maximum, mean CPM values, the sum of CPM values, the variance of CPM 
#' values, maximum CPM value, mean RPKM values, the sum of RPKM values, the variance of RPKM values, or 
#' the maximum RPKM value. The filtering criteria used may be for a given cutoff value, a number of genes, 
#' or a given quantile value.
#' 
#' 
#' @export
#' 
#' @param x A numeric matrix or data.frame representing the counts of dimension (\emph{g} x \emph{n}), 
#' for \emph{g} genes in \emph{n} samples, a \code{DGEList} object, a 
#' \code{DGEExact} object, a \code{DGEGLM} object, a \code{DGELRT} object, or a \code{DESeqDataSet} object.
#' 
#' @param method  Basic filtering method to be used: \dQuote{mean}, \dQuote{sum}, \dQuote{rpkm}, 
#' \dQuote{variance}, \dQuote{cpm}, \dQuote{max}, \dQuote{cpm.mean}, \dQuote{cpm.sum}, \dQuote{cpm.variance}, 
#' \dQuote{cpm.max}, \dQuote{rpkm.mean}, \dQuote{rpkm.sum}, \dQuote{rpkm.variance}, or \dQuote{rpkm.max}
#' 
#' @param cutoff.type Type of cutoff to be used: a numeric value indicating the number of samples to be 
#' used for filtering (when \code{method} = \dQuote{cpm} or \dQuote{rpkm}), or one of \dQuote{value}, 
#' \dQuote{number}, or \dQuote{quantile}
#' 
#' @param cutoff Cutoff to be used for chosen filter
#' 
#' @param length Optional vector of length \emph{n} containing the lengths of each gene in \code{x}; 
#' optional except in the case of \code{method} = \dQuote{rpkm}
#' 
#' @param normalization Normalization method to be used to correct for differences in library sizes, with 
#' choices  \dQuote{TMM} (Trimmed Mean of M-values), \dQuote{DESeq} (normalization method proposed in the
#' DESeq package), \dQuote{pseudo.counts} (pseudo-counts obtained via quantile-quantile normalization in 
#' the edgeR package, only available for objects of class \code{DGEList} and \code{DGEExact}), and 
#' \dQuote{none} (to be used only if user is certain no normalization is required, or if data have already 
#' been pre-normalized by an alternative method)
#' 
#' @param pAdjustMethod The method used to adjust p-values, see \code{?p.adjust}
#' 
#' @param ... Additional optional arguments
#' 
#' @return 
#' \itemize{
#'  \item{filteredData }{An object of the same class as \code{x} containing the data that passed the filter}
#' 
#'  \item{on }{A binary vector of length \emph{g}, where 1 indicates a gene with normalized expression
#'   greater than the optimal filtering threshold \code{s.optimal} in at least one sample (irrespective of 
#'   condition labels), and 0 indicates a gene with normalized expression less than or equal to the optimal 
#'   filtering threshold in all samples}
#' 
#'  \item{normFactor }{A vector of length \emph{n} giving the estimated library sizes estimated by the
#'   normalization method specified in \code{normalization}}
#' 
#'  \item{removedData }{A matrix containing the filtered data}
#' 
#'  \item{filterCrit }{A vector or matrix containing the criteria used to perform filtering}
#'  }
#' 
#' @references 
#' R. Bourgon, R. Gentleman, and W. Huber. (2010) Independent filtering increases detection power for high-
#' throughput experiments. \emph{PNAS} \bold{107}(21):9546-9551.
#' 
#' A. Rau, M. Gallopin, G. Celeux, F. Jaffrezic (2013). Data-based filtering 
#' for replicated high-throughput transcriptome sequencing experiments. \emph{Bioinformatics},
#' doi: 10.1093/bioinformatics/btt350.
#' 
#' @author 
#' Andrea Rau, Melina Gallopin, Gilles Celeux, and Florence Jaffrezic
#' 
#' @example /inst/examples/HTSBasicFilter.R
#' 
#' @keywords methods
#' @importMethodsFrom DESeq2 counts
#' @importFrom stats loess predict quantile var
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot lines points abline legend
#' @importFrom edgeR DGEList calcNormFactors
setGeneric(
  name = "HTSBasicFilter",
  def = function(x, ...) 
    {standardGeneric("HTSBasicFilter")}
)


#' @rdname HTSBasicFilter
## matrix
setMethod(
  f= "HTSBasicFilter",
  signature = signature(x="matrix"),
  definition = function(x, method, cutoff.type="value", cutoff=10, 
                        length=NA, normalization=c("TMM", "DESeq", "none")) 
  {
    normalization <- match.arg(normalization)
    data <- x
    
    ## Run filter
    filter <- .HTSBasicFilterBackground(data=data, method=method,
                                        cutoff.type=cutoff.type, cutoff=cutoff, length=length,
                                        normalization=normalization)
    
    ## Return various results
    filter.results <- list(filteredData =  filter$filteredData,
                           on = filter$on, normFactor = filter$normFactor,
                           removedData = filter$removedData, filterCrit = filter$filterCrit)
    return(filter.results)
  }
)

#' @rdname HTSBasicFilter
## data.frame
setMethod(
  f= "HTSBasicFilter",
  signature = signature(x="data.frame"),
  definition = function(x, method, cutoff.type="value", cutoff=10, 
                        length=NA, normalization=c("TMM", "DESeq", "none")) 
    
  {
    normalization <- match.arg(normalization)
    data <- as.matrix(x)
    if(is.character(data) == TRUE) {
      stop("Character values detected in data.frame x.\nPlease check that ID names are not included in one of the columns.")
    }
    
    ## Run filter
    filter <- .HTSBasicFilterBackground(data=data, method=method,
                                        cutoff.type=cutoff.type, cutoff=cutoff, length=length,
                                        normalization=normalization)
    
    ## Return various results
    filter.results <- list(filteredData =  filter$filteredData,
                           on = filter$on, normFactor = filter$normFactor,
                           removedData = filter$removedData, filterCrit = filter$filterCrit)
    return(filter.results)
  }
)



#' @rdname HTSBasicFilter
## DGEList
setMethod(
  f="HTSBasicFilter",
  signature = signature(x="DGEList"),
  definition = function(x, method, cutoff.type="value", cutoff=10, 
                        length=NA, normalization=c("TMM", "DESeq", "pseudo.counts", "none")) 
  {
    
    if(!require(edgeR)) {
      stop("edgeR library must be installed.")
    }
    if(is.null(x$common.dispersion) == FALSE & is.null(x$tagwise.dispersion == TRUE)) {
      stop("Filtering must be performed either before calling estimateCommonDisp, or after calling estimateTagwiseDisp.")
    }
    
    normalization <- match.arg(normalization)
    data <- x$counts
    conds <- x$samples$group
    if(normalization == "pseudo.counts") {
      if(is.null(x$pseudo.counts) == TRUE) {
        stop(paste("To use pseudo.counts for filtering, you must first call estimateCommonDisp."))
      }
      data <- x$pseudo.counts
      normalization <- "none"
    }
    
    ## Run filter
    filter <- .HTSBasicFilterBackground(data=data, method=method,
                                        cutoff.type=cutoff.type, cutoff=cutoff, length=length,
                                        normalization=normalization)
    on <- filter$on
    on.index <- which(on == 1)
    filteredData <- x
    
    ## Create a new DGEList
    filteredData <- x
    filteredData$counts <- x$counts[on.index,]
    filteredData$pseudo.counts <- x$pseudo.counts[on.index,]
    ## logCPM removed April 3, 2015 because no longer included in DGEList
    ## filteredData$logCPM <- x$logCPM[on.index]
    filteredData$tagwise.dispersion <- x$tagwise.dispersion[on.index]
    ## Added August 8, 2013: thanks to Marie-Laure Endale for catching this one!
    filteredData$genes <- x$genes[on.index]
    filteredData$AveLogCPM <- x$AveLogCPM[on.index]
    filteredData$trended.dispersion <- x$trended.dispersion[on.index]
    filteredData$offset <- x$offset[on.index,]
    
    ## Reset library sizes if filtering before estimating dispersion parameters 
    if(is.null(x$common.dispersion) == TRUE) {
      filteredData$samples$lib.size = colSums(filteredData$counts)
    }
    
    ## Return various results
    nf <- filter$norm.factor
    if(normalization == "pseudo.counts") nf <- "pseudo.counts"
    
    filter.results <- list(filteredData =  filteredData,
                           on = filter$on, normFactor = nf,
                           removedData = filter$removedData, filterCrit = filter$filterCrit)
    
    return(filter.results)
  }
)


#' @rdname HTSBasicFilter
## DGEExact
setMethod(
  f="HTSBasicFilter",
  signature = signature(x="DGEExact"),
  definition = function(x, method, cutoff.type="value", cutoff=10, 
                        length=NA, normalization=c("TMM", "DESeq", "pseudo.counts", "none")) 
  {
    
    if(!require(edgeR)) {
      stop("edgeR library must be installed.")
    }
    normalization <- match.arg(normalization)
    data <- DGEList$counts
    conds <- DGEList$samples$group
    if(normalization == "pseudo.counts") {
      if(is.null(DGEList$pseudo.counts) == TRUE) {
        stop(paste("To use pseudo.counts for filtering, you must first call estimateCommonDisp."))
      }
      data <- DGEList$pseudo.counts
      normalization <- "none"
    }
    
    ## Run filter
    filter <- .HTSBasicFilterBackground(data=data, method=method,
                                        cutoff.type=cutoff.type, cutoff=cutoff, length=length,
                                        normalization=normalization)
    on <- filter$on
    on.index <- which(on == 1) 
    
    ## Create a new DGEExact
    filteredData <- x
    filteredData$table <- x$table[on.index,]
    filteredData$genes <- x$genes[on.index]
    
    ## Return various results
    nf <- filter$norm.factor
    if(normalization == "pseudo.counts") nf <- "pseudo.counts"
    filter.results <- list(filteredData =  filteredData,
                           on = filter$on, normFactor = nf,
                           removedData = filter$removedData, filterCrit = filter$filterCrit)
    return(filter.results)
  }
)

#' @rdname HTSBasicFilter
## DGEGLM
setMethod(
  f="HTSBasicFilter",
  signature = signature(x="DGEGLM"),
  definition = function(x, method, cutoff.type="value", cutoff=10, 
                        length=NA, normalization=c("TMM", "DESeq", "none"))
  {
    if(!require(edgeR)) {
      stop("edgeR library must be installed.")
    }
    normalization <- match.arg(normalization)
    data <- x$counts
    
    ## Run filter
    filter <- .HTSBasicFilterBackground(data=data, method=method,
                                        cutoff.type=cutoff.type, cutoff=cutoff, length=length,
                                        normalization=normalization)
    on <- filter$on
    on.index <- which(on == 1) 
    
    ## Create a new DGEGLM
    filteredData <- x
    filteredData$counts <- x$counts[on.index,]
    filteredData$coefficients <- x$coefficients[on.index,]
    filteredData$df.residual <- x$df.residual[on.index]
    filteredData$deviance <- x$deviance[on.index]
    filteredData$genes <- x$genes[on.index]
    filteredData$dispersion <- x$dispersion[on.index]
    filteredData$weights <- x$weights[on.index]
    filteredData$fitted.values <- x$fitted.values[on.index,]
    filteredData$abundance <- x$abundance[on.index]
    filteredData$offset <- x$offset[on.index,]
    ## Added August 8, 2013: thanks to Marie-Laure Endale for catching this one!
    filteredData$AveLogCPM <- x$AveLogCPM[on.index]
    ## Added April 3, 2015
    filteredData$unshrunk.coefficients <- x$unshrunk.coefficients[on.index,]
    
    ## Return various results
    nf <- filter$norm.factor
    filter.results <- list(filteredData =  filteredData,
                           on = filter$on, normFactor = nf,
                           removedData = filter$removedData, filterCrit = filter$filterCrit)
    return(filter.results)
  }
)

#' @rdname HTSBasicFilter
## DGELRT
setMethod(
  f="HTSBasicFilter",
  signature = signature(x="DGELRT"),
  definition = function(x, method, cutoff.type="value", cutoff=10, 
                        length=NA, normalization=c("TMM", "DESeq", "none"))
  {
    
    if(!require(edgeR)) {
      stop("edgeR library must be installed.")
    }
    normalization <- match.arg(normalization)
    data <- x$counts
    
    ## Run filter
    filter <- .HTSBasicFilterBackground(data=data, method=method,
                                        cutoff.type=cutoff.type, cutoff=cutoff, length=length,
                                        normalization=normalization)
    on <- filter$on
    on.index <- which(on == 1) 
    
    ## Create a new DGELRT
    filteredData <- x
    filteredData$table <- x$table[on.index,]
    filteredData$coefficients <- x$coefficients[on.index,]
    ## Added April 3, 2015:
    filteredData$unshrunk.coefficients <- x$unshrunk.coefficients[on.index,]
    filteredData$df.test <- x$df.test[on.index]
    ## Removed April 3, 2015:
    ## filteredData$genes <- x$genes[on.index]
    ## filteredData$weights <- x$weights[on.index,]
    ## filteredData$abundance <- x$abundance[on.index]
    filteredData$df.residual <- x$df.residual[on.index]
    filteredData$dispersion <- x$dispersion[on.index]
    filteredData$fitted.values <- x$fitted.values[on.index,]
    filteredData$deviance <- x$deviance[on.index]
    filteredData$offset <- x$offset[on.index,]
    ## Added August 8, 2013: thanks to Marie-Laure Endale for catching this one!
    filteredData$AveLogCPM <- x$AveLogCPM[on.index]
    
    ## Return various results
    nf <- filter$norm.factor
    filter.results <- list(filteredData =  filteredData,
                           on = filter$on, normFactor = nf,
                           removedData = filter$removedData, filterCrit = filter$filterCrit)
    return(filter.results)
  }
)


#' @rdname HTSBasicFilter
## DESeqDataSet
setMethod(
  f= "HTSBasicFilter",
  signature = signature(x="DESeqDataSet"),
  definition = function(x, method, cutoff.type="value", cutoff=10,
                        length=NA, normalization=c("DESeq", "TMM", "none"), pAdjustMethod="BH")
  {
    if(!require(DESeq2)) {
      stop("DESeq2 library must be installed.")
    }          
    normalization <- match.arg(normalization)
    data <- counts(x, normalized=FALSE)
    
    ## Run filter
    filter <- .HTSBasicFilterBackground(data=data, method=method,
                                        cutoff.type=cutoff.type, cutoff=cutoff, length=length,
                                        normalization=normalization)
    on <- filter$on
    on.index <- which(on == 1) 
    filteredData <- x[on.index,]
    
    ## April 3, 2015: remove this as p-values are not adjusted in DESeqDataSet			
    #		if(length(colnames(mcols(filteredData))) > 0) {
    #                ## Re-adjust p-values, if they exist
    #                nm <- strsplit(colnames(mcols(filteredData)), split="_", fixed=TRUE)
    #                Waldindex <- which(unlist(lapply(nm, function(yy) yy[1]))=="WaldPvalue")
    #                LRTindex <- which(unlist(lapply(nm,  function(yy) yy[1]))=="LRTPvalue")
    #                # Wald p-values
    #                if(length(Waldindex) > 0 ) {
    #                  for(j in Waldindex) {
    #                    look <- substr(colnames(mcols(filteredData))[j], 12, 100)
    #                    find <- which(substr(colnames(mcols(filteredData)), 12+3, 100) == look)
    #                    find <- find[which(find > j)]
    #                    mcols(filteredData)[,find] <- p.adjust(mcols(filteredData)[,j], method=pAdjustMethod)
    #                  }
    #                }
    #                # LRT p-values
    #                if(length(LRTindex) > 0 ) {
    #                  for(j in LRTindex) {
    #                    look <- substr(colnames(mcols(filteredData))[j], 11, 100)
    #                    find <- which(substr(colnames(mcols(filteredData)), 11+3, 100) == look)
    #                    find <- find[which(find > j)]
    #                    mcols(filteredData)[,find] <- p.adjust(mcols(filteredData)[,j], method=pAdjustMethod)
    #                  }
    #                }
    #		}
    
    ## Return various results
    filter.results <- list(filteredData = filteredData,
                           on = filter$on, normFactor = filter$norm.factor,
                           removedData = data[which(filter$on == 0),], filterCrit = filter$filterCrit)
    return(filter.results)
  }
)

