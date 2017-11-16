#' Calculate data-based filtering threshold for replicated transcriptome sequencing data.
#' 
#' Calculate a data-based filtering threshold for replicated transcriptome
#' sequencing data through the pairwise Jaccard similarity index between pairs
#' of replicates within each experimental condition.
#' 
#' The Jaccard similarity index, which measures the overlap of two sets, is calculated as follows. 
#' Given two binary vectors, each of length \emph{n}, we define the following values:
#'   \itemize{
#'       \item \emph{a} = the number of attributes with a value of 1 in both vectors
#'       \item \emph{b} = the number of attributes with a value of 1 in the first vector and 0 in the second
#'       \item \emph{c} = the number of attributes with a value of 0 in the first vector and 1 in the second
#'       \item \emph{d} = the number of attributes with a value of 0 in both vectors
#'         }
#' We note that all attributes fall into one of these four quantities, so \eqn{a+b+c+d=n}. Given these
#' quantities, we may calculate the Jaccard similarity index between the two vectors as follows:
#'   \deqn{J = \frac{a}{a+b+c}.}{J = a/(a+b+c).}
#' 
#' @aliases 
#' HTSFilter
#' HTSFilter-methods
#' HTSFilter,matrix-method
#' HTSFilter,data.frame-method
#' HTSFilter,DGEList-method
#' HTSFilter,DGEExact-method
#' HTSFilter,DGEGLM-method
#' HTSFilter,DGELRT-method
#' HTSFilter,DESeqDataSet-method
#' 
#' @export
#' 
#' @param x A numeric matrix or data.frame representing the counts of dimension (\emph{g} x \emph{n}), 
#' for \emph{g} genes in \emph{n} samples, a \code{DGEList} object, a 
#' \code{DGEExact} object, a \code{DGEGLM} object, a \code{DGELRT} object, or a \code{DESeqDataSet} object.
#' 
#' @param conds  Vector of length \emph{n} identifying the experimental condition of each of the \emph{n} samples; required when sQuote(x)
#' is a numeric matrix. In the case of objects of class \code{DGEList}, 
#' \code{DGEExact}, \code{DGEGLM}, \code{DGELRT}, or \code{DESeqDataSet}, the design matrix is automatically
#’ detected from the object; if an alternative design should be used, this may optionally be provided.
#’ Note that in multi-factor designs in which one or more combinations of factors are unreplicated, the bonds
#’ vector may be used to specify only the primary factor of interest for use in HTSFilter.
#' 
#' @param s.min Minimum value of filtering threshold to be considered, with default value equal to 1
#' 
#' @param s.max Maximum value of filtering threshold to be considered, with default value equal to 200
#' 
#' @param s.len Length of sequence of filtering thresholds to be considered (from \code{s.min} to \code{s.max}) 
#'   for the calculation of the global similarity index
#' 
#' @param loess.span Span of the loess curve to be fitted to the filtering thresholds and corresponding global similarity
#'   indices, with default value equal to 0.3
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
#' @param plot If \dQuote{TRUE}, produce a plot of the calculated global similarity indices against the
#' filtering threshold with superimposed loess curve
#' 
#' @param plot.name  If \code{plot} = \dQuote{TRUE}, the name of the PDF file to be saved to the current working directory.
#'   If \code{plot.name} = NA, the plot is drawn in the current window.
#'   
#' @param DGEList Object of class DGEList, to be used when filtering objects of class DGEExact
#' 
#' @param DGEGLM Object of class DGEGLM, to be used when filtering objects of class DGELRT
#' 
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel 
#' execution using BiocParallel (see next argument \code{BPPARAM}). A note on running 
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects 
#' from the current R environment before calling the function, as it is possible that R's 
#' internal garbage collection will copy these files while running on worker nodes.
#' 
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when 
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register} 
#' will be used.
#' @param ... Additional optional arguments
#' 
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
#'  \item{s }{The optimal filtering threshold as identified by the global similarity index}
#'  
#'  \item{indexValues }{A matrix of dimension (\code{s.len} x 2) giving the tested filtering thersholds and the
#'  corresponding global similarity indices. Note that the threshold values are equally spaced on the \emph{log}
#'  scale, and thus unequally spaced on the count scale (i.e., we test more threshold values at very low levels
#'  of expression, and fewer at very high levels of expression).}
#'  
#'  \item{normFactor }{A vector of length \emph{n} giving the estimated library sizes estimated by the
#'   normalization method specified in \code{normalization}}
#' 
#'  \item{removedData }{A matrix containing the filtered data}
#' }
#' 
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
#' @author 
#' Andrea Rau, Melina Gallopin, Gilles Celeux, and Florence Jaffrezic
#' 
#' @example /inst/examples/HTSFilter-package.R
#' @keywords methods
#' @rdname HTSFilter
#' @importFrom BiocParallel bpmapply bpparam
#' @import methods

setGeneric(
  name = "HTSFilter",
  def = function(x, ...) 
    {standardGeneric("HTSFilter")}
)



#' @rdname HTSFilter
## matrix
setMethod(
  f= "HTSFilter",
  signature = signature(x="matrix"),
  definition = function(x, conds, s.min=1, s.max=200, s.len=100, 
                        loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
                        plot=TRUE, plot.name=NA, parallel=FALSE, BPPARAM=bpparam()) 
    
  {
    normalization <- match.arg(normalization)
    data <- x
    
    ## Multiple factor designs
    if(is.null(dim(conds)) == FALSE) conds <- do.call(paste, as.data.frame(conds))
    
    ## Run filter
    filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
                                   s.max=s.max, s.len=s.len, loess.span=loess.span,
                                   normalization=normalization, plot=plot, 
                                   plot.name=plot.name, parallel=parallel, BPPARAM=BPPARAM)
    
    filteredData <- filter$data.filter
    if(length(rownames(data))) rownames(filteredData) <- rownames(data)[which(filter$on == 1)]
    if(length(colnames(data))) colnames(filteredData) <- colnames(data)
    
    ## Return various results
    filter.results <- list(filteredData = filteredData,
                           on = filter$on, s = filter$s.optimal,
                           indexValues = filter$index.values, normFactor = filter$norm.factor,
                           removedData = data[which(filter$on == 0),])
    
    return(filter.results)
  }
)

#' @rdname HTSFilter
## data.frame
setMethod(
  f= "HTSFilter",
  signature = signature(x="data.frame"),
  definition = function(x, conds, s.min=1, s.max=200, s.len=100, 
                        loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
                        plot=TRUE, plot.name=NA, parallel=FALSE, BPPARAM=bpparam()) 
    
  {
    normalization <- match.arg(normalization)
    data <- as.matrix(x)
    if(is.character(data) == TRUE) {
      stop("Character values detected in data.frame x.\nPlease check that ID names are not included in one of the columns.")
    }
    
    ## Multiple factor designs
    if(is.null(dim(conds)) == FALSE) conds <- do.call(paste, as.data.frame(conds))
    
    ## Run filter
    filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
                                   s.max=s.max, s.len=s.len, loess.span=loess.span,
                                   normalization=normalization, plot=plot, 
                                   plot.name=plot.name, parallel=parallel, BPPARAM=BPPARAM)
    
    filteredData <- filter$data.filter
    if(length(rownames(data))) rownames(filteredData) <- rownames(data)[which(filter$on == 1)]
    if(length(colnames(data))) colnames(filteredData) <- colnames(data)
    
    
    ## Return various results
    filter.results <- list(filteredData =  filteredData,
                           on = filter$on, s = filter$s.optimal,
                           indexValues = filter$index.values, normFactor = filter$norm.factor,
                           removedData = data[which(filter$on == 0),])
    
    return(filter.results)
  }
)


#' @rdname HTSFilter
## DGEList
setMethod(
  f="HTSFilter",
  signature = signature(x="DGEList"),
  definition = function(x, s.min=1, s.max=200, s.len=100,
                        loess.span=0.3, normalization=c("TMM","DESeq","pseudo.counts","none"),
                        plot=TRUE, plot.name=NA, parallel=FALSE, BPPARAM=bpparam(), conds)
  {
    
    if(missing(conds)) conds <- NULL
    if(!require(edgeR)) {
      stop("edgeR library must be installed.")
    }
    
    if(is.null(x$common.dispersion) == FALSE & is.null(x$tagwise.dispersion == TRUE)) {
      stop("Filtering must be performed either before calling estimateCommonDisp, or after calling estimateTagwiseDisp.")
    }
    
    normalization <- match.arg(normalization)
    data <- x$counts
    if(!length(conds))  conds <- x$samples$group
    if(normalization == "pseudo.counts") {
      if(is.null(x$pseudo.counts) == TRUE) {
        stop(paste("To use pseudo.counts for filtering, you must first call estimateCommonDisp."))
      }
      data <- x$pseudo.counts
      normalization <- "none"
    }
    
    ## Run filter
    filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
                                   s.max=s.max, s.len=s.len, loess.span=loess.span,
                                   normalization=normalization, plot=plot, 
                                   plot.name=plot.name,  parallel=parallel, BPPARAM=BPPARAM)
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
    
    filter.results <- list(filteredData = filteredData,
                           on = filter$on, s = filter$s.optimal,
                           indexValues = filter$index.values, 
                           normFactor = nf,
                           removedData = x$counts[-on.index,])
    
    return(filter.results)
  }
)

#' @rdname HTSFilter
## DGEExact
setMethod(
  f="HTSFilter",
  signature = signature(x="DGEExact"),
  definition = function(x, DGEList, s.min=1, s.max=200, s.len=100,
                        loess.span=0.3, normalization=c("TMM","DESeq","pseudo.counts","none"),
                        plot=TRUE, plot.name=NA, parallel=FALSE, BPPARAM=bpparam(), conds)
  {
    
    if(missing(conds)) conds <- NULL

    if(!require(edgeR)) {
      stop("edgeR library must be installed.")
    }
    normalization <- match.arg(normalization)
    data <- DGEList$counts
    if(!length(conds)) conds <- DGEList$samples$group
    if(normalization == "pseudo.counts") {
      if(is.null(DGEList$pseudo.counts) == TRUE) {
        stop(paste("To use pseudo.counts for filtering, you must first call estimateCommonDisp."))
      }
      data <- DGEList$pseudo.counts
      normalization <- "none"
    }
    
    ## Run filter
    filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
                                   s.max=s.max, s.len=s.len, loess.span=loess.span,
                                   normalization=normalization, plot=plot, 
                                   plot.name=plot.name,  parallel=parallel, BPPARAM=BPPARAM)
    on <- filter$on
    on.index <- which(on == 1) 
    
    ## Create a new DGEExact
    filteredData <- x
    filteredData$table <- x$table[on.index,]
    filteredData$genes <- x$genes[on.index]
    
    ## Return various results
    nf <- filter$norm.factor
    if(normalization == "pseudo.counts") nf <- "pseudo.counts"
    filter.results <- list(filteredData = filteredData,
                           on = filter$on, s = filter$s.optimal,
                           indexValues = filter$index.values, 
                           normFactor = nf,
                           removedData = x$counts[-on.index,])
    
    return(filter.results)
  }
)



#' @rdname HTSFilter
## DGEGLM
setMethod(
  f="HTSFilter",
  signature = signature(x="DGEGLM"),
  definition = function(x, s.min=1, s.max=200, s.len=100,
                        loess.span=0.3, normalization=c("TMM","DESeq","none"),
                        plot=TRUE, plot.name=NA, parallel=FALSE, BPPARAM=bpparam(), conds)
  {
    
    if(missing(conds)) conds <- NULL

    if(!require(edgeR)) {
      stop("edgeR library must be installed.")
    }
    normalization <- match.arg(normalization)
    data <- x$counts
    if(!length(conds)) {
      conds <- x$design
      ## Remove intercept if present
      if(colnames(conds)[1] == "(Intercept)") conds <- conds[,-1]
      ## Multiple factor designs
      if(is.null(dim(conds)) == FALSE) conds <- do.call(paste, as.data.frame(conds))
    }
    
    ## Run filter
    filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
                                   s.max=s.max, s.len=s.len, loess.span=loess.span,
                                   normalization=normalization, plot=plot, 
                                   plot.name=plot.name,  parallel=parallel, BPPARAM=BPPARAM)
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
    filter.results <- list(filteredData = filteredData,
                           on = filter$on, s = filter$s.optimal,
                           indexValues = filter$index.values, 
                           normFactor = nf,
                           removedData = x$counts[-on.index,])
    
    return(filter.results)
  }
)

#' @rdname HTSFilter
## DGELRT
setMethod(
  f="HTSFilter",
  signature = signature(x="DGELRT"),
  definition = function(x, DGEGLM, s.min=1, s.max=200, s.len=100,
                        loess.span=0.3, normalization=c("TMM","DESeq","none"),
                        plot=TRUE, plot.name=NA, parallel=FALSE, BPPARAM=bpparam(), conds)
  {
    
    if(missing(conds)) conds <- NULL

    if(!require(edgeR)) {
      stop("edgeR library must be installed.")
    }
    normalization <- match.arg(normalization)
    data <- DGEGLM$counts
    if(!length(conds)) {
      conds <- x$design
      ## Remove intercept if present
      if(colnames(conds)[1] == "(Intercept)") conds <- conds[,-1]
      ## Multiple factor designs
      if(is.null(dim(conds)) == FALSE) conds <- do.call(paste, as.data.frame(conds))
    }
    
    ## Run filter
    filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
                                   s.max=s.max, s.len=s.len, loess.span=loess.span,
                                   normalization=normalization, plot=plot, 
                                   plot.name=plot.name,  parallel=parallel, BPPARAM=BPPARAM)
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
    filter.results <- list(filteredData = filteredData,
                           on = filter$on, s = filter$s.optimal,
                           indexValues = filter$index.values, 
                           normFactor = nf,
                           removedData = x$counts[-on.index,])
    
    return(filter.results)
  }
)



#' @rdname HTSFilter
## DESeqDataSet
setMethod(
  f= "HTSFilter",
  signature = signature(x="DESeqDataSet"),
  definition = function(x, s.min=1, s.max=200, s.len=100, 
                        loess.span=0.3, normalization=c("DESeq", "TMM", "none"), 
                        plot=TRUE, plot.name=NA, pAdjustMethod="BH", 
                        parallel=FALSE, BPPARAM=bpparam(), conds) 
    
  {
    if(missing(conds)) conds <- NULL

    if(!require(DESeq2)) {
      stop("DESeq2 library must be installed.")
    }          
    normalization <- match.arg(normalization)
    data <- counts(x, normalized=FALSE)
    
    ## November 14, 2017: adjust for more complicated designs in DESeq2
    ## especially when not all columns of colData are used in design
    ## Thanks to Stephanie Durand for finding this bug!
    if(!length(conds)) {
      condsMat <- x@colData
      des <- x@design
      desfactors <- strsplit(as.character(x@design)[-1], split=" + ", fixed=TRUE)[[1]]
      indfactors <- strsplit(as.character(x@design)[-1], split=" * ", fixed=TRUE)[[1]]
      desfactors <- unique(c(desfactors, indfactors))
      index <- which(colnames(condsMat) %in% desfactors)
      condsMat <- condsMat[,index]
      ## Multiple factor designs
      if(is.null(dim(condsMat)) == FALSE) condsMat <- do.call(paste, as.data.frame(condsMat))
      conds <- condsMat
    }
    if(min(table(conds)) < 2) {
      print(table(conds))
      stop(paste("Watch out! One of the combinations of factors in your design has a single replicate.
Use the conds argument to specify a condition vector in which all levels are replicated (for instance, a single factor)"))
    }
    
    ## Run filter
    filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
                                   s.max=s.max, s.len=s.len, loess.span=loess.span,
                                   normalization=normalization, plot=plot, 
                                   plot.name=plot.name,  parallel=parallel, BPPARAM=BPPARAM)
    on <- filter$on
    on.index <- which(on == 1) 
    filteredData <- x[on.index,]
    
    ## April 3, 2015: remove this as p-values are not adjusted in DESeqDataSet
    #		if(length(colnames(mcols(filteredData))) > 0) {
    #                ## Re-adjust p-values, if they exist
    #                nm <- strsplit(colnames(mcols(filteredData)), split="_", fixed=TRUE)
    #                Waldindex <- which(unlist(lapply(nm, function(yy) yy[1]))=="WaldPvalue")
    #                LRTindex <- which(unlist(lapply(nm,  function(yy) yy[1]))=="LRTPvalue")
    #               # Wald p-values
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
                           on = filter$on, s = filter$s.optimal,
                           indexValues = filter$index.values, normFactor = filter$norm.factor,
                           removedData = data[which(filter$on == 0),])
    
    return(filter.results)
  }
)



#' Deprecated functions in package \sQuote{HTSFilter}
#' 
#' These functions are provided for compatibility with older versions of \sQuote{HTSFilter} only,
#' and will be defunct at the next release.
#' 
#' The following functions are deprecated and will be made defunct. Objects of type
#' \sQuote{CountDataSet} from the original \sQuote{DESeq} package will no longer be supported;
#' users should make use of \sQuote{DESeqDataSet} objects from the \sQuote{DESeq2} package:
#' \itemize{
#'   \item{HTSFilter.CountDataSet: \code{\link{HTSFilter}}}
#'   \item{HTSBasicFilter.CountDataSet: \code{\link{HTSBasicFilter}}} 
#' }
#' 
#' @name HTSFilter-deprecated
#' @aliases 
#' HTSFilter-deprecated
#' HTSFilter,CountDataSet-method
#' HTSBasicFilter,CountDataSet-method
#' 
#' @export
#' 
#' @param x A numeric matrix or data.frame representing the counts of dimension (\emph{g} x \emph{n}), 
#' for \emph{g} genes in \emph{n} samples, a \code{DGEList} object, a 
#' \code{DGEExact} object, a \code{DGEGLM} object, a \code{DGELRT} object, or a \code{DESeqDataSet} object.
#' 
#' @param conds  Vector of length \emph{n} identifying the experimental condition of each of the \emph{n} samples; required when sQuote(x)
#' is a numeric matrix.
#' 
#' @param s.min Minimum value of filtering threshold to be considered, with default value equal to 1
#' 
#' @param s.max Maximum value of filtering threshold to be considered, with default value equal to 200
#' 
#' @param s.len Length of sequence of filtering thresholds to be considered (from \code{s.min} to \code{s.max}) 
#'   for the calculation of the global similarity index
#' 
#' @param loess.span Span of the loess curve to be fitted to the filtering thresholds and corresponding global similarity
#'   indices, with default value equal to 0.3
#' 
#' @param normalization Normalization method to be used to correct for differences in library sizes, with 
#' choices  \dQuote{TMM} (Trimmed Mean of M-values), \dQuote{DESeq} (normalization method proposed in the
#' DESeq package), \dQuote{pseudo.counts} (pseudo-counts obtained via quantile-quantile normalization in 
#' the edgeR package, only available for objects of class \code{DGEList} and \code{DGEExact}), and 
#' \dQuote{none} (to be used only if user is certain no normalization is required, or if data have already 
#' been pre-normalized by an alternative method)
#'
#' @param plot If \dQuote{TRUE}, produce a plot of the calculated global similarity indices against the
#' filtering threshold with superimposed loess curve
#' 
#' @param plot.name  If \code{plot} = \dQuote{TRUE}, the name of the PDF file to be saved to the current working directory.
#'   If \code{plot.name} = NA, the plot is drawn in the current window
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
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel 
#' execution using BiocParallel (see next argument \code{BPPARAM}). A note on running 
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects 
#' from the current R environment before calling the function, as it is possible that R's 
#' internal garbage collection will copy these files while running on worker nodes.
#' 
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when 
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register} 
#' will be used.
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
#'  \item{s }{The optimal filtering threshold as identified by the global similarity index}
#'  
#'  \item{indexValues }{A matrix of dimension (\code{s.len} x 2) giving the tested filtering thersholds and the
#'  corresponding global similarity indices. Note that the threshold values are equally spaced on the \emph{log}
#'  scale, and thus unequally spaced on the count scale (i.e., we test more threshold values at very low levels
#'  of expression, and fewer at very high levels of expression).}
#'  
#'  \item{normFactor }{A vector of length \emph{n} giving the estimated library sizes estimated by the
#'   normalization method specified in \code{normalization}}
#' 
#'  \item{removedData }{A matrix containing the filtered data}
#' }
#' 
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
#' @author 
#' Andrea Rau, Melina Gallopin, Gilles Celeux, and Florence Jaffrezic
#' @keywords methods
#' @importClassesFrom DESeq CountDataSet
#' @importMethodsFrom DESeq counts sizeFactors
#' @importMethodsFrom Biobase pData
## CountDataSet
setMethod(
  f="HTSFilter",
  signature = signature(x="CountDataSet"),
  definition = function(x, conds=NA, s.min=1, s.max=200, s.len=100,
                        loess.span=0.3, normalization=c("DESeq","TMM","none"),
                        plot=TRUE, plot.name=NA, parallel=FALSE, BPPARAM=bpparam())
  {
    .Deprecated("HTSFilter", 
                msg="HTSFilter will soon no longer support CountDataSet objects from the original DESeq package.")
    
    norm <- normalization <- match.arg(normalization)
    if(!require(DESeq)) {
      stop("DESeq library must be installed.")
    }
    
    if(is.na(conds)[1] == TRUE) {
      conds <- pData(x)[,-which(colnames(pData(x)) == "sizeFactor")];
      ## Multiple factor designs
      if(is.null(dim(conds)) == FALSE) conds <- do.call(paste, as.data.frame(conds))
    }
    if(is.na(conds)[1] == FALSE & is.null(dim(conds)) == FALSE) {
      ## Multiple factor designs
      conds <- do.call(paste, as.data.frame(conds))
    }
    
    ## What if alternative normalization is desired in filter?
    if(norm == "TMM") {
      message("NOTE: use of TMM normalization in filter for S4 object of class CountDataSet.")
    }
    if(norm == "none") {
      message("NOTE: use of no normalization in filter for S4 object of class CountDataSet.")
    }
    
    opt <- nf <- NA
    ## Option 1: filter, lib size, disp
    if(sum(is.na(sizeFactors(x))) > 0) {
      data <- counts(x)
      opt <- 1
    }
    ## Option 2 or 3: lib size, filter, disp / lib size, disp, filter
    if(sum(is.na(sizeFactors(x))) == 0) {
      if(norm == "DESeq") {
        data <- counts(x, normalized = TRUE)
        normalization <- "none"
        nf <- sizeFactors(x)
      }
      if(norm != "DESeq") {
        data <- counts(x)
      }
      opt <- 2
      if(length(ls(x@fitInfo)) > 0) opt <- 3
    }
    
    ## Run filter
    filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
                                   s.max=s.max, s.len=s.len, loess.span=loess.span,
                                   normalization=normalization, plot=plot, plot.name=plot.name,
                                   parallel=parallel, BPPARAM=BPPARAM)
    on <- filter$on
    on.index <- which(on == 1) 
    filteredData <- x[on.index,]
    
    #		## Option 1: filter, lib size, disp
    #		## Option 2: lib size, filter, disp
    #		if(opt == 1 | opt == 2) {
    #			filteredData <- assayDataElementReplace(filteredData, "counts", 
    #				assayData(x)[["counts"]][on.index,])
    #			featureData(filteredData)@data <- featureData(x)@data[on.index,]
    #			## Sanity check
    #			if(validObject(filteredData)!=TRUE) {
    #				stop(paste(sQuote("filteredData"), 
    #					"is not a valid CountDataSet object."))
    #			}
    #		}
    #
    #		## Option 3: lib size, disp, filter
    #		if(opt == 3) {
    #			filteredData <- assayDataElementReplace(filteredData, "counts", 
    #				assayData(x)[["counts"]][on.index,])		
    #			disp <- data.frame(featureData(x)@data[on.index,])
    #			rownames(disp) <- rownames(featureData(x)@data)[on.index]
    #			colnames(disp) <- colnames(featureData(x)@data)
    #			featureData(filteredData)@data <- disp
    #			## Sanity check
    #			if(validObject(filteredData)!=TRUE) {
    #				stop(paste(sQuote("filteredData"), 
    #					"is not a valid CountDataSet object."))
    #			}
    #		}
    
    nf <- filter$norm.factor
    if(opt != 1 & norm == "DESeq") nf <- sizeFactors(x)
    ## Return various results
    filter.results <- list(filteredData = filteredData,
                           on = filter$on, s = filter$s.optimal,
                           indexValues = filter$index.values, 
                           normFactor = nf,
                           removedData = x[-on.index,])
    
    return(filter.results)
  }
)


#' @rdname HTSFilter-deprecated
## CountDataSet
setMethod(
  f="HTSBasicFilter",
  signature = signature(x="CountDataSet"),
  definition = function(x, method, cutoff.type="value", cutoff=10, 
                        length=NA, normalization=c("DESeq", "TMM", "none")) 
  {
    norm <- normalization <- match.arg(normalization)
    if(!require(DESeq)) {
      stop("DESeq library must be installed.")
    }
    
    ## What if alternative normalization is desired in filter?
    if(norm == "TMM") {
      message("NOTE: use of TMM normalization in filter for S4 object of class CountDataSet.")
    }
    if(norm == "none") {
      message("NOTE: use of no normalization in filter for S4 object of class CountDataSet.")
    }
    
    opt <- nf <- NA
    ## Option 1: filter, lib size, disp
    if(sum(is.na(sizeFactors(x))) > 0) {
      data <- counts(x)
      opt <- 1
    }
    ## Option 2 or 3: lib size, filter, disp / lib size, disp, filter
    if(sum(is.na(sizeFactors(x))) == 0) {
      if(norm == "DESeq") {
        data <- counts(x, normalized = TRUE)
        normalization <- "none"
        nf <- sizeFactors(x)
      }
      if(norm != "DESeq") {
        data <- counts(x)
      }
      opt <- 2
      if(length(ls(x@fitInfo)) > 0) opt <- 3
    }
    
    
    ## Run filter
    filter <- .HTSBasicFilterBackground(data=data, method=method,
                                        cutoff.type=cutoff.type, cutoff=cutoff, length=length,
                                        normalization=normalization)
    on <- filter$on
    on.index <- which(on == 1)
    filteredData <- x[on.index,]
    
    #		## Option 1: filter, lib size, disp
    #		## Option 2: lib size, filter, disp
    #		if(opt == 1 | opt == 2) {
    #			filteredData <- assayDataElementReplace(filteredData, "counts", 
    #				assayData(x)[["counts"]][on.index,])
    #			featureData(filteredData)@data <- featureData(x)@data[on.index,]
    #			## Sanity check
    #			if(validObject(filteredData)!=TRUE) {
    #				stop(paste(sQuote("filteredData"), 
    #					"is not a valid CountDataSet object."))
    #			}
    #		}
    #
    #		## Option 3: lib size, disp, filter
    #		if(opt == 3) {
    #			filteredData <- assayDataElementReplace(filteredData, "counts", 
    #				assayData(x)[["counts"]][on.index,])		
    #			disp <- data.frame(featureData(x)@data[on.index,])
    #			rownames(disp) <- rownames(featureData(x)@data)[on.index]
    #			colnames(disp) <- colnames(featureData(x)@data)
    #			featureData(filteredData)@data <- disp
    #			## Sanity check
    #			if(validObject(filteredData)!=TRUE) {
    #				stop(paste(sQuote("filteredData"), 
    #					"is not a valid CountDataSet object."))
    #			}
    #		}
    
    nf <- filter$norm.factor
    if(opt != 1 & norm == "DESeq") nf <- sizeFactors(x)
    
    ## Return various results
    filter.results <- list(filteredData =  filteredData,
                           on = filter$on, normFactor = nf,
                           removedData = x[-on.index,], filterCrit = filter$filterCrit)
    return(filter.results)
  }
)

