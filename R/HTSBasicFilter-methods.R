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
		filteredData$logCPM <- x$logCPM[on.index]
		filteredData$tagwise.dispersion <- x$tagwise.dispersion[on.index]

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
		
		## Return various results
		nf <- filter$norm.factor
		filter.results <- list(filteredData =  filteredData,
			on = filter$on, normFactor = nf,
			removedData = filter$removedData, filterCrit = filter$filterCrit)
		return(filter.results)
	}
)


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
		data <- DGEGLM$counts
	
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
		filteredData$genes <- x$genes[on.index]
		filteredData$weights <- x$weights[on.index,]
		filteredData$df.residual <- x$df.residual[on.index]
		filteredData$dispersion <- x$dispersion[on.index]
		filteredData$fitted.values <- x$fitted.values[on.index,]
		filteredData$deviance <- x$deviance[on.index]
		filteredData$abundance <- x$abundance[on.index]
		filteredData$offset <- x$offset[on.index,]
		
		## Return various results
		nf <- filter$norm.factor
		filter.results <- list(filteredData =  filteredData,
			on = filter$on, normFactor = nf,
			removedData = filter$removedData, filterCrit = filter$filterCrit)
		return(filter.results)
	}
)



## DESeqDataSet
setMethod(
	f= "HTSBasicFilter",
	signature = signature(x="DESeqDataSet"),
	definition = function(x, method, cutoff.type="value", cutoff=10,
	length=NA, normalization=c("DESeq", "TMM", "none"))
	{
		if(!require(DESeq2)) {
			stop("DESeq2 library must be installed.")
		}          
		normalization <- match.arg(normalization)
                data <- counts(x, normalized=FALSE)
                if(ncol(colData(x)) > 2) {
                  stop("HTSFilter currently only supports a single factor when working with DESeq2.")
                }

		## Run filter
		filter <- .HTSBasicFilterBackground(data=data, method=method,
			cutoff.type=cutoff.type, cutoff=cutoff, length=length,
			normalization=normalization)
		on <- filter$on
		on.index <- which(on == 1) 
		filteredData <- x[on.index,]
                
                ## Re-adjust p-values
                nm <- strsplit(colnames(mcols(filteredData)), split="_", fixed=TRUE)
                Waldindex <- which(unlist(lapply(nm, function(yy) yy[1]))=="WaldPvalue")
                LRTindex <- which(unlist(lapply(nm,  function(yy) yy[1]))=="LRTPvalue")
                message("Note: BH correction of p-values used in HTSFilter.")
                # Wald p-values
                if(length(Waldindex) > 0 ) {
                  for(j in Waldindex) {
                    look <- substr(colnames(mcols(filteredData))[j], 12, 100)
                    find <- which(substr(colnames(mcols(filteredData)), 12+3, 100) == look)
                    find <- find[which(find > j)]
                    mcols(filteredData)[,find] <- p.adjust(mcols(filteredData)[,j], method="BH")
                  }
                }
                # LRT p-values
                if(length(LRTindex) > 0 ) {
                  for(j in LRTindex) {
                    look <- substr(colnames(mcols(filteredData))[j], 11, 100)
                    find <- which(substr(colnames(mcols(filteredData)), 11+3, 100) == look)
                    find <- find[which(find > j)]
                    mcols(filteredData)[,find] <- p.adjust(mcols(filteredData)[,j], method="BH")
                  }
                }
                
		## Return various results
		filter.results <- list(filteredData = filteredData,
			on = filter$on, normFactor = filter$norm.factor,
			removedData = data[which(filter$on == 0),], filterCrit = filter$filterCrit)
		return(filter.results)
	}
)



