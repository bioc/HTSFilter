## matrix
setMethod(
	f= "HTSFilter",
	signature = signature(x="matrix"),
	definition = function(x, conds, s.min=1, s.max=200, s.len=100, 
		loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
		plot=TRUE, plot.name=NA) 

	{
		normalization <- match.arg(normalization)
		data <- x

		## Multiple factor designs
		if(is.null(dim(conds)) == FALSE) conds <- do.call(paste, as.data.frame(conds))

		## Run filter
		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
			s.max=s.max, s.len=s.len, loess.span=loess.span,
			normalization=normalization, plot=plot, plot.name=plot.name)

		## Return various results
		filter.results <- list(filteredData =  filter$data.filter,
			on = filter$on, s = filter$s.optimal,
			indexValues = filter$index.values, normFactor = filter$norm.factor,
			removedData = data[which(filter$on == 0),])

		return(filter.results)
	}
)

## data.frame
setMethod(
	f= "HTSFilter",
	signature = signature(x="data.frame"),
	definition = function(x, conds, s.min=1, s.max=200, s.len=100, 
		loess.span=0.3, normalization=c("TMM", "DESeq", "none"), 
		plot=TRUE, plot.name=NA) 

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
			normalization=normalization, plot=plot, plot.name=plot.name)

		## Return various results
		filter.results <- list(filteredData =  filter$data.filter,
			on = filter$on, s = filter$s.optimal,
			indexValues = filter$index.values, normFactor = filter$norm.factor,
			removedData = data[which(filter$on == 0),])

		return(filter.results)
	}
)


## CountDataSet
setMethod(
	f="HTSFilter",
	signature = signature(x="CountDataSet"),
	definition = function(x, conds=NA, s.min=1, s.max=200, s.len=100,
 		loess.span=0.3, normalization=c("DESeq","TMM","none"),
 		plot=TRUE, plot.name=NA)
 	{
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
 			normalization=normalization, plot=plot, plot.name=plot.name)
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

## DGEList
setMethod(
	f="HTSFilter",
	signature = signature(x="DGEList"),
	definition = function(x, s.min=1, s.max=200, s.len=100,
		loess.span=0.3, normalization=c("TMM","DESeq","pseudo.counts","none"),
		plot=TRUE, plot.name=NA)
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
		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
			s.max=s.max, s.len=s.len, loess.span=loess.span,
			normalization=normalization, plot=plot, plot.name=plot.name)
		on <- filter$on
		on.index <- which(on == 1) 
		filteredData <- x

		## Create a new DGEList
		filteredData <- x
		filteredData$counts <- x$counts[on.index,]
		filteredData$pseudo.counts <- x$pseudo.counts[on.index,]
		filteredData$logCPM <- x$logCPM[on.index]
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


## DGEExact
setMethod(
	f="HTSFilter",
	signature = signature(x="DGEExact"),
	definition = function(x, DGEList, s.min=1, s.max=200, s.len=100,
		loess.span=0.3, normalization=c("TMM","DESeq","pseudo.counts","none"),
		plot=TRUE, plot.name=NA)
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
		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
			s.max=s.max, s.len=s.len, loess.span=loess.span,
			normalization=normalization, plot=plot, plot.name=plot.name)
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


## DGEGLM
setMethod(
	f="HTSFilter",
	signature = signature(x="DGEGLM"),
	definition = function(x, s.min=1, s.max=200, s.len=100,
		loess.span=0.3, normalization=c("TMM","DESeq","none"),
		plot=TRUE, plot.name=NA)
	{

		if(!require(edgeR)) {
			stop("edgeR library must be installed.")
		}
		normalization <- match.arg(normalization)
		data <- x$counts
		conds <- x$design
		## Remove intercept if present
		if(colnames(conds)[1] == "(Intercept)") conds <- conds[,-1]
		## Multiple factor designs
		if(is.null(dim(conds)) == FALSE) conds <- do.call(paste, as.data.frame(conds))
	
		## Run filter
		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
			s.max=s.max, s.len=s.len, loess.span=loess.span,
			normalization=normalization, plot=plot, plot.name=plot.name)
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


## DGELRT
setMethod(
	f="HTSFilter",
	signature = signature(x="DGELRT"),
	definition = function(x, DGEGLM, s.min=1, s.max=200, s.len=100,
		loess.span=0.3, normalization=c("TMM","DESeq","none"),
		plot=TRUE, plot.name=NA)
	{

		if(!require(edgeR)) {
			stop("edgeR library must be installed.")
		}
		normalization <- match.arg(normalization)
		data <- DGEGLM$counts
		conds <- x$design
		## Remove intercept if present
		if(colnames(conds)[1] == "(Intercept)") conds <- conds[,-1]
		## Multiple factor designs
		if(is.null(dim(conds)) == FALSE) conds <- do.call(paste, as.data.frame(conds))
	
		## Run filter
		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
			s.max=s.max, s.len=s.len, loess.span=loess.span,
			normalization=normalization, plot=plot, plot.name=plot.name)
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




## DESeqDataSet
setMethod(
	f= "HTSFilter",
	signature = signature(x="DESeqDataSet"),
	definition = function(x, s.min=1, s.max=200, s.len=100, 
		loess.span=0.3, normalization=c("DESeq", "TMM", "none"), 
		plot=TRUE, plot.name=NA, pAdjustMethod="BH") 

	{
		if(!require(DESeq2)) {
			stop("DESeq2 library must be installed.")
		}          
		normalization <- match.arg(normalization)
                data <- counts(x, normalized=FALSE)

		conds <- colData(x)
		if("sizeFactor" %in% colnames(conds)) {
			conds <- conds[,-which(colnames(conds) == "sizeFactor")]
		}		
		## Multiple factor designs
		if(is.null(dim(conds)) == FALSE) conds <- do.call(paste, as.data.frame(conds))

		## Run filter
		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
			s.max=s.max, s.len=s.len, loess.span=loess.span,
			normalization=normalization, plot=plot, plot.name=plot.name)
		on <- filter$on
		on.index <- which(on == 1) 
		filteredData <- x[on.index,]
                
		if(length(colnames(mcols(filteredData))) > 0) {
                	## Re-adjust p-values, if they exist
                	nm <- strsplit(colnames(mcols(filteredData)), split="_", fixed=TRUE)
                	Waldindex <- which(unlist(lapply(nm, function(yy) yy[1]))=="WaldPvalue")
                	LRTindex <- which(unlist(lapply(nm,  function(yy) yy[1]))=="LRTPvalue")
                	# Wald p-values
                	if(length(Waldindex) > 0 ) {
                  	for(j in Waldindex) {
                    		look <- substr(colnames(mcols(filteredData))[j], 12, 100)
                    		find <- which(substr(colnames(mcols(filteredData)), 12+3, 100) == look)
                    		find <- find[which(find > j)]
                    		mcols(filteredData)[,find] <- p.adjust(mcols(filteredData)[,j], method=pAdjustMethod)
                  	}
                	}
                	# LRT p-values
                	if(length(LRTindex) > 0 ) {
                  	for(j in LRTindex) {
                    		look <- substr(colnames(mcols(filteredData))[j], 11, 100)
                    		find <- which(substr(colnames(mcols(filteredData)), 11+3, 100) == look)
                    		find <- find[which(find > j)]
                    		mcols(filteredData)[,find] <- p.adjust(mcols(filteredData)[,j], method=pAdjustMethod)
                  	}
                	}
		}
                
		## Return various results
		filter.results <- list(filteredData = filteredData,
			on = filter$on, s = filter$s.optimal,
			indexValues = filter$index.values, normFactor = filter$norm.factor,
			removedData = data[which(filter$on == 0),])

		return(filter.results)
	}
)





## SeqExpressionSet
#setMethod(
#	f="HTSFilter",
#	signature = signature(x="SeqExpressionSet"),
#	definition = function(x, conds=NA, s.min=1, s.max=200, s.len=100,
# 		loess.span=0.3, normalization=c("TMM","DESeq","none"),
#		plot=TRUE, plot.name=NA)
#	{
#		normalization <- match.arg(normalization)
#		data <- exprs(x)
#		## TO DO: WHAT TO DO IF MORE THAN ONE FACTOR
# 		if(is.na(conds)[1] == TRUE) conds <- pData(x);
# 
#		## Run filter
# 		filter <- .HTSFilterBackground(data=data, conds=conds, s.min=s.min,
# 			s.max=s.max, s.len=s.len, loess.span=loess.span,
# 			normalization=normalization, plot=plot, plot.name=plot.name)
# 
#		## Assign filtered data to matrix
#		filteredData <- newSeqExpressionSet(filter$data.filter,
# 			 phenoData = phenoData(x))
# 
# 		## Return various results
# 		filter.results <- list(filteredData = filteredData,
#			on = filter$on, s = filter$s.optimal,
# 			indexValues = filter$index.values, normFactor = filter$norm.factor,
# 			removedData = data[which(filter$on == 0),])
# 
#		return(filter.results)
# 	}
#)



