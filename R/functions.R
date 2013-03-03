.perConditionSimilarityIndex <-
function(data.norm.perCondition, log.s){

	countsNormLog <- log(data.norm.perCondition)
	countsNormLog_binary <- countsNormLog > log.s
	sum <- 0
	nbindiv <- ncol(countsNormLog)
	calc <- mapply(function(i,j) {
		x <- countsNormLog_binary[,i]
		y <- countsNormLog_binary[,j]
		a <- length(which(x == 1 & y == 1))
		b <- length(which(x == 1 & y == 0))
		c <- length(which(x == 0 & y == 1))
		d <- length(which(x == 0 & y == 0))
		## Sanity check
		if(sum(c(a,b,c,d)) != length(x) & sum(c(a,b,c,d)) != length(y)) {
			stop("Something's wrong.")
		}
		a / (a+b+c)
		}, 
		rep(1:(nbindiv), times = c((nbindiv-1):0)),
		unlist(lapply(2:nbindiv, function(x) seq(x,nbindiv))) )
	## AR (1/15/2013): use mean rather than sum in case of unbalanced experiments
	return(mean(calc))
}


.HTSFilterBackground <- 
function(data, conds, s.min, s.max, s.len, 
	loess.span, normalization,  plot, plot.name) {

	## Sanity checks: replicated data?
	if(min(table(conds) == 1)) {
		stop("All conditions must have at least TWO replicates.")
	}
	if(s.min <= 0) stop("s.min must be > 0.");
	if(s.min > s.max) stop(paste(sQuote(s.min), "must be <", sQuote(s.max)));
	if(dim(data)[2] != length(conds)) {
		stop(paste(sQuote(conds), " must be of length ", dim(data)[2], 
		" (the number of columns of", sQuote(data), ")", sep = ""))
	}

	## Normalization (calculated on full dataset)
	norm <- normalizeData(data = data,  
		normalization = normalization)
	data.norm <- norm$data.norm
	norm.factor <- norm$norm.factor

	## Calculate index for each threshold value
	s.test <- seq(log(s.min), log(s.max), length = s.len)
	## AR (1/17/2013): use mean rather than sum so max value is 1
	index <- rowMeans(matrix(mapply(function(condition, s) {
		.perConditionSimilarityIndex(data.norm.perCondition = data.norm[,which(conds == unique(conds)[condition])], 
		log.s = s)}, 
		rep(1:length(unique(conds)), times = length(s.test)), 
		rep(s.test, each = length(unique(conds)))), ncol = length(unique(conds)), 
		byrow = TRUE))

	## Choosing s.optimal: Fit a loess curve to data
	index[which(is.nan(index) == TRUE)] <- NA
	p <- predict(loess(index ~ exp(s.test), span = loess.span), 
		seq(s.min, s.max, by = 0.001))
	s.optimal <- seq(s.min, s.max, by = 0.001)[which(p == max(p, na.rm = TRUE))]
	if(plot == TRUE) {
		if(is.na(plot.name) != TRUE) pdf(plot.name, width = 6, height = 6);
		plot(exp(s.test), index, log = "x", xlab = "Threshold", 
			ylab = "Global Jaccard index")
		lines(seq(s.min, s.max, by = 0.001), p, col = "blue", lwd = 2)
		points(s.optimal, max(p, na.rm = TRUE), col = "red", pch = "X", cex= 1.5)
		abline(v = s.optimal, lty = 2, col= "red")
		legend("bottomleft", pch = "X", col="red", cex = 1.25, paste("s =", s.optimal), bty = "n")
		if(is.na(plot.name) != TRUE) dev.off();
	}

	## Calculate which genes are present in at least one sample with s.optimal
	data.binary <- matrix(NA, nrow = length(conds), ncol = dim(data)[1])
	data.binary <- data.norm > s.optimal
	on <- ifelse(rowSums(data.binary) == 0, 0, 1)

	data.filter <- data[which(on == 1),]
	colnames(data.binary) <- colnames(data)
	rownames(data.binary) <- rownames(data)

	filter.results <- list(on = on, s.optimal = s.optimal, data.binary = data.binary,
		index.values = data.frame(threshold = exp(s.test), index = index), 
		data.filter = data.filter, norm.factor = norm.factor)

	return(filter.results)
}


.HTSBasicFilterBackground <- function(data, method, cutoff.type="value", cutoff=10, 
	length=NA, normalization) {

	## Sanity checks
	check <- method %in% c("mean", "sum", "rpkm", "variance", "cpm",
		"max", "cpm.mean", "cpm.sum", "cpm.variance", "cpm.max",
		"rpkm.mean", "rpkm.sum", "rpkm.variance", "rpkm.max")
	if(check != TRUE) stop("Only the following basic filters are currently supported:
		mean, sum, rpkm, variance, cpm, max")
	if(class(cutoff.type) != "character" & class(cutoff.type) != "numeric" & 
		length(cutoff.type) != 1)
		stop(paste("cutoff.type must be equal to a numeric value, or one of the following:\n",
			dQuote("value"), dQuote("number"), dQuote("quantile")))  

	## Normalize data
	x <- data
	norm <- normalizeData(x, normalization)
	x.norm <- norm$data.norm
	norm.factor <- norm$norm.factor

	## Prep filtering criteria
	if(method == "mean") crit <- apply(x.norm, 1, mean);
	if(method == "sum") crit <- apply(x.norm, 1, sum);
	if(method == "variance") crit <- apply(x.norm, 1, var);
	if(method == "max") crit <- apply(x.norm, 1, max);
	if(method == "cpm" | method == "cpm.mean" | method == "cpm.sum" | 
		method == "cpm.variance" | method == "cpm.max") {
		if(normalization != "TMM") message("Note that TMM normalization is used for cpm filter.")
		dge <- DGEList(counts=x)
		dge <- calcNormFactors(dge)
		crit <- .cpmDGEList(dge, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
		if(method == "cpm.mean") crit <- apply(crit, 1, mean)
		if(method == "cpm.sum") crit <- apply(crit, 1, sum)
		if(method == "cpm.variance") crit <- apply(crit, 1, var)
		if(method == "cpm.max") crit <- apply(crit, 1, max)
	}
	if(method == "rpkm" | method == "rpkm.mean" | method == "rpkm.sum" | 
		method == "rpkm.variance" | method == "rpkm.max") {
		if(is.vector(length) == FALSE | length(length) != nrow(x)) 
			stop(paste("length needed for rpkm filter"))
		if(normalization != "TMM") message("Note that TMM normalization is used for rpkm filter.")
		dge <- DGEList(counts=x)
		dge <- calcNormFactors(dge)
		crit <- .rpkmDGEList(dge, length, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25)
		if(method == "rpkm.mean") crit <- apply(crit, 1, mean)
		if(method == "rpkm.sum") crit <- apply(crit, 1, sum)
		if(method == "rpkm.variance") crit <- apply(crit, 1, var)
		if(method == "rpkm.max") crit <- apply(crit, 1, max)
	}

	## Apply filters
	if(method == "mean" | method == "sum" | method == "variance" | method == "max" |
		method == "cpm.mean" | method == "cpm.sum" | method == "cpm.variance" | 
		method == "cpm.max" | method == "rpkm.mean" | method == "rpkm.sum" | 
		method == "rpkm.variance" | method == "rpkm.max") {
		if(class(cutoff.type) == "numeric")
			stop(paste("cutoff.type must be equal to one of the following:\n",
			dQuote("value"), dQuote("number"), dQuote("quantile")))  

		if(cutoff.type == "value") {
			on.index <- which(crit > cutoff)
		}
		if(cutoff.type == "number") {
			o <- order(crit, decreasing=TRUE)
			on.index <- which(o <= cutoff)
		}
		if(cutoff.type == "quantile") {
			q <- quantile(crit, cutoff, na.rm=TRUE)
			on.index <- which(crit > q)
		}
	}
	if(method == "cpm" | method == "rpkm") {
		if(class(cutoff.type) == "character")
			stop(paste("cutoff.type must be numeric.")) 
		on.index <- rowSums(crit>cutoff) >= cutoff.type
		on.index <- which(on.index == TRUE)
	}
	if(method =="rpkm.mean" | method == "rpkm.sum" | method == "rpkm.variance" |
		method == "rpkm.max") {
		no.length <- which(is.na(crit) == TRUE)
		on.index <- sort(c(on.index, no.length))	
	}
	if(method == "rpkm") {
		no.length <- which(is.na(crit[,1]) == TRUE)
		on.index <- sort(c(on.index, no.length))	
	}

	## Return filter results
	filteredData <- x; removedData <- NA;
	if(length(on.index) > 0) filteredData <- x[on.index,]
	if(length(on.index) < nrow(x)) removedData <- x[-on.index,]
	on <- rep(0, nrow(x)); on[on.index] <- 1
	filter.results <- list(filteredData =  filteredData,
		on = on, normFactor = norm.factor, removedData = removedData, 
		filterCrit = crit)
	return(filter.results)
} 


## RPKM function taken from edgeR version 3.1.3
.rpkmDGEList <- function(x, gene.length, normalized.lib.sizes=TRUE, log=FALSE, 
    prior.count = 0.25) 
{
    y <- .cpmDGEList(x, normalized.lib.sizes=normalized.lib.sizes, 
        log=log, prior.count = prior.count)
    if (log) 
        y - log2(gene.length) + log2(1000)
    else y/(gene.length/1000)
}

## CPM functions taken from edgeR version 3.1.3
.cpmDGEList <- function (x, normalized.lib.sizes=TRUE, log=FALSE, prior.count=0.25) 
{
    lib.size <- x$samples$lib.size
    if (normalized.lib.sizes) 
        lib.size <- lib.size * x$samples$norm.factors
    .cpm(x$counts, lib.size=lib.size, log=log, prior.count=prior.count)
}
.cpm <- function (x, lib.size=NULL, log=FALSE, prior.count=0.25) 
{
    x <- as.matrix(x)
    if (is.null(lib.size)) 
        lib.size <- colSums(x)
    if (log) {
        prior.count.scaled <- lib.size/mean(lib.size) * prior.count
        lib.size <- lib.size + prior.count.scaled
    }
    lib.size <- 1e-06 * lib.size
    if (log) 
        log2(t((t(x) + prior.count.scaled)/lib.size))
    else t(t(x)/lib.size)
}
