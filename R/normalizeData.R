normalizeData <-
function(data, normalization) {
	if(normalization == "TMM") {
		N <- colSums(data)
		f <- calcNormFactors(data,method="TMM")
		TMM <- N*f / mean(N*f)
		norm.factor <- TMM
		data.norm <- scale(data, center=FALSE, scale=TMM)


	}
	if(normalization == "DESeq") {
		## Code taken from DESeq (v1.8.3)
		## estimateSizeFactorsForMatrix() function:
    		loggeomeans <- rowMeans(log(data))
   		deseq <- apply(data, 2, function(cnts) exp(median((log(cnts) - 
       			loggeomeans)[is.finite(loggeomeans)])))
#		deseq <- estimateSizeFactorsForMatrix(data, locfunc = median)
		norm.factor <- deseq
		data.norm <- scale(data, center=FALSE, scale=deseq)
	}
	if(normalization == "none") {
		data.norm <- data
		norm.factor <- NA
	}
	return(list(data.norm = data.norm, norm.factor = norm.factor))
}
