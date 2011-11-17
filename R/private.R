.covMatConstruction <- function(sigma2s, rho, type="unif"){
	###########################################################
	# function that calculates the covariance matrix 
	###########################################################
	if (type=="unif"){
		cov.mat <- matrix(rho, nrow=length(sigma2s), ncol=length(sigma2s)) + diag(1-rho, length(sigma2s), length(sigma2s))
		cov.mat <- diag(sqrt(sigma2s)) %*% cov.mat %*% diag(sqrt(sigma2s))
	} 
	if (type=="ar1"){
		rho.seq <- numeric()
		for (j in 1:length(sigma2s)){
			rho.seq <- c(rho.seq, rho^(j-1))
		}
		cov.mat <- rbind(as.numeric(rho.seq), t(sapply(c(2:length(sigma2s)), function(j, x){ c(x[(j):2], x[1:(length(x)-j+1)]) }, x=rho.seq)))
		cov.mat <- diag(sqrt(sigma2s)) %*% cov.mat %*% diag(sqrt(sigma2s))
	} 
	return(cov.mat)
}


.makeSegments <- function(data) {
	previous    <- 2000
	values      <- c()
	start       <- c()
	end         <- c()
	for (i in 1:length(data)) {
		if (data[i] != previous) {
			start   <- c(start, i)
			last    <- i - 1
			if (last > 0) end <- c(end, last)
			values  <- c(values, data[i])
       		}
		previous  <- data[i]
	}
	end     <- c(end, length(data))
	result  <- cbind(values, start, end)
	result
}


.minDistCghFeature_CGHcall2ExpressionSet <- function(GEfeature.ann, CNann){
	CNfeaturesKept.temp <- which(CNann[,1] == as.numeric(GEfeature.ann[1]))
	if (length(CNfeaturesKept.temp) >= 1){
		first.feature <- min(CNfeaturesKept.temp) - 1
		CNann <- CNann[CNfeaturesKept.temp, , drop=FALSE]
		CNfeature.closest <- which.min(abs(CNann[,2] - GEfeature.ann[2]))[1] + first.feature
		return(CNfeature.closest)
	} else {
		return(NA)
	}
}


.overlapPercentage_CGHcall2ExpressionSet <- function(GEfeature.ann, CNann, reference=1){
	CNfeaturesKept.temp <- which(CNann[,1] == as.numeric(GEfeature.ann[1]))
	first.feature <- min(CNfeaturesKept.temp) - 1
	CNann <- CNann[CNfeaturesKept.temp, , drop=FALSE]
	junk <- sapply(CNann[,3], function(x,y){ min(x,y) }, y=GEfeature.ann[3]) - sapply(CNann[,2], function(x,y){ max(x,y) }, y=GEfeature.ann[2]) + 1
	if (reference==2){
		percOverlap <- as.numeric(sapply(junk, function(x){ max(x, 0) }) / (GEfeature.ann[3] - GEfeature.ann[2] +1))
	}
	if (reference==1){
		percOverlap <- as.numeric(sapply(junk, function(x){ max(x, 0) }) / (CNann[,3] - CNann[,2] +1))
	}
	candidate <- as.numeric(which.max(percOverlap)[1])
	if (percOverlap[candidate] == 0){
		return(NA)
	} else {
		return(candidate + first.feature)
	}
}


.overlapPlusInterpolation_CGHcall2ExpressionSet <- function(GEfeature.ann, CNann, CNsegments, reference=1){
	CNfeaturesKept.temp <- which(CNann[,1] == as.numeric(GEfeature.ann[1]))
	first.feature <- min(CNfeaturesKept.temp) - 1
	CNann <- CNann[CNfeaturesKept.temp, , drop=FALSE]
	CNsegments <- CNsegments[CNfeaturesKept.temp, , drop=FALSE]
	junk <- sapply(CNann[,3], function(x,y){ min(x,y) }, y=as.numeric(GEfeature.ann[3])) - sapply(CNann[,2], function(x,y){ max(x,y) }, y=as.numeric(GEfeature.ann[2])) + 1
	if (reference==2){
		percOverlap <- as.numeric(sapply(junk, function(x){ max(x, 0) }) / (as.numeric(GEfeature.ann[3]) - as.numeric(GEfeature.ann[2]) +1))
	}
	if (reference==1){
		percOverlap <- as.numeric(sapply(junk, function(x){ max(x, 0) }) / (CNann[,3] - CNann[,2] +1))
	}
	percOverlap <- as.numeric(sapply(junk, function(x){ max(x, 0) }) / (as.numeric(GEfeature.ann[3]) - as.numeric(GEfeature.ann[2]) +1))
	candidate <- as.numeric(which.max(percOverlap)[1])
	if (percOverlap[candidate] == 0){
		dists <- cbind(CNann[,3]-as.numeric(GEfeature.ann[3]), as.numeric(GEfeature.ann[2]) - CNann[,2])
		ids.temp <- which(dists[,1] > 0)
		if (length(ids.temp) > 0){ id1 <- ids.temp[which.min(dists[ids.temp,1])] } else { id1 <- which.max(dists[,1]) }
		ids.temp <- which(dists[,2] > 0)
		if (length(ids.temp) > 0){ id2 <- ids.temp[which.min(dists[ids.temp,2])] } else { id2 <- which.max(dists[,2]) }
		ids <- sort(c(id1, id2))
		if ( all(CNsegments[ids,][1,] == CNsegments[ids,][2,]) ){ 
			return(ids[which.min(apply(abs(dists[ids,]), 1, mean))[1]] + first.feature)
		} else {
			return(NA)
		}
	} else {
		return(candidate + first.feature)
	}
}


.minDistFeature <- function(feature2ann, ann1temp, maxDist){
	featuresKept1temp <- which(ann1temp[,1] == as.numeric(feature2ann[1]))
	if (length(featuresKept1temp) >= 1){
		ann1temp <- ann1temp[featuresKept1temp, , drop=FALSE]
		dists <- abs(ann1temp[,2] - feature2ann[2])
		candidates <- which(dists < maxDist)
		if (length(candidates) > 0){
			slh <- cbind(feature2ann[3], ann1temp[candidates, 3], dists[candidates])
			colnames(slh) <- c("platform2", "platform1", "distance")
			rownames(slh) <- NULL
			return(slh)
		} else {
			return(NULL)
		}
	} else {
		return(NULL)
	}
}


.overlapPercentage <- function(feature2ann, ann1temp, minPerc, reference=1){
	featuresKept1temp <- which(ann1temp[,1] == as.numeric(feature2ann[1]))
	ann1temp <- ann1temp[featuresKept1temp, , drop=FALSE]
	junk <- sapply(ann1temp[,3], function(x,y){ min(x,y) }, y=feature2ann[3]) - sapply(ann1temp[,2], function(x,y){ max(x,y) }, y=feature2ann[2]) + 1
	if (reference==2){
		percOverlap <- as.numeric(sapply(junk, function(x){ max(x, 0) }) / (feature2ann[3] - feature2ann[2] +1))
	}
	if (reference==1){
		percOverlap <- as.numeric(sapply(junk, function(x){ max(x, 0) }) / (ann1temp[,3] - ann1temp[,2] +1))
	}
	candidates <- which(percOverlap > minPerc)
	if (length(candidates) > 0){
			slh <- cbind(feature2ann[4], ann1temp[candidates, 4], percOverlap[candidates])
			colnames(slh) <- c("platform2", "platform1", "overlap")
			rownames(slh) <- NULL
			return(slh)
	} else {
		return(NULL)
	}
}


.slhFuncWeightedCGHcallExpressionSetData <- function(featuresWithWeights, copynumberMat){
	return(as.numeric(matrix(featuresWithWeights[,2,drop=TRUE], nrow=1) %*% copynumberMat[featuresWithWeights[,1,drop=TRUE], , drop=FALSE]) / sum(featuresWithWeights[,2,drop=TRUE]))
}


.slhFuncWeightedCGHcallExpressionSetAnnotation <- function(featuresWithWeights, CNann, chr, bpstart, bpend){
	fAnn <- CNann[featuresWithWeights[, 1, drop=TRUE], , drop=FALSE]
	return(as.numeric(c(fAnn[1, chr], min(fAnn[, bpstart]), max(fAnn[, bpend]))))
}


.slhFuncCombineFeatureNames <- function(featuresWithWeights, fNames){
	return(as.character(paste(fNames[featuresWithWeights[,1,drop=TRUE]], collapse="_")))
}


.slhFuncMaxCGHcallData <- function(featuresWithWeights, segmentedMat){
	return(as.numeric(apply(abs(segmentedMat[featuresWithWeights[,1,drop=TRUE], , drop=FALSE]), 2, which.max)))
}


.slhFunctMaxFeatures <- function(j, featuresAndWeights, idMat, segmentedMat){
	return(as.numeric(segmentedMat[featuresAndWeights[[j]][,1], , drop=FALSE][cbind(idMat[j,], c(1:ncol(segmentedMat)))]))
}



.matrixRows2list <- function(X){
	subList <- sapply(1:nrow(X), function(u, Z){ Z[u, , drop=FALSE] }, Z=X, simplify=FALSE)
	names(subList) <- rep(X[1, 1], nrow(X))
	return(subList)
}
