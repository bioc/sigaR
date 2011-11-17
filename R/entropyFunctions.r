hdEntropy <- function(Y, method="normal", k=1, center=TRUE, indKnn=TRUE){
	###########################################################
	# calculate entropy for high-dimensional matrix
	###########################################################

	normEntropy <- function(Y){
		###########################################################
		# calculate the entropy (assuming a multivar normal dist)
		###########################################################

		# determine number of samples and features
		n <- nrow(Y)
		p <- ncol(Y)

		# calculate entropy (part 1)
		part1 <- sum(log((apply(Y, 2, sd) * sqrt((n-1) / n))))
		Y <- sweep(Y, 2, apply(Y, 2, sd), "/") * sqrt((n / (n-1)))
 
		# determine shrinkage parameter
		w <- rep(1/n, n)
		lambda <- as.numeric(pvt.get.lambda(x=wt.scale(Y, w, center = TRUE, scale = TRUE), lambda=-1, w=w, verbose=FALSE, type = "correlation", target=0)$lambda, 0)
		print(paste("shrinkage factor:", round(lambda, d=5), sep=" "))
 
		# calculate entropy (part 2)
		part2.slh <- rep(0, p)
		part2.slh[1:(min(dim(Y)))] <- svd(Y %*% t(Y) / dim(Y)[1])$d[1:(min(n, p))]
		part2 <- sum(log((1-lambda)*part2.slh + lambda))/2
		return(part1 + part2 + (1 + log(2*pi)) * p/2)
	}

	knnEntropy <- function(Y, k, indKnn){
		###########################################################
		# calculate the knn entropy (per sample)
		###########################################################
		distMat <- as.matrix(dist(Y))
		rankDistMat <- t(apply(distMat, 1, rank))
		knns <- as.numeric(apply(rankDistMat, 1, function(R, k){ which(R == (k+1)) }, k))
		knnDists <- distMat[cbind(c(1:length(knns)), knns)]
		if (!indKnn){
			return(mean(- (log(k) - log(dim(Y)[1]-1) + lgamma(dim(Y)[2]/2+1) - (dim(Y)[2]/2)*log(pi) - dim(Y)[2]*log(knnDists) ) ))
			} else {
			return(- (log(k) - log(dim(Y)[1]-1) + lgamma(dim(Y)[2]/2+1) - (dim(Y)[2]/2)*log(pi) - dim(Y)[2]*log(knnDists) ) )
		}
	}

	# input checks
	if (as.character(class(Y)) != "matrix"){ stop("Input (Y) is of wrong class.") }
	if (sum(is.na(Y)) != 0){ stop("Matrix Y contains missings.") }
	if (!(method %in% c("normal", "knn"))){ stop("method parameter ill-specified.") }
	if (sum(c(as.character(class(k)) != "numeric", as.character(class(k)) != "integer")) == 0){ stop("Input (k) is of wrong class.") }
	if (k < 1){ stop("k smaller than 1.") }
	if (as.character(class(center)) != "logical"){ stop("center of wrong class.") }
	if (as.character(class(indKnn)) != "logical"){ stop("center of wrong class.") }

	# center expression of all genes around zero
	if(center){ Y <- sweep(Y, 2, apply(Y, 2, mean)) }

	# calculate entropy
	if (method=="normal"){ return(normEntropy(Y)) }
	if (method=="knn"){ return(knnEntropy(Y, k, indKnn)) }
}
	
hdMI <- function(Y, X, method="normal", k=1, center=TRUE, rescale=TRUE){
	###########################################################
	# calculate mutual information for two high-dimensional matrices
	###########################################################

	normMI <- function(Y, X, center){
		###########################################################
		# calculate the entropy (assuming a multivar normal dist)
		###########################################################

		# determine number of samples and features
		n <- nrow(Y)
		pY <- ncol(Y)
		pX <- ncol(X)

		# center data around zero (per row, both matrices)
		if(center){ Y <- sweep(Y, 2, apply(Y, 2, mean)) }
		if(center){ X <- sweep(X, 2, apply(X, 2, mean)) }
 
		# determine shrinkage parameter of joint covariance matrix
		w <- rep(1/n, n)
		lambda <- as.numeric(pvt.get.lambda(x=wt.scale(cbind(Y, X), w, center = TRUE, scale = TRUE), lambda=-1, w=w, verbose=FALSE, type = "correlation", target=0)$lambda, 0)
		print(paste("shrinkage factor:", round(lambda, d=5), sep=" "))
 
		# calculate joint entropy
		part2.slh <- rep(0, pY+pX)
		part2.slh[1:(min(dim(cbind(Y, X))))] <- svd(cbind(Y, X) %*% t(cbind(Y, X)) / n)$d[1:(min(n, pY+pX))]
		part2 <- sum(log((1-lambda)*part2.slh + lambda))/2
		jointEnt <- part2 + (1 + log(2*pi)) * (pY+pX)/2

		# calculate entropy I
		part2.slh <- rep(0, pY)
		part2.slh[1:(min(dim(Y)))] <- svd(Y %*% t(Y) / n)$d[1:(min(n, pY))]
		part2 <- sum(log((1-lambda)*part2.slh + lambda))/2
		EntI <- part2 + (1 + log(2*pi)) * pY/2

		# calculate entropy II
		part2.slh <- rep(0, pX)
		part2.slh[1:(min(dim(X)))] <- svd(X %*% t(X) / n)$d[1:(min(n, pX))]
		part2 <- sum(log((1-lambda)*part2.slh + lambda))/2
		EntII <- part2 + (1 + log(2*pi)) * pX/2

		return(EntI + EntII - jointEnt)
	}

	knnMI <- function(Y, X, k, rescale){
		###########################################################
		# calculate the knn entropy (per sample)
		###########################################################

		# determine number of features
		pY <- dim(Y)[2]
		pX <- dim(X)[2]

		# rescale both inputs
		if (rescale){
			Y <- Y / mean(apply(Y, 2, mad))
			X <- X / mean(apply(X, 2, mad))
		}

		# determine marginal and joint distance matrices
		distMatY <- as.matrix(dist(Y, method="minkowski", p=2))
		distMatX <- as.matrix(dist(X, method="minkowski", p=2))
		rm(Y, X)
		gc(verbose=FALSE)	
		distMatJoint <- distMatX
		distMatJoint[(distMatY - distMatX > 0)] <- distMatY[(distMatY - distMatX > 0)]

		# calculate joint entropy
		rankDistMat <- t(apply(distMatJoint, 1, rank, ties.method="random"))
		knns <- as.numeric(apply(rankDistMat, 1, function(R, k){ which(R == (k+1)) }, k))
		knnDists <- distMatJoint[cbind(c(1:length(knns)), knns)]
		# jointEnt <- -log(k) + log(dim(distMatY)[2]-1) - lgamma((pY + pX)/2+1) + 0.5*(pY + pX)*log(pi) + (pY + pX)*log(knnDists)
		rm(distMatJoint, rankDistMat, knns)
		gc(verbose=FALSE)	

		# calculate marginal entropy for Y
		kY <- numeric()
		knnDistsY <- numeric()
		for (i in 1:dim(distMatY)[2]){
			kY <- c(kY, length(which(knnDists[i] >= distMatY[i,-i])))
			knnDistsY <- c(knnDistsY, sort((distMatY[i,-i]))[kY[i]])
		}
		# EntY <- -digamma(kY) + log(dim(distMatY)[2]-1) - lgamma(pY/2+1) + 0.5*pY*log(pi) + pY*log(knnDistsY)
		rm(distMatY)
		gc(verbose=FALSE)	

		# calculate marginal entropy for X
		kX <- numeric()
		knnDistsX <- numeric()
		for (i in 1:dim(distMatX)[2]){
			kX <- c(kX, length(which(knnDists[i] >= distMatX[i,-i])))
			knnDistsX <- c(knnDistsX, sort((distMatX[i,-i]))[kX[i]])
		}
		# EntX <- -digamma(kX) + log(dim(distMatX)[2]-1) - lgamma(pX/2+1) + 0.5*pX*log(pi) + pX*log(knnDistsX)

		# return mutual information
		return(digamma(k) + digamma(dim(distMatX)[2]) - mean(digamma(kX) + digamma(kY)))
	}	


	# input checks
	if (as.character(class(Y)) != "matrix"){ stop("Input (Y) is of wrong class.") }
	if (sum(is.na(Y)) != 0){ stop("Matrix Y contains missings.") }
	if (as.character(class(X)) != "matrix"){ stop("Input (X) is of wrong class.") }
	if (sum(is.na(X)) != 0){ stop("Matrix X contains missings.") }
	if (dim(X)[1] != dim(Y)[1]){ stop("Dimension mismatch between matrices X and Y.") }
	if (!(method %in% c("normal", "knn"))){ stop("method parameter ill-specified.") }
	if (sum(c(as.character(class(k)) != "numeric", as.character(class(k)) != "integer")) == 0){ stop("Input (k) is of wrong class.") }
	if (k < 1){ stop("k smaller than 1.") }
	if (as.character(class(center)) != "logical"){ stop("center of wrong class.") }
	if (as.character(class(rescale)) != "logical"){ stop("center of wrong class.") }

	# calculate entropy
	if (method=="normal"){ return(normMI(Y, X, center)) }
	if (method=="knn"){ return(knnMI(Y, X, k, rescale)) }
}

entropyTest <- function(Y, id, nPerm=1000, method="normal", k0=1, k1=1, center=TRUE, lowCiThres=0.10, ncpus=1, verbose=FALSE){
	###########################################################
	# two-sample test on equal entropy
	###########################################################

	# input checks
	if (as.character(class(Y)) != "matrix"){ stop("Input (Y) is of wrong class.") }
	if (sum(c(as.character(class(id)) != "numeric", as.character(class(id)) != "integer")) == 0){ stop("Input (id) is of wrong class.") }
	if (any(sort(unique(id)) != c(0,1))){ stop("Input (id) has wrong entries.") }
	if (length(id) != nrow(Y)){ stop("Dimensions of Y and id do no match.") }
	if (sum(is.na(Y)) != 0){ stop("Matrix Y contains missings.") }
	if (!(method %in% c("normal", "knn"))){ stop("method parameter ill-specified.") }
	if (sum(c(as.character(class(k0)) != "numeric", as.character(class(k0)) != "integer")) == 0){ stop("Input (k0) is of wrong class.") }
	if (sum(c(as.character(class(k1)) != "numeric", as.character(class(k1)) != "integer")) == 0){ stop("Input (k1) is of wrong class.") }
	if (k0 < 1){ stop("k0 smaller than 1.") }
	if (k1 < 1){ stop("k1 smaller than 1.") }
	if (as.character(class(center)) != "logical"){ stop("center of wrong class.") }

	# calculate test statistic
	entObs <- hdEntropy(Y[id==1,], method=method, k=k1, center=TRUE, indKnn=FALSE) - hdEntropy(Y[id==0,], method=method, k=k0, center=TRUE, indKnn=FALSE)	

	# testing through resampling
	steps <- sort(unique(c(0,25,50,100,150,200, seq(from=250,to=2750,by=250), seq(from=3000,to=10000,by=500), seq(from=11000,to=50000,by=1000), nPerm)))
	steps <- steps[steps <= nPerm]
	nullDist <- numeric()
	if (ncpus > 1){
		sfInit(parallel=TRUE, cpus=ncpus)
		sfLibrary(sigaR, verbose=FALSE)
		sfLibrary(corpcor, verbose=FALSE)
	}
	for(j in 1:(length(steps)-1)){
		if (verbose){ cat(paste(steps[j]," of ", steps[length(steps)]," permutations done, and counting...", sep=""), "\n") }	
		if (ncpus == 1){ nullDistPart <- sapply(c((steps[j]+1):steps[j+1]), function(i, Y, id, method, k1, k0, center){
					id <- id[sample(1:length(id), length(id))];
					ENT1 <- hdEntropy(Y[id==1,], method=method, k=k1, center=center, indKnn=FALSE);
					ENT0 <- hdEntropy(Y[id==0,], method=method, k=k0, center=center, indKnn=FALSE);
					return(ENT1-ENT0) }, Y=Y, id=id, method=method, k0=k0, k1=k1, center=center) }
		if (ncpus > 1){	nullDistPart <- sfSapply(c((steps[j]+1):steps[j+1]), function(i, Y, id, method, k1, k0, center){
					id <- id[sample(1:length(id), length(id))];
					ENT1 <- hdEntropy(Y[id==1,], method=method, k=k1, center=center, indKnn=FALSE);
					ENT0 <- hdEntropy(Y[id==0,], method=method, k=k0, center=center, indKnn=FALSE);
					return(ENT1-ENT0) }, Y=Y, id=id, method=method, k0=k0, k1=k1, center=center) }
		nullDist <- c(nullDist, nullDistPart)
		rm(nullDistPart)
		gc()

		# compute pval bound and delete row from data set and permutation set when 0.001 < lower bound
		pVal <- sum(nullDist >= as.numeric(entObs)) / steps[j+1]
		pBound <- pVal - sqrt(pVal*(1-pVal) / steps[j+1]) * 3.09
		significanceUnlikely <- (pBound > lowCiThres)

		# if probability of becoming significant is small, stop resampling.
		if (significanceUnlikely){ pVal <- 1; break }

		# status report
		if (verbose){ cat(paste(steps[j+1], "of", steps[length(steps)], " permutations done", sep=" "), "\n") }
    	}
	if (ncpus > 1){ sfStop() }

	# return output
	return(new("entTest", statistic=entObs, p.value=pVal, null.dist=nullDist, nperm=nPerm, remark=if(significanceUnlikely){ "resampling terminated prematurely: unlikely significance" } else { "none" }))
}

mutInfTest <- function(Y, X, nPerm=1000, method="normal", k=1, center=TRUE, rescale=TRUE, lowCiThres=0.10, ncpus=1, verbose=FALSE){
	###########################################################
	# one-sample test on zero mutual information
	###########################################################

	# input checks
	if (as.character(class(Y)) != "matrix"){ stop("Input (Y) is of wrong class.") }
	if (sum(is.na(Y)) != 0){ stop("Matrix Y contains missings.") }
	if (as.character(class(X)) != "matrix"){ stop("Input (X) is of wrong class.") }
	if (sum(is.na(X)) != 0){ stop("Matrix X contains missings.") }
	if (dim(X)[1] != dim(Y)[1]){ stop("Dimension mismatch between matrices X and Y.") }
	if (sum(c(as.character(class(nPerm)) != "numeric", as.character(class(nPerm)) != "integer")) == 0){ stop("nPerm is wrongly specified.") }
	if (!(method %in% c("normal", "knn"))){ stop("method parameter ill-specified.") }
	if (sum(c(as.character(class(k)) != "numeric", as.character(class(k)) != "integer")) == 0){ stop("Input (k) is of wrong class.") }
	if (k < 1){ stop("k smaller than 1.") }
	if (as.character(class(center)) != "logical"){ stop("center of wrong class.") }
	if (as.character(class(rescale)) != "logical"){ stop("center of wrong class.") }

	# calculate test statistic
	miObs <- hdMI(Y, X, method=method, k=k, center=center, rescale=rescale)

	# testing through resampling
	steps <- sort(unique(c(0,25,50,100,150,200, seq(from=250,to=2750,by=250), seq(from=3000,to=10000,by=500), seq(from=11000,to=50000,by=1000), nPerm)))
	steps <- steps[steps <= nPerm]
	nullDist <- numeric()
	if (ncpus > 1){
		sfInit(parallel=TRUE, cpus=ncpus)
		sfLibrary(sigaR, verbose=FALSE)
		sfLibrary(corpcor, verbose=FALSE)
	}
	for(j in 1:(length(steps)-1)){
		if (verbose){ cat(paste(steps[j]," of ", steps[length(steps)]," permutations done, and counting...", sep=""), "\n") }	
		if (ncpus == 1){ nullDistPart <- sapply(c((steps[j]+1):steps[j+1]), function(i, Y, Z, method, k, center, rescale){
					shuffle <- sample(1:nrow(Y), nrow(Y), replace=FALSE);
					MI <- hdMI(Y[shuffle,], Z, method=method, k=k, center=center, rescale=rescale); 
					return(MI)}, Y=Y, Z=X, method=method, k=k, center=center, rescale=rescale) }
		if (ncpus > 1){ nullDistPart <- sfSapply(c((steps[j]+1):steps[j+1]), function(i, Y, Z, method, k, center, rescale){
					shuffle <- sample(1:nrow(Y), nrow(Y), replace=FALSE);
					MI <- hdMI(Y[shuffle,], Z, method=method, k=k, center=center, rescale=rescale);
					return(MI) }, Y=Y, Z=X, method=method, k=k, center=center, rescale=rescale) }
		nullDist <- c(nullDist, nullDistPart)
		rm(nullDistPart)
		gc()

		# compute pval bound and delete row from data set and permutation set when 0.001 < lower bound
		pVal <- sum(nullDist >= as.numeric(miObs)) / steps[j+1]
		pBound <- pVal - sqrt(pVal*(1-pVal) / steps[j+1]) * 3.09
		significanceUnlikely <- (pBound > lowCiThres)

		# if probability of becoming significant is small, stop resampling.
		if (significanceUnlikely){ pVal <- 1; break }
    	}
	if (ncpus > 1){ sfStop() }

	# return output
	return(new("miTest", statistic=miObs, p.value=pVal, null.dist=nullDist, nperm=nPerm, remark=if(significanceUnlikely){ "resampling terminated prematurely: unlikely significance" } else { "none" }))
}
