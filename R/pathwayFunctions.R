pathway2sample <- function(Y, X, id, lambda1=1, lambdaF=1, method="FL", constr=TRUE, startCis=numeric(), startTrans1=matrix(), startTrans2=matrix(), epsilon=0, verbose=FALSE){
	############################################################
	# fit simultaneous equations model via penalized OLS
	############################################################

	softThresholding <- function(x, threshold){
		############################################################
		# apply soft-thresholding to x with threshold 
		############################################################
		ids1 <- which(x > threshold)
		ids2 <- which(abs(x) <= threshold)
		ids3 <- which(x < -threshold)
		if (length(ids1) > 0){ x[ids1] <- x[ids1] - threshold[ids1] }
		if (length(ids2) > 0){ x[ids2] <- 0 }
		if (length(ids3) > 0){ x[ids3] <- x[ids3] + threshold[ids3] }
		return(x)
	}

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(Y)) != "matrix"){ stop("Input (Y) is of wrong class.") }
	if (sum(is.na(Y)) != 0){ stop("Matrix Y contains missings.") }
	if (as.character(class(X)) != "matrix"){ stop("Input (X) is of wrong class.") }
	if (sum(is.na(X)) != 0){ stop("Matrix X contains missings.") }
	if (dim(X)[1] != dim(Y)[1]){ stop("Dimension mismatch in columns between matrices X and Y.") }
	if (dim(X)[2] != dim(Y)[2]){ stop("Dimension mismatch in rows between matrices X and Y.") }
	if (sum(c(as.character(class(id)) != "numeric", as.character(class(id)) != "integer")) == 0){ stop("Input (id) is of wrong class.") }
	if (any(sort(unique(id)) != c(0,1))){ stop("Input (id) has wrong entries.") }
	if (length(id) != nrow(Y)){ stop("Number of rows of Y and id do no match.") }
	if (!(as.character(class(lambda1)) == "matrix" | as.character(class(lambda1)) == "numeric" | as.character(class(lambda1)) == "integer")){ stop("Input (lambda1) is of wrong class.") }
	if (as.character(class(lambda1)) == "numeric" | as.character(class(lambda1)) == "integer"){ if (length(lambda1) != 1){ stop("lambda1 should have length equal to one.") } }
	if (as.character(class(lambda1)) == "matrix"){ if (nrow(lambda1) != ncol(Y) | ncol(lambda1) != ncol(Y)){ stop("number of rows and columns of lambda1 should  equal the number of columns of Y.") } }
	if (!(as.character(class(lambdaF)) == "matrix" | as.character(class(lambdaF)) == "numeric" | as.character(class(lambdaF)) == "integer")){ stop("Input (lambdaF) is of wrong class.") }
	if (as.character(class(lambdaF)) == "numeric" | as.character(class(lambdaF)) == "integer"){ if (length(lambdaF) != 1){ stop("lambdaF should have length equal to one.") } }
	if (as.character(class(lambdaF)) == "matrix"){ if (nrow(lambdaF) != ncol(Y) | ncol(lambdaF) != ncol(Y)){ stop("number of rows and columns of lambdaF should  equal the number of columns of Y.") } }
	if (as.character(class(method)) != "character"){ stop("Input (method) is of wrong class.") }
	if (!(method %in% c("FL", "FLs"))){ stop("method parameter ill-specified.") }
	if (as.character(class(constr)) != "logical"){ stop("constr of wrong class.") }
	if (sum(c(as.character(class(epsilon)) != "numeric", as.character(class(epsilon)) != "integer")) == 0){ stop("Input (epsilon) is of wrong class.") }
	if (length(epsilon) != 1){ stop("Input (epsilon) of wrong length.") }
	if (epsilon < 0){ stop("Input (epsilon) should be positive.") }
	if (as.character(class(verbose)) != "logical"){ stop("verbose of wrong class.") }
	if (length(startCis) == 0 & ncol(startTrans1) == 1 & ncol(startTrans1) == 1 & ncol(startTrans2) == 1 & ncol(startTrans2) == 1){
		startValues <- FALSE
	} else {
		if (length(startCis) != ncol(Y)){ stop("startCis of wrong length, does not match number of columns of X and Y.") } 
		if (ncol(startTrans1) != ncol(Y)){ stop("startTrans1 of wrong dimensions, number of columns does not match that of Y and X.") } 
		if (nrow(startTrans1) != ncol(Y)){ stop("startTrans1 of wrong dimensions, number of rows does not match the number of columns of Y and X.") } 
		if (ncol(startTrans2) != ncol(Y)){ stop("startTrans2 of wrong dimensions, number of columns does not match that of Y and X.") } 
		if (nrow(startTrans2) != ncol(Y)){ stop("startTrans2 of wrong dimensions, number of rows does not match the number of columns of Y and X.") } 
		startValues <- TRUE
	}

	# define parameters
	transHatG1 <- diag(1, ncol(Y))
	transHatG2 <- diag(1, ncol(Y))
	SigmaHat <- numeric(length=ncol(Y))
	betaHat <- numeric(length=ncol(Y))

	# define objects for analysis
	Z <- rbind(cbind(Y[id==0,], matrix(0, ncol=ncol(Y), nrow=sum(id==0))), cbind(matrix(0, ncol=ncol(Y), nrow=sum(id==1)), Y[id==1,]))
	Dtilde <- rbind(cbind(diag(rep(1, ncol(Y))), diag(rep(-1, ncol(Y)))), cbind(diag(rep(1, ncol(Y))), diag(rep(1, ncol(Y)))))
	Y <- rbind(Y[id==0,], Y[id==1,])
	X <- rbind(X[id==0,], X[id==1,])
	id <- c(id[id==0], id[id==1])
	if (method == "FL"){ if (constr){ posVec <- c(T, rep(F, ncol(Y)-1)) } else { posVec <- rep(F, ncol(Y)) } }
	if (method == "FLs"){ if (constr){ posVec <- c(T, rep(F, 2*(ncol(Y)-1))) } else { posVec <- rep(F, 2*ncol(Y)-1) } }
	lambda1vec <- rep(lambda1, ncol(Y)-1)
	lambdaFvec <- c(0, rep(lambdaF, ncol(Y)-1))

	if (epsilon > 0 & method=="FL"){
		Y <- rbind(Y, matrix(0, ncol=ncol(Y), nrow=ncol(Z)))
		X <- rbind(X, matrix(0, ncol=ncol(X), nrow=ncol(Z)))
		Z <- rbind(Z, epsilon * diag(ncol(Z)))
	}

	# penalized estimation gene-by-gene
	for (j in 1:ncol(Y)){
		# in case regularization varies per parameter
		if (as.character(class(lambda1)) == "matrix"){ lambda1vec <- lambda1[j, -j] }
		if (as.character(class(lambdaF)) == "matrix"){ lambdaFvec <- c(0, lambdaF[j, -j]) }

		if (method=="FL"){
			# define objects for analysis
			Z1 <- (Z[,-c(j, j+ncol(Y))] %*% t(Dtilde[-c(j, j+ncol(Y)), -c(j, j+ncol(Y))])[,c(1:(ncol(Y)-1))])/2
			Z2 <- (Z[,-c(j, j+ncol(Y))] %*% t(Dtilde[-c(j, j+ncol(Y)), -c(j, j+ncol(Y))])[,-c(1:(ncol(Y)-1))])/2
			P <- Z2 %*% solve(t(Z2) %*% Z2) %*% t(Z2)

			# estimation with fusion penalty
		        # penFit <- penalized((diag(rep(1, nrow(P))) - P) %*% Y[,j], cbind(0* X[,j], ((diag(rep(1, nrow(P))) - P)) %*% Z1), unpenalized=~0, lambda1=lambdaFvec, positive=posVec, trace=verbose)
			if (startValues){
		        	penFit <- penalized::penalized((diag(rep(1, nrow(P))) - P) %*% Y[,j], cbind(X[,j], ((diag(rep(1, nrow(P))) - P)) %*% Z1), unpenalized=~0, lambda1=lambdaFvec, positive=posVec, startbeta=c(startCis[j], startTrans1[j, -j] - startTrans2[j, -j]), trace=verbose)
			} else {
		        	penFit <- penalized::penalized((diag(rep(1, nrow(P))) - P) %*% Y[,j], cbind(X[,j], ((diag(rep(1, nrow(P))) - P)) %*% Z1), unpenalized=~0, lambda1=lambdaFvec, positive=posVec, trace=verbose)
			}

			# fabric and format fused lasso estimates
			xi <- coefficients(penFit, "all")[-1]
			R <- matrix(c(1, rep(0, ncol(Z2))), nrow=1)
			QPpars <- solve.QP(R %*% ginv(t(cbind(X[,j], Z2)) %*% cbind(X[,j], Z2)) %*% t(R),  -R %*% ginv(cbind(X[,j], Z2)) %*% Y[,j, drop=FALSE], diag(1))$solution
			zeta <- ginv(t(cbind(X[,j], Z2)) %*% cbind(X[,j], Z2)) %*% t(cbind(X[,j], Z2)) %*% Y[,j] + QPpars * ginv(t(cbind(X[,j], Z2)) %*% cbind(X[,j], Z2)) %*% t(R)
			betaHat[j] <- zeta[1]
			zeta <- zeta[-1]
			estF <- solve(Dtilde[-c(j, j+ncol(Y)), -c(j, j+ncol(Y))]) %*% matrix(c(xi, zeta), ncol=1)
			estFL <- softThresholding(estF, c(lambda1vec, lambda1vec))
			transHatG1[j, -j] <- -estFL[c(1:(ncol(Y)-1))]
			transHatG2[j, -j] <- -estFL[-c(1:(ncol(Y)-1))]
		}
		if (method=="FLs"){
			# penalized estimation
			if (startValues){
				penFit <- penalized::penalized(Y[,j], cbind(X[,j], (Z %*% t(Dtilde))[,-c(j, j+ncol(Y))]), unpenalized=~0, lambda1=c(lambdaFvec, lambda1vec), positive=posVec, startbeta=c(startCis[j], startTrans1[j, -j] - startTrans2[j, -j]), trace=verbose)
			} else {
				penFit <- penalized::penalized(Y[,j], cbind(X[,j], (Z %*% t(Dtilde))[,-c(j, j+ncol(Y))]), unpenalized=~0, lambda1=c(lambdaFvec, lambda1vec), positive=posVec, trace=verbose)
			}
			betaHat[j] <- coefficients(penFit, "all")[1]
			xi <- coefficients(penFit, "all")[2:ncol(Y)]
			zeta <- coefficients(penFit, "all")[-c(1:ncol(Y))]
			estR <- solve(Dtilde[-c(j, j+ncol(Y)), -c(j, j+ncol(Y))]) %*% matrix(c(xi, zeta), ncol=1)
			transHatG1[j, -j] <- -estR[c(1:(ncol(Y)-1))]
			transHatG2[j, -j] <- -estR[-c(1:(ncol(Y)-1))]
		}
		SigmaHat[j] <- var(c(Y[id==0, j, drop=FALSE] - Y[id==0, -j] %*% t(transHatG1[j, -j, drop=FALSE]) - betaHat[j] * X[id==0, j, drop=FALSE], Y[id==1, j, drop=FALSE] - Y[id==1, -j] %*% t(transHatG2[j, -j, drop=FALSE]) - betaHat[j] * X[id==1, j, drop=FALSE]))
	}


	return(new("pathwayFit", Cis=betaHat, Trans=matrix(), Trans1=transHatG1, Trans2=transHatG2, Sigma=SigmaHat, lambda1=data.matrix(lambda1), lambdaF=data.matrix(lambdaF), epsilon=epsilon, method=method, constr=constr)) 
}





pathway1sample <- function(Y, X, lambda1=1, constr=TRUE, startCis=numeric(), startTrans=matrix(), verbose=FALSE){
	########################################################################################################################
	# fit simultaneous equations model via penalized OLS
	########################################################################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(Y)) != "matrix"){ stop("Input (Y) is of wrong class.") }
	if (sum(is.na(Y)) != 0){ stop("Matrix Y contains missings.") }
	if (as.character(class(X)) != "matrix"){ stop("Input (X) is of wrong class.") }
	if (sum(is.na(X)) != 0){ stop("Matrix X contains missings.") }
	if (dim(X)[1] != dim(Y)[1]){ stop("Dimension mismatch in columns between matrices X and Y.") }
	if (dim(X)[2] != dim(Y)[2]){ stop("Dimension mismatch in rows between matrices X and Y.") }
	if (!(as.character(class(lambda1)) == "matrix" | as.character(class(lambda1)) == "numeric" | as.character(class(lambda1)) == "integer")){ stop("Input (lambda1) is of wrong class.") }
	if (as.character(class(lambda1)) == "numeric" | as.character(class(lambda1)) == "integer"){ if (length(lambda1) != 1){ stop("lambda1 should have length equal to one.") } }
	if (as.character(class(lambda1)) == "matrix"){ if (nrow(lambda1) != ncol(Y) | ncol(lambda1) != ncol(Y)){ stop("number of rows and columns of lambda1 should  equal the number of columns of Y.") } }
	if (as.character(class(constr)) != "logical"){ stop("constr of wrong class.") }
	if (as.character(class(verbose)) != "logical"){ stop("verbose of wrong class.") }
	if (as.character(class(startCis)) != "numeric"){ stop("startCis of wrong class.") } 
	if (as.character(class(startTrans)) != "matrix"){ stop("startTrans of wrong class.") } 
	if (length(startCis) == 0 & ncol(startTrans) == 1 & ncol(startTrans) == 1){
		startValues <- FALSE
	} else {
		if (length(startCis) != ncol(Y)){ stop("startCis of wrong length, does not match number of columns of X and Y.") } 
		if (ncol(startTrans) != ncol(Y)){ stop("startTrans of wrong dimensions, number of columns does not match that of Y and X.") } 
		if (nrow(startTrans) != ncol(Y)){ stop("startTrans of wrong dimensions, number of rows does not match the number of columns of Y and X.") } 
		startValues <- TRUE
	}

	# define parameters
	transHat <- diag(1, ncol(Y))
	SigmaHat <- numeric(length=ncol(Y))
	betaHat <- numeric(length=ncol(Y))

	# estimation
	if (constr){ posVec <- c(T, rep(F, ncol(Y)-1)) } else { posVec <- rep(F, ncol(Y)) }
	lambda1vec <- c(0, rep(lambda1, ncol(Y)-1))
	for (j in 1:ncol(Y)){ 
		# in case regularization varies per parameter
		if (as.character(class(lambda1)) == "matrix"){ lambda1vec <- c(0, lambda1[j, -j]) }
	
		# lasso estimation
		if (startValues){
			penFit <- penalized::penalized(Y[,j], cbind(X[,j,drop=FALSE], Y[,-j]), unpenalized=~0, lambda1=lambda1vec, positive=posVec, startbeta=c(startCis[j], startTrans[j, -j]), trace=verbose)
		} else {
			penFit <- penalized::penalized(Y[,j], cbind(X[,j,drop=FALSE], Y[,-j]), unpenalized=~0, lambda1=lambda1vec, positive=posVec, trace=verbose)
		}

		# extract lasso estimates
		transHat[j, -j] <- -coefficients(penFit, "all")[-1]
		betaHat[j] <- coefficients(penFit, "all")[1]
		SigmaHat[j] <- var(residuals(penFit))
		print(j)
	}
	return(new("pathwayFit", Cis=betaHat, Trans=transHat, Trans1=matrix(), Trans2=matrix(), Sigma=SigmaHat, lambda1=data.matrix(lambda1), lambdaF=matrix(), epsilon=numeric(), method=character(), constr=constr)) 
}


pathwayPlot <- function(pFit, directed=TRUE, tWidth=1, cWidth=1, geWidth=10, cnWidth=10, circleDist=1.5, gNames=NULL, main="", remove=FALSE){
	############################################################
	# plot pathway with circular layout
	############################################################

	# check input
	if (as.character(class(pFit)) != "pathwayFit"){ stop("Input (pFit) is of wrong class.") }
	if (sum(is.na(pFit@Trans)) != 0){ stop("Matrix Trans contains missings.") }
	if (length(pFit@Cis) > 0){
		if (sum(is.na(pFit@Cis)) != 0){ stop("Vector Cis contains missings.") }
		if (nrow(pFit@Trans) != length(pFit@Cis)){ stop("Dimension mismatch between Cis and Trans.") }
		if (ncol(pFit@Trans) != length(pFit@Cis)){ stop("Dimension mismatch between Cis and Trans.") }
	}
	if (as.character(class(directed)) != "logical"){ stop("directed of wrong class.") }
	if (as.character(class(tWidth)) != "numeric"){ stop("Input (tWidth) is of wrong class.") }
	if (length(tWidth) != 1){ stop("Length tWidth must equal one.") }
	if (tWidth <= 0){ stop("tWidth must be positive.") }
	if (as.character(class(cWidth)) != "numeric"){ stop("Input (cWidth) is of wrong class.") }
	if (length(cWidth) != 1){ stop("Length cWidth must equal one.") }
	if (cWidth <= 0){ stop("cWidth must be positive.") }
	if (as.character(class(circleDist)) != "numeric"){ stop("Input (circleDist) is of wrong class.") }
	if (length(circleDist) != 1){ stop("Length circleDist must equal one.") }
	if (circleDist <= 0){ stop("circleDist must be positive.") }
	if (!is.null(gNames)){
		if (as.character(class(gNames)) != "character"){ stop("Input (gNames) is of wrong class.") }
		if (length(gNames) != nrow(pFit@Trans)){ stop("gNames is of wrong length.") }
	}

	diag(pFit@Trans) <- 0
	if (remove){
		idsRemove <- which(rowSums(abs(pFit@Trans)) == 0 & colSums(abs(pFit@Trans)) == 0)
		if (length(idsRemove) > 0){
			pFit@Trans <- pFit@Trans[-idsRemove, -idsRemove]
			if (!is.null(pFit@Cis)){ pFit@Cis <- pFit@Cis[-idsRemove] }
			if (!is.null(gNames)){ gNames <- gNames[-idsRemove] }
		}
	}
	pFit@Trans <- tWidth * pFit@Trans
	if (directed){
		gObj <- graph.adjacency(abs(pFit@Trans), "directed", weighted=TRUE)
	} else {
		gObj <- graph.adjacency(abs(pFit@Trans), "undirected", weighted=TRUE)
	}
	if (length(pFit@Cis) == 0){
		plot.igraph(gObj, layout=layout.circle(gObj), vertex.size=geWidth, vertex.label.font=2,
		vertex.color="lightskyblue2", vertex.frame.color=NA, vertex.label=gNames, vertex.label.cex=0.9,
		vertex.label.color="white", edge.width=E(gObj)$weight, vertex.label.family="sans", main=main)
	}
	if (length(pFit@Cis) > 0){
		nGenes <- nrow(pFit@Trans)
		vertexCoord <- layout.circle(gObj)
		vertexCoord <- rbind(vertexCoord, circleDist * vertexCoord)
		if (directed){
			Trans <- cbind(rbind(pFit@Trans, cWidth * diag(pFit@Cis)), rbind(0 * diag(pFit@Cis), matrix(0, ncol=nGenes, nrow=nGenes)))
			gObj <- graph.adjacency(abs(Trans), "directed", weighted=TRUE)
		} else {
			Trans <- cbind(rbind(pFit@Trans, cWidth * diag(pFit@Cis)), rbind(cWidth * diag(pFit@Cis), matrix(0, ncol=length(pFit@Cis), nrow=length(pFit@Cis))))
			gObj <- graph.adjacency(abs(Trans), "undirected", weighted=TRUE)
		}
		plot.igraph(gObj, layout=vertexCoord, vertex.size=c(rep(geWidth, nGenes), rep(cnWidth, nGenes)), vertex.label.font=2,
		vertex.color=c(rep("lightskyblue2", nGenes), rep("orangered4", nGenes)) , 
		vertex.frame.color=NA, vertex.label=c(rep("", nGenes), gNames), vertex.label.cex=0.75,
		vertex.label.color="white", edge.width=E(gObj)$weight, vertex.label.family="sans", main=main)
		legend("bottomright", c("DNA", "RNA"), pch=20, col=c("orangered4", "lightskyblue2"), cex=1.3, pt.cex=4) 
	}
}

