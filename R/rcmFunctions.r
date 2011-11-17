RCMtest <- function(Y, X, R, testType="I", nBoot=100, lowCiThres=0.10, shrinkType="none", estType="normal", corType="unif",  maxNoIt=100, minSuccDist=0.005, returnNullDist=FALSE, ncpus=1, verbose=FALSE){
	######################################################################################################################
	# function that evaluates significance of parameters of multivariate random coefficients model
	######################################################################################################################

	# input checks
	if (as.character(class(Y)) != "matrix"){ stop("Input (Y) is of wrong class.") }
	if (sum(is.na(Y)) != 0){ stop("Y contains missings.") }
	if (as.character(class(X)) != "matrix"){ stop("Input (X) is of wrong class.") }
	if (sum(is.na(X)) != 0){ stop("X contains missings.") }
	if (as.character(class(R)) != "matrix"){ stop("Input (R) is of wrong class.") }
	if (sum(is.na(R)) != 0){ stop("R contains missings.") }
	if ( !(as.character(class(maxNoIt)) == "numeric" | as.character(class(maxNoIt)) == "integer") ){ stop("Input (maxNoIt) is of wrong class.") }
	if (as.character(class(minSuccDist)) != "numeric"){ stop("Input (minSuccDist) is of wrong class.") }
	if (as.character(class(corType)) != "character"){ stop("Input (corType) is of wrong class.") }
	if (as.character(class(estType)) != "character"){ stop("Input (estType) is of wrong class.") }
	if (as.character(class(shrinkType)) != "character"){ stop("Input (shrinkType) is of wrong class.") }
	if (dim(Y)[1] < 2){ stop("This is a multivariate procedure, provide a multivariate data set.") }
	if (!(testType %in% c("I", "II", "III"))){ stop("testType parameter ill-specified.") }
	if (!(corType %in% c("unif", "ar1"))){ stop("corType parameter ill-specified.") }
	if (!(estType %in% c("normal", "robust"))){ stop("estType parameter ill-specified.") }
	if (!(shrinkType %in% c("none", "opt", "full"))){ stop("shrinkType parameter ill-specified.") }
	if (dim(Y)[1] != dim(X)[1]){ stop("Dimension mismatch between expression and design matrix.") }
	if (dim(X)[2] != dim(R)[2]){ stop("Dimension mismatch between constraint and design matrix.") }
	if (minSuccDist <= 0){ stop("Stopping criterion incorrect.") }
	if (maxNoIt < 1){ stop("Maximum number of iterations smaller than 1.") }
	if (as.character(class(lowCiThres)) != "numeric"){ stop("Input (lowCiThres) is of wrong class.") }
	if ( !(lowCiThres > 0 | lowCiThres <= 1) ){ stop("lowCiThres out of range.") }
	if (as.character(class(verbose)) != "logical"){ stop("verbose of wrong class.") }
	if (as.character(class(returnNullDist)) != "logical"){ stop("verbose of wrong class.") }

	# fit the null model
	if (testType=="I"){
		RCMresultH0 <- RCMestimation(Y, X, R, hypothesis="H0", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
		RCMresultH2 <- RCMestimation(Y, X, R, hypothesis="H2", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
		likRatio <- -2*(RCMresultH0@loglik - RCMresultH2@loglik)
	}
	if (testType=="II"){
		RCMresultH0 <- RCMestimation(Y, X, R, hypothesis="H0", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
		RCMresultH1 <- RCMestimation(Y, X, R, hypothesis="H1", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
		likRatio <- -2*(RCMresultH0@loglik - RCMresultH1@loglik)
	}	
	if (testType=="III"){
		RCMresultH1 <- RCMestimation(Y, X, R, hypothesis="H1", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
		RCMresultH2 <- RCMestimation(Y, X, R, hypothesis="H2", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
		likRatio <- -2*(RCMresultH1@loglik - RCMresultH2@loglik)
	}

	# testing through resampling
	likRatioNull <- numeric()
	steps <- sort(unique(c(0,25,50,100,150,200, seq(from=250,to=2750,by=250), seq(from=3000,to=10000,by=500), seq(from=11000,to=50000,by=1000), nBoot)))
	steps <- steps[steps <= nBoot]
	if (ncpus > 1){
		sfInit(parallel=TRUE, cpus=ncpus)
		sfLibrary(sigaR, verbose=FALSE)
		sfLibrary(quadprog, verbose=FALSE)
		sfLibrary(MASS, verbose=FALSE)
	}
	for(j in 1:(length(steps)-1)){
		if (verbose){ cat(paste(steps[j]," of ", steps[length(steps)]," bootstraps done, and counting...", sep=""), "\n") }
		if (testType=="I"){
			if (ncpus == 1){ likRatioNullPart <- sapply(c((steps[j]+1):steps[j+1]), function(i, RCMresultH0, Z, R, shrinkType, estType, corType, maxNoIt, minSuccDist, verbose){
							Yboot <- RCMrandom(RCMresultH0)
							RCMresultH0temp <- RCMestimation(Yboot, Z, R, hypothesis="H0", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							RCMresultH2temp <- RCMestimation(Yboot, Z, R, hypothesis="H2", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							return(-2*(RCMresultH0temp@loglik - RCMresultH2temp@loglik))
						}, RCMresultH0=RCMresultH0, Z=X, R=R, shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose) }
			if (ncpus > 1){ likRatioNullPart <- sfSapply(c((steps[j]+1):steps[j+1]), function(i, RCMresultH0, Z, R, shrinkType, estType, corType, maxNoIt, minSuccDist, verbose){
							Yboot <- RCMrandom(RCMresultH0)
							RCMresultH0temp <- RCMestimation(Yboot, Z, R, hypothesis="H0", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							RCMresultH2temp <- RCMestimation(Yboot, Z, R, hypothesis="H2", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							return(-2*(RCMresultH0temp@loglik - RCMresultH2temp@loglik))
						}, RCMresultH0=RCMresultH0, Z=X, R=R, shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose) }
		}
		if (testType=="II"){
			if (ncpus == 1){ likRatioNullPart <- sapply(c((steps[j]+1):steps[j+1]), function(i, RCMresultH0, Z, R, shrinkType, estType, corType, maxNoIt, minSuccDist, verbose){
							Yboot <- RCMrandom(RCMresultH0)
							RCMresultH0temp <- RCMestimation(Yboot, Z, R, hypothesis="H0", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							RCMresultH1temp <- RCMestimation(Yboot, Z, R, hypothesis="H1", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							return(-2*(RCMresultH0temp@loglik - RCMresultH1temp@loglik))
						}, RCMresultH0=RCMresultH0, Z=X, R=R, shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose) }
			if (ncpus > 1){ likRatioNullPart <- sfSapply(c((steps[j]+1):steps[j+1]), function(i, RCMresultH0, Z, R, shrinkType, estType, corType, maxNoIt, minSuccDist, verbose){
							Yboot <- RCMrandom(RCMresultH0)
							RCMresultH0temp <- RCMestimation(Yboot, Z, R, hypothesis="H0", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							RCMresultH1temp <- RCMestimation(Yboot, Z, R, hypothesis="H1", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							return(-2*(RCMresultH0temp@loglik - RCMresultH1temp@loglik))
						}, RCMresultH0=RCMresultH0, Z=X, R=R, shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose) }
		}
		if (testType=="III"){
			if (ncpus == 1){ likRatioNullPart <- sapply(c((steps[j]+1):steps[j+1]), function(i, RCMresultH1, Z, R, shrinkType, estType, corType, maxNoIt, minSuccDist, verbose){
							Yboot <- RCMrandom(RCMresultH1)
							RCMresultH1temp <- RCMestimation(Yboot, Z, R, hypothesis="H1", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							RCMresultH2temp <- RCMestimation(Yboot, Z, R, hypothesis="H2", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							return(-2*(RCMresultH1temp@loglik - RCMresultH2temp@loglik))
						}, RCMresultH1=RCMresultH1, Z=X, R=R, shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose) }
			if (ncpus > 1){ likRatioNullPart <- sfSapply(c((steps[j]+1):steps[j+1]), function(i, RCMresultH1, Z, R, shrinkType, estType, corType, maxNoIt, minSuccDist, verbose){
							Yboot <- RCMrandom(RCMresultH1)
							RCMresultH1temp <- RCMestimation(Yboot, Z, R, hypothesis="H1", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							RCMresultH2temp <- RCMestimation(Yboot, Z, R, hypothesis="H2", shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose)
							return(-2*(RCMresultH1temp@loglik - RCMresultH2temp@loglik))
						}, RCMresultH1=RCMresultH1, Z=X, R=R, shrinkType=shrinkType, estType=estType, corType=corType,  maxNoIt=maxNoIt, minSuccDist=minSuccDist, verbose=verbose) }
		}
		likRatioNull <- c(likRatioNull, likRatioNullPart)
		rm(likRatioNullPart)
		gc()

		# compute pval bound and delete row from data set and permutation set when 0.001 < lower bound
		pVal <- sum(likRatioNull >= as.numeric(likRatio)) / steps[j+1]
		pBound <- pVal - sqrt(pVal*(1-pVal) / steps[j+1]) * 3.09
		significanceUnlikely <- (pBound > lowCiThres)

		# if probability of becoming significant is small, stop resampling.
		if (significanceUnlikely){ pVal <- 1; break }

    	}
	if (ncpus > 1){ sfStop() }
	
	# return rcmTest object
	if (testType=="I"){
		return(new("rcmTest", statistic=likRatio, p.value=pVal, betas=RCMresultH2@betas, tau2s=RCMresultH2@tau2s, sigma2s=RCMresultH2@sigma2s, rho=RCMresultH2@rho, av.sigma2s=RCMresultH2@av.sigma2s, shrinkage=RCMresultH2@shrinkage, loglik=RCMresultH2@loglik, nBoot=nBoot, corType=corType, null.dist=if(returnNullDist){ likRatioNull } else {numeric() }, remark=if(significanceUnlikely){ "resampling terminated prematurely: unlikely significance" } else { "none" }))
	}
	if (testType=="II"){
		return(new("rcmTest", statistic=likRatio, p.value=pVal, betas=RCMresultH1@betas, tau2s=RCMresultH1@tau2s, sigma2s=RCMresultH1@sigma2s, rho=RCMresultH1@rho, av.sigma2s=RCMresultH1@av.sigma2s, shrinkage=RCMresultH1@shrinkage, loglik=RCMresultH1@loglik, nBoot=nBoot, corType=corType, null.dist=if(returnNullDist){ likRatioNull } else {numeric() }, remark=if(significanceUnlikely){ "resampling terminated prematurely: unlikely significance" } else { "none" }))
	}
	if (testType=="III"){
		return(new("rcmTest", statistic=likRatio, p.value=pVal, betas=RCMresultH2@betas, tau2s=RCMresultH2@tau2s, sigma2s=RCMresultH2@sigma2s, rho=RCMresultH2@rho, av.sigma2s=RCMresultH2@av.sigma2s, shrinkage=RCMresultH2@shrinkage, loglik=RCMresultH2@loglik, nBoot=nBoot, corType=corType, null.dist=if(returnNullDist){ likRatioNull } else {numeric() }, remark=if(significanceUnlikely){ "resampling terminated prematurely: unlikely significance" } else { "none" }))
	}
}



RCMestimation <- function(Y, X, R, hypothesis="H2", shrinkType="none", estType="normal", corType="unif",  maxNoIt=100, minSuccDist=0.005, verbose=FALSE){
	######################################################################################################################
	# function that performs constraint estimation in the multivariate random coefficients model
	######################################################################################################################

	RCMmlH2 <- function(Y, X, R, maxNoIt, minSuccDist, shrinkType="none", estType="normal", corType="unif"){
		######################################################################################################################
		# function that performs constraint estimation in the multivariate random coefficients model under H2
		######################################################################################################################

		QPest.RCM <- function(cov.mat, ed.Om2, X, Y.as.vec, R, R.dual, weights, np){
			###########################################################
			# function that estimates the betas-bars
			###########################################################

			mat.calc.1 <- function(ev.cov, ed.Om2, X, weights){
				QPmat1 <- matrix(0, ncol=dim(X)[2], nrow=dim(X)[2])
				for (i in 1:length(ed.Om2$values)){
					QPmat1 <- QPmat1 + (1 / (ev.cov[1] + ed.Om2$values[i])) * (t((matrix(1, nrow=length(ev.cov[-1]), ncol=1) %x% X)) %*% (matrix(ev.cov[-1] %x% ed.Om2$vectors[,i], ncol=1) * matrix(weights, ncol=1)) %*% (t((matrix(ev.cov[-1] %x% ed.Om2$vectors[,i], ncol=1) * matrix(weights, ncol=1))) %*% ((matrix(1, nrow=length(ev.cov[-1]), ncol=1) %x% X))))
				}
				return(QPmat1)
			}

			mat.calc.2 <- function(ev.cov, ed.Om2, X, Y.as.vec, weights){
				QPmat2 <- matrix(0, ncol=1, nrow=dim(X)[2])
				for (i in 1:length(ed.Om2$values)){
					QPmat2 <- QPmat2 + (1 / (ev.cov[1] + ed.Om2$values[i])) * (t((matrix(1, nrow=length(ev.cov[-1]), ncol=1) %x% X)) %*% (matrix(ev.cov[-1] %x% ed.Om2$vectors[,i], ncol=1) * matrix(weights, ncol=1))) %*% (t(matrix(ev.cov[-1] %x% ed.Om2$vectors[,i], ncol=1) * matrix(weights, ncol=1)) %*% Y.as.vec)
				}
				return(QPmat2)
			}

			# calculation of necessary matrices and vectors
			ed.cov <- eigen(cov.mat, symmetric=TRUE)
			ed.cov <- rbind(ed.cov$values, ed.cov$vectors)
			if (np == 1){
				slh.mat1 <- matrix(mean(apply(ed.cov, 2, mat.calc.1, ed.Om2, X, weights)), ncol=dim(X)[2])
				slh.mat2 <- matrix(mean(apply(ed.cov, 2, mat.calc.2, ed.Om2, X, Y.as.vec, weights)), ncol=1)
			}
			if (np > 1){
				slh.mat1 <- matrix(apply(t(apply(ed.cov, 2, mat.calc.1, ed.Om2, X, weights)), 2, mean), ncol=dim(X)[2])
				slh.mat2 <- matrix(apply(t(apply(ed.cov, 2, mat.calc.2, ed.Om2, X, Y.as.vec, weights)), 2, mean), ncol=1)
			}
			Dmat <- R %*% solve(slh.mat1) %*% t(R)
			dvec <- R %*% solve(slh.mat1) %*% slh.mat2

			# maximize lagrangian by means of quadratic programming
			qp.sol <- solve.QP(Dmat, -dvec, R.dual, bvec)
			lambda <- matrix(qp.sol$solution, ncol=1)
			betas.conc <- solve(slh.mat1) %*% slh.mat2 + solve(slh.mat1) %*% t(R) %*% lambda

			return(as.numeric(betas.conc))
		}

		# estimate covariance matrix parameters
		cov.mat.par.ests <- cov.par.estimation(projectY2Xcomp(Y, X), dim(Y)[1], dim(Y)[2], estType, corType)

		# specify initial parameter values from estimates of SURE model with regression coefficients identical for all genes
		betas <- rep(0, dim(X)[2])
		sigma2s <- cov.mat.par.ests$sigma2s
		shrink.par <- shrink.par.calculation(sigma2s, dim(Y)[2], shrinkType)
		sigma2s <- (1 - shrink.par) * sigma2s + shrink.par * mean(sigma2s)
		tau2s <- rep(0, dim(X)[2])
		rho <- cov.mat.par.ests$rho
		if (rho < 0){ warning("The estimate of rho is negative. To ensure the positive definiteness of the covariance matrix, it is set to zero in the remainder." ); rho <- 0 	}

		# put data, design and constraints in right format
		X.circ <- matrix(rep(1, dim(Y)[1]), ncol=1) %x% X
		# X.star <- diag(rep(1, dim(Y)[1])) %x% X 
		Y.as.vec <- matrix(as.numeric(t(Y)), ncol=1)

		# matrices and vectors for the lagrangian (part 1)
		R.dual <- diag(rep(1, dim(R)[1]))
		bvec <- rep(0, dim(R)[1])

		# start of ML estimation by iterative procedure
		for (i1 in 1:maxNoIt){

			# report status
			if (verbose){ cat(paste("ML estimation: iteration ", i1, " of ", maxNoIt, sep=""), "\n") }

			# estimate betas
			cov.mat <- .covMatConstruction(sigma2s, rho, corType)
			tau.mat <- matrix(0, nrow=length(as.numeric(tau2s)), ncol=length(as.numeric(tau2s)))
			diag(tau.mat) <- as.numeric(tau2s)
			Omega.mat.2 <- X %*% tau.mat %*% t(X)
			ed.Om2 <- eigen(Omega.mat.2, symmetric=TRUE)
			betas.old <- betas
			weights <- matrix(1, ncol=1, nrow=length(Y.as.vec))
			betas <- QPest.RCM(cov.mat, ed.Om2, X, Y.as.vec, R, R.dual, weights, dim(X)[2])

			# estimate tau's
			Y.res <- matrix(Y.as.vec - X.circ %*% matrix(betas, ncol=1), nrow=dim(Y)[1], byrow=TRUE)
			Y.res <- Y.res - matrix(apply(Y.res, 1, mean), ncol=dim(Y.res)[2], nrow=dim(Y.res)[1], byrow=FALSE)
			# tau.estimator <- function(gen.res, Xmat, ns){ solve(t(Xmat) %*% Xmat) %*% t(Xmat) %*% matrix(gen.res, ncol=1) %*% t(matrix(gen.res, ncol=1)) %*% Xmat %*% solve(t(Xmat) %*% Xmat) }
			tau.estimator <- function(gen.res, Xmat){ solve(t(Xmat) %*% Xmat) %*% t(Xmat) %*% matrix(gen.res, ncol=1) %*% t(matrix(gen.res, ncol=1)) %*% Xmat %*% solve(t(Xmat) %*% Xmat) }
			tau2s.old <- tau2s
			if (dim(X)[2] == 1){
				# tau2s <- apply(Y.res, 1, tau.estimator, Xmat=X, ns=dim(Y)[2])
				tau2s <- apply(Y.res, 1, tau.estimator, Xmat=X)
				tau2s <- tau2s - sigma2s * solve(t(X) %*% X)
				tau2s <- mean(tau2s)
				if (tau2s < 0){ tau2s <- 0 }
			} else {
				# tau2s <- apply(Y.res, 1, tau.estimator, Xmat=X, ns=dim(Y)[2])
				tau2s <- apply(Y.res, 1, tau.estimator, Xmat=X)
				tau2s.per.gene <- function(st, Xmat, np){ return(diag(matrix(st[-1], ncol=np) - st[1] * solve(t(Xmat) %*% Xmat))) }
				tau2s <- apply(rbind(sigma2s, tau2s), 2, tau2s.per.gene, Xmat=X, np=dim(X)[2])
				tau2s <- apply(tau2s, 1, mean)
				tau2s[(tau2s < 0)] <- 0 
			}

			# robust estimation
			if (estType == "robust"){
				weights <- matrix(as.numeric(t(obs.weights(Y.res, dim(Y)[1], dim(Y)[2], shrink.par))), ncol=1)
				rm(Y.res)

				# specify covariance matrix
				cov.mat <- .covMatConstruction(sigma2s, rho, corType)
				tau.mat <- matrix(0, nrow=length(as.numeric(tau2s)), ncol=length(as.numeric(tau2s)))
				diag(tau.mat) <- as.numeric(tau2s)
				Omega.mat.2 <- X %*% tau.mat %*% t(X)
				ed.Om2 <- eigen(Omega.mat.2, symmetric=TRUE)
				betas <- QPest.RCM(cov.mat, ed.Om2, X, Y.as.vec, R, R.dual, weights, dim(X)[2])

				# estimate tau's
				Y.res <- matrix(Y.as.vec - X.circ %*% matrix(betas, ncol=1), nrow=dim(Y)[1], byrow=TRUE)
				Y.res <- Y.res - matrix(apply(Y.res, 1, mean), ncol=dim(Y.res)[2], nrow=dim(Y.res)[1], byrow=FALSE)
				tau.estimator <- function(gen.res, Xmat, ns){ solve(t(Xmat) %*% Xmat) %*% t(Xmat) %*% matrix(sqrt(gen.res[(ns+1):(2*ns)]) * gen.res[1:ns], ncol=1) %*% t(matrix(sqrt(gen.res[(ns+1):(2*ns)]) * gen.res[1:ns], ncol=1)) %*% Xmat %*% solve(t(Xmat) %*% Xmat) }
				# tau.estimator <- function(gen.res, Xmat){ solve(t(Xmat) %*% Xmat) %*% t(Xmat) %*% matrix(sqrt(gen.res[(ns+1):(2*ns)]) * gen.res[1:ns], ncol=1) %*% t(matrix(sqrt(gen.res[(ns+1):(2*ns)]) * gen.res[1:ns], ncol=1)) %*% Xmat %*% solve(t(Xmat) %*% Xmat) }
				weights <- matrix(weights, nrow=dim(Y)[1], byrow=TRUE)
				betasPerGene <- apply(Y, 1, lm.uc, Xmat=X, estType="robust", mad.times=4)
				betasResiduals <- betasPerGene - betas
				if (dim(X)[2] == 1){
					weights <- obs.weights.betas(betasResiduals, mad.times=5)
					tau2s <- sum(weights*(betasResiduals)^2)/ length(sigma2s) - sum(sigma2s * solve(t(X) %*% X)) / length(sigma2s)
					tau2s[(tau2s < 0)] <- 0 
				} else {
					tau2s <- apply(cbind(Y.res, weights), 1, tau.estimator, Xmat=X, dim(Y)[2])
					# tau2s <- apply(cbind(Y.res, weights), 1, tau.estimator, Xmat=X)
					tau2s.per.gene <- function(st, Xmat, np){ return(diag(matrix(st[-1], ncol=np) - st[1] * solve(t(Xmat) %*% Xmat))) }
					tau2s <- apply(rbind(sigma2s, tau2s), 2, tau2s.per.gene, Xmat=X, np=dim(X)[2])
					tau2s <- apply(tau2s, 1, mean)
					tau2s[(tau2s < 0)] <- 0 
				}
			}
	
			# early stopping if convergence has been reached
			if ( max(abs(c(as.numeric(betas), sqrt(tau2s)) - c(as.numeric(betas.old), sqrt(tau2s.old)))) < minSuccDist){ break }
		}

		estRes <- new("rcmFit", betas=as.numeric(round(betas, digits=5)), tau2s=round(tau2s, digits=10), sigma2s=as.numeric(round(sigma2s, digits=5)), rho=round(rho, digits=5), shrinkage=round(shrink.par, digits=5), av.sigma2s=round(mean(sigma2s), digits=5), loglik=1, corType=corType, X=X)
		estRes@loglik <- .RCMloss(estRes, Y)
		return(estRes)		
	}

	RCMmlH0orH1 <- function(Y, X, R, shrinkType="none", estType, corType, hypothesis){
		######################################################################################################################
		# function that performs estimation in the multivariate SURE model with linear contraints and all genes sharing the same regression coeffients
		######################################################################################################################
	
		# estimate covariance matrix parameters
		cov.mat.par.ests <- cov.par.estimation(projectY2Xcomp(Y, X), dim(Y)[1], dim(Y)[2], estType, corType)

		# specify initial parameter values from estimates of SURE model with regression coefficients identical for all genes
		betas <- rep(0, dim(X)[2])
		sigma2s <- cov.mat.par.ests$sigma2s
		shrink.par <- shrink.par.calculation(sigma2s, dim(Y)[2], shrinkType)
		sigma2s <- (1 - shrink.par) * sigma2s + shrink.par * mean(sigma2s)
		rho <- cov.mat.par.ests$rho

		# put data, design and constraints in right format
		X.circ <- matrix(rep(1, dim(Y)[1]), ncol=1) %x% X
		# X.star <- diag(rep(1, dim(Y)[1])) %x% X 
		Y.as.vec <- matrix(as.numeric(t(Y)), ncol=1)

		# matrices and vectors for the lagrangian (part 1)
		R.star.dual <- diag(rep(1, dim(R)[1]))
		bvec <- rep(0, dim(R)[1])

		if (hypothesis=="H0"){
			# calculate betas
			cov.mat.inv <- solve(.covMatConstruction(sigma2s, rho, corType))
			slh.mat1 <- (matrix(1, ncol=dim(Y)[1], nrow=1) %*% cov.mat.inv %*% matrix(1, nrow=dim(Y)[1], ncol=1)) %x% (t(X) %*% X)
			slh.mat2 <- ((matrix(1, ncol=dim(Y)[1], nrow=1) %*% cov.mat.inv) %x% t(X)) %*% Y.as.vec
			beta.unconstr <- solve(slh.mat1) %*% slh.mat2
			beta.eqconstr <- beta.unconstr - solve(slh.mat1) %*% (t(R) %*% solve(R %*% solve(slh.mat1) %*% t(R)) %*% R %*% beta.unconstr)
			betas <- matrix(beta.eqconstr, ncol=dim(X)[2], byrow=TRUE)
		} 

		if (hypothesis=="H1"){
			# matrices and vectors for the lagrangian (part 2)
			cov.mat.inv <- solve(.covMatConstruction(sigma2s, rho, corType))
			slh.mat <- solve((matrix(1, ncol=dim(Y)[1], nrow=1) %*% cov.mat.inv %*% matrix(1, nrow=dim(Y)[1], ncol=1)) %x% (t(X) %*% X))
			Dmat <- R %*% slh.mat %*% t(R)
			dvec <- -R %*% slh.mat %*% ((matrix(1, ncol=dim(Y)[1], nrow=1) %*% cov.mat.inv) %x% t(X)) %*% Y.as.vec
	
			# maximize lagrangian by means of quadratic programming
			qp.sol <- solve.QP(Dmat, dvec, R.star.dual, bvec)
			lambda <- matrix(qp.sol$solution, ncol=1)

			# calculate inequality constrained parameter estimates and the residual error
			betas <- slh.mat %*% ((matrix(1, ncol=dim(Y)[1], nrow=1) %*% cov.mat.inv) %x% t(X)) %*% Y.as.vec + slh.mat %*% t(R) %*% lambda
			betas <- matrix(betas, ncol=dim(X)[2], byrow=TRUE)
		}

		# robust estimation
		if (estType=="robust"){
			res.mat.star <- matrix(Y.as.vec - X.circ %*% as.numeric(t(betas)), nrow=dim(Y)[1], byrow=TRUE)
			weights <- matrix(as.numeric(t(obs.weights(res.mat.star, dim(Y)[1], dim(Y)[2], shrink.par))), ncol=1)
			X.circ <- matrix(sqrt(weights), ncol=dim(X)[2], nrow=dim(Y)[2]*dim(Y)[1]) * X.circ

			# eigendecompostion and evaluation of quadratic loss term in dual function
			ed.cov <- eigen(cov.mat.inv, symmetric=TRUE)
			ed.cov <- rbind(ed.cov$values, ed.cov$vectors)
			WY <- sqrt(weights) * Y.as.vec

			if (dim(X)[2] == 1){
				slh.mat1 <- matrix(sum(apply(ed.cov, 2, slh.mat.weighted.one, WX.circ=X.circ, ns=dim(Y)[2])), ncol=dim(X)[2])
				if (hypothesis=="H1"){ slh.mat1 <- solve(slh.mat1) }
				slh.mat2 <- matrix(sum(apply(ed.cov, 2, slh.mat.weighted.two, WX.circ=X.circ, WY=WY, ns=dim(Y)[2])), nrow=dim(X)[2])
			} else {
				slh.mat1 <- matrix(rowSums(apply(ed.cov, 2, slh.mat.weighted.one, WX.circ=X.circ, ns=dim(Y)[2])), ncol=dim(X)[2])
				if (hypothesis=="H1"){ slh.mat1 <- solve(slh.mat1) }
				slh.mat2 <- matrix(rowSums(apply(ed.cov, 2, slh.mat.weighted.two, WX.circ=X.circ, WY=WY, ns=dim(Y)[2])), nrow=dim(X)[2])
			}
			rm(ed.cov)

			if (hypothesis=="H0"){ 
				beta.unconstr <- solve(slh.mat1) %*% slh.mat2
				beta.eqconstr <- beta.unconstr - solve(slh.mat1) %*% (t(R) %*% solve(R %*% solve(slh.mat1) %*% t(R)) %*% R %*% beta.unconstr)
				betas <- matrix(beta.eqconstr, ncol=dim(X)[2], byrow=TRUE)			
			}

			if (hypothesis=="H1"){ 
				# matrices and vectors for the lagrangian (part 2)
				Dmat <- R %*% slh.mat1 %*% t(R)
				dvec <- -R %*% slh.mat1 %*% slh.mat2
	
				# maximize lagrangian by means of quadratic programming
				qp.sol <- solve.QP(Dmat, dvec, R.star.dual, bvec)
				lambda <- matrix(qp.sol$solution, ncol=1)

				# calculate inequality constrained parameter estimates and the residual error
				betas <- slh.mat1 %*% slh.mat2 + slh.mat1 %*% t(R) %*% lambda
				betas <- matrix(betas, ncol=dim(X)[2], byrow=TRUE)
			}
		}
		estRes <- new("rcmFit", betas=as.numeric(round(betas, digits=5)), tau2s=rep(0, length(betas)), sigma2s=as.numeric(round(sigma2s, digits=5)), rho=round(rho, digits=5), shrinkage=round(shrink.par, digits=5), av.sigma2s=round(mean(sigma2s), digits=5), loglik=1, corType=corType, X=X)
		estRes@loglik <- .RCMloss(estRes, Y)
		return(estRes)
	}

	projectY2Xcomp <- function(Y, X){
			###########################################################
			# function that projects Y onto the orthogonal complement of X
			###########################################################
			Uorth <- Null(svd(X)$u)
			Porth <- Uorth %*% solve(t(Uorth) %*% Uorth) %*% t(Uorth)
			return(t(Porth %*% t(Y)))
	}

	slh.mat.weighted.one <- function(ev.cov, WX.circ, ns){
		###########################################################
		# function that calculates the loss of the SURE model
		###########################################################
		junk <- 0
		for (k in 1:ns){
			ev.unit <- matrix(0, ncol=1, nrow=ns)
			ev.unit[k,1] <- 1
			junk <- junk + ev.cov[1] * t(WX.circ) %*% (ev.cov[-1] %x% ev.unit) %*% t(t(WX.circ) %*% (ev.cov[-1] %x% ev.unit))
		}
		return(junk)
	}

	slh.mat.weighted.two <- function(ev.cov, WX.circ, WY, ns){
		###########################################################
		# function that calculates the loss of the SURE model
		###########################################################
		junk <- 0
		for (k in 1:ns){
			ev.unit <- matrix(0, ncol=1, nrow=ns)
			ev.unit[k,1] <- 1
			junk <- junk + ev.cov[1] * t(WX.circ) %*% (ev.cov[-1] %x% ev.unit) %*% (t(ev.cov[-1] %x% ev.unit) %*% WY)
		}
		return(junk)
	}

	shrink.par.calculation <- function(sigma2s, ns, shrinkType){
		###########################################################
		# function that calculates the shrinkage parameter
		###########################################################
		if (shrinkType=="none"){
			shrink.par <- 0
		}	
		if (shrinkType=="full"){
			shrink.par <- 1
		}
		if (shrinkType=="opt"){
			teller <- (2 * length(sigma2s) - 3) * sum(sigma2s^2) / ((ns-1) * length(sigma2s))
			noemer <- sum((sigma2s - mean(sigma2s))^2) + (2 * length(sigma2s) - 4) * sum(sigma2s^2) / ((ns-1) * length(sigma2s))
			shrink.par <- max(0, min(1, teller/noemer))
		}
		return(shrink.par)
	}

	cov.par.estimation <- function(res.mat, ng, ns, estType="robust", corType="unif"){
		###########################################################
		# function that estimates the covariance matrix parameters
		###########################################################
		if (corType=="unif"){
			if (estType=="normal"){
				res.mat <- res.mat - matrix(apply(res.mat, 1, mean), ncol=dim(res.mat)[2], nrow=dim(res.mat)[1], byrow=FALSE)
				res.inprods.mat <- res.mat %*% t(res.mat)
				rho <- (sum(res.inprods.mat[upper.tri(res.inprods.mat)]) / (ng * (ng - 1) / 2) / (sum(diag(res.inprods.mat)) / ng))
				sigma2s <- diag(res.inprods.mat)/(ns)
			} 
			if (estType=="robust"){
				sigma2s <- apply(res.mat, 1, mad)^2
				cov.sum <- function(j, ng, X){ sum(apply(X[(j+1):ng, , drop=FALSE], 1, function(x, y){ (mad(x + y))^2 - (mad(x - y))^2 }, X[j,])) }
				cov <- sum(sapply(c(1:(ng-1)), cov.sum, ng=ng, res.mat))
				rho <- (2 * cov / (4 * ng * (ng - 1))) / mean(sigma2s)
			} 
		}
		if (corType=="ar1"){
			if (estType=="normal"){
				res.mat <- res.mat - matrix(apply(res.mat, 1, mean), ncol=dim(res.mat)[2], nrow=dim(res.mat)[1], byrow=FALSE)
				res.inprods.mat <- res.mat %*% t(res.mat)
				sigma2s <- diag(res.inprods.mat)/(ns)
				noemer <- (sum(diag(res.inprods.mat)) / ng)
				res.inprods.mat <- res.inprods.mat[-ng, -1]
				rho <- (sum(diag(res.inprods.mat)) / (ng - 1)  / noemer)
			} 
			if (estType=="robust"){
				sigma2s <- apply(res.mat, 1, mad)^2
				cov.sum <- function(j, ng, X){ sum(apply(X[j+1, , drop=FALSE], 1, function(x, y){ (mad(x + y))^2 - (mad(x - y))^2 }, X[j,])) }
				cov <- sum(sapply(c(1:(ng-1)), cov.sum, ng=ng, res.mat))
				rho <- (cov / (4 * (ng - 1))) / mean(sigma2s)
			} 
		}
		cov.pars <- list()
		cov.pars$rho <- rho
		cov.pars$sigma2s <- sigma2s
		return(cov.pars)
	}

	obs.weights.betas <- function(res.mat, mad.times=5){
		###########################################################
		# function that calculates observation weights
		# in accordance with Tukey's biweight function
		###########################################################
		res.mat <- res.mat - median(res.mat)
		mads <- mad(res.mat)
		res.mat <- res.mat * 1/(mad.times*mads)
		res.mat[abs(res.mat) > 1] <- 1
		weights <- (1 - res.mat^2)^2
		return(weights)
	}

	obs.weights <- function(res.mat, ng, ns, shrink.par=0, mad.times=5){
		###########################################################
		# function that calculates observation weights
		# in accordance with Tukey's biweight function
		###########################################################
		res.mat <- res.mat - apply(res.mat, 1, median)
		mads <- apply(res.mat, 1, mad)
		mads <- (1-shrink.par) * mads + shrink.par * mean(mads)
		res.mat <- t(apply(cbind(res.mat, mads), 1, function(x, mad.times){ spread <- x[length(x)]; x <- x[-length(x)]; x[abs(x) > mad.times*spread] <- mad.times*spread; return(x) }, mad.times))
		res.mat <- res.mat * matrix(1/(mad.times*mads), ncol=ns, nrow=ng, byrow=FALSE)
		weights <- (1 - res.mat^2)^2
		return(weights)
	}

	lm.uc <- function(Y, Xmat, estType, mad.times=3){
		###########################################################
		# function that fits the univariate regression model
		###########################################################
		# specify covariance matrix and calculate (or specify) its inverse
		sigma <- 1
		cov.mat.inv <- 1/sigma*diag(rep(1, length(Y)))

		# calculate parameter estimates
		betas <- solve(t(Xmat) %*% cov.mat.inv %*% Xmat) %*% t(Xmat) %*% cov.mat.inv %*% Y
		if (estType == "robust"){
			residuals <- Y - X %*% betas
			res.mad <- mad(residuals)
			residuals[residuals > mad.times*res.mad] <- mad.times*res.mad
			weights <- (1 - (residuals/(mad.times*mad(residuals)))^2)^2
			betas <- solve(t(Xmat) %*% diag(sqrt(as.numeric(weights))) %*% cov.mat.inv %*% diag(sqrt(as.numeric(weights))) %*% Xmat) %*% t(Xmat) %*% cov.mat.inv %*% Y
		}

		# return unconstrained parameter estimates
		return(betas)
	}

	# input checks
	if (as.character(class(Y)) != "matrix"){ stop("Input (Y) is of wrong class.") }
	if (sum(is.na(Y)) != 0){ stop("Y contains missings.") }
	if (as.character(class(X)) != "matrix"){ stop("Input (X) is of wrong class.") }
	if (sum(is.na(X)) != 0){ stop("X contains missings.") }
	if (as.character(class(R)) != "matrix"){ stop("Input (R) is of wrong class.") }
	if (sum(is.na(R)) != 0){ stop("R contains missings.") }
	if ( !(as.character(class(maxNoIt)) == "numeric" | as.character(class(maxNoIt)) == "integer") ){ stop("Input (maxNoIt) is of wrong class.") }
	if (as.character(class(minSuccDist)) != "numeric"){ stop("Input (minSuccDist) is of wrong class.") }
	if (as.character(class(corType)) != "character"){ stop("Input (corType) is of wrong class.") }
	if (as.character(class(estType)) != "character"){ stop("Input (estType) is of wrong class.") }
	if (as.character(class(shrinkType)) != "character"){ stop("Input (shrinkType) is of wrong class.") }
	if (as.character(class(hypothesis)) != "character"){ stop("Input (hypothesis) is of wrong class.") }
	if (dim(Y)[1] < 2){ stop("This is a multivariate procedure, provide a multivariate data set.") }
	if (!(corType %in% c("unif", "ar1"))){ stop("corType parameter ill-specified.") }
	if (!(estType %in% c("normal", "robust"))){ stop("estType parameter ill-specified.") }
	if (!(shrinkType %in% c("none", "opt", "full"))){ stop("shrinkType parameter ill-specified.") }
	if (!(hypothesis %in% c("H0", "H1", "H2"))){ stop("hypothesis parameter ill-specified.") }
	if (dim(Y)[1] != dim(X)[1]){ stop("Dimension mismatch between expression and design matrix.") }
	if (dim(X)[2] != dim(R)[2]){ stop("Dimension mismatch between constraint and design matrix.") }
	if (minSuccDist <= 0){ stop("Stopping criterion incorrect.") }
	if (maxNoIt < 1){ stop("Maximum number of iterations smaller than 1.") }

	# transpose Y for compatibility
	Y <- t(Y)
	
	if ( hypothesis=="H0" ){ return(RCMmlH0orH1(Y, X, R, shrinkType=shrinkType, estType=estType, corType=corType, hypothesis=hypothesis)) }
	if ( hypothesis=="H1" ){ return(RCMmlH0orH1(Y, X, R, shrinkType=shrinkType, estType=estType, corType=corType, hypothesis=hypothesis)) }
	if ( hypothesis=="H2" ){ return(RCMmlH2(Y, X, R, maxNoIt, minSuccDist, shrinkType=shrinkType, estType=estType, corType=corType)) }
}
