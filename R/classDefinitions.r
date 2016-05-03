setClass("pathwayFit", representation(Cis="numeric", Trans="matrix", Trans1="matrix", Trans2="matrix", Sigma="numeric", lambda1="matrix", lambdaF="matrix", constr="logical", epsilon="numeric", method="character"))

setClass("cisTest", representation(geneInfo="data.frame", geneId="numeric", comparison="numeric", av.prob1="numeric", av.prob2="numeric", effectSize="numeric", R2="numeric", regId="numeric", beginReg="numeric", endReg="numeric", shrinkage="numeric", p.value="numeric", adjP.value="numeric", analysisType="character", testStatistic="character", nPerm="numeric"))

setClass("entTest", representation(statistic="numeric", p.value="numeric", null.dist="numeric", nperm="numeric", remark="character"))

setClass("miTest", representation(statistic="numeric", p.value="numeric", null.dist="numeric", nperm="numeric", remark="character"))

setClass("rcmFit", representation(betas="numeric", tau2s="numeric", sigma2s="numeric", rho="numeric", av.sigma2s="numeric", shrinkage="numeric", loglik="numeric", corType="character", X="matrix"))

setClass("rcmTest", representation(statistic="numeric", p.value="numeric", betas="numeric", tau2s="numeric", sigma2s="numeric", rho="numeric", av.sigma2s="numeric", shrinkage="numeric", loglik="numeric", nBoot="numeric", corType="character", null.dist="numeric", remark="character"))

setGeneric(name="RCMrandom", def=function(object){standardGeneric("RCMrandom")})

setGeneric(name=".RCMloss", def=function(object, Y){standardGeneric(".RCMloss")})

setGeneric("summary")

setMethod("summary", "rcmTest", function(object){
	###########################################################
	# summarize clm object
	###########################################################
	cat(paste("Coefficients:", sep=""), "\n")
	print(matrix(round(object@betas, digits=3), nrow=1))
	cat(paste("Random effects (as variances):", sep=""), "\n")
	print(matrix(round(object@tau2s, digits=3), nrow=1))
	cat(paste("Variance (average): ", round(object@av.sigma2s, digits=3), sep=""), "\n")
	cat(paste("Correlation (", object@corType, "): ", round(object@rho, digits=3), sep=""), "\n")
	cat(paste("Shrinkage: ", round(object@shrinkage, digits=3), sep=""), "\n")
	cat(paste("Log-likelihood: ", round(object@loglik, digits=3), sep=""), "\n")
	cat("\n")
	cat(paste("Test statistic: ", round(object@statistic, digits=3), ",  p-value: ", round(object@p.value, digits=3), sep=""), "\n")
	cat(paste("Remarks: ", object@remark, sep=""), "\n")
})


setMethod("summary", "rcmFit", function(object){
	###########################################################
	# summarize rcmFit object
	###########################################################
	cat(paste("Coefficients:", sep=""), "\n")
	print(matrix(round(object@betas, digits=3), nrow=1))
	cat(paste("Random effects (as variances):", sep=""), "\n")
	print(matrix(round(object@tau2s, digits=3), nrow=1))
	cat(paste("Variance (average): ", round(object@av.sigma2s, digits=3), sep=""), "\n")
	cat(paste("Correlation (", object@corType, "): ", round(object@rho, digits=3), sep=""), "\n")
	cat(paste("Shrinkage: ", round(object@shrinkage, digits=3), sep=""), "\n")
	cat(paste("Log-likelihood: ", round(object@loglik, digits=3), sep=""), "\n")
})

setMethod("summary", "entTest", function(object){
	###########################################################
	# summarize entTest object
	###########################################################
	cat(paste("Two-sample entropy test:", sep=""), "\n")
	cat(paste("Test statistic: ", round(object@statistic, digits=3), ",  p-value: ", round(object@p.value, digits=3), sep=""), "\n")
	cat(paste("Remarks: ", object@remark, sep=""), "\n")
})

setMethod("summary", "miTest", function(object){
	###########################################################
	# summarize miTest object
	###########################################################
	cat(paste("Mutual information test:", sep=""), "\n")
	cat(paste("Test statistic: ", round(object@statistic, digits=3), ",  p-value: ", round(object@p.value, digits=3), sep=""), "\n")
	cat(paste("Remarks: ", object@remark, sep=""), "\n")
})

setMethod("RCMrandom", "rcmFit", function(object){
	###########################################################
	# generate random data
	###########################################################

	# construct parameters of multivariate normal from which is drawn
	Sigma <- .covMatConstruction(object@sigma2s, object@rho, object@corType)
	In <- matrix(0, dim(object@X)[1], dim(object@X)[1])
	diag(In) <- 1
	Ip <- matrix(0, length(object@sigma2s), length(object@sigma2s))
	diag(Ip) <- 1
	if (length(object@betas) == 1){
		Slh <- (matrix(object@X, ncol=1) %*% matrix(object@tau2s, ncol=1) %*% matrix(object@X, nrow=1))
	} else {
		Tau <- matrix(0, ncol=length(object@tau2s), nrow=length(object@tau2s))
		diag(Tau) <- object@tau2s
		Slh <- object@X %*% Tau %*% t(object@X)
	}
	gge <- matrix(object@X %*% object@betas, ncol=dim(object@X)[1], nrow=length(object@sigma2s), byrow=TRUE)
	em.part1 <- t(rmvnorm(dim(object@X)[1], mean=rep(0, length(object@sigma2s)), sigma=Sigma, method="chol"))
	em.part2 <- t(rmvnorm(length(object@sigma2s), mean=rep(0, dim(object@X)[1]), sigma=Slh, method="svd"))
	em <- em.part1 + t(em.part2)
	return(t(em + gge))
})

setMethod(".RCMloss", "rcmFit", function(object, Y){
	###########################################################
	# function that calculates the log-likelihood of the rcm
	###########################################################
	loss.part2 <- function(ev.cov, ed.Om2, residuals){
		loss <- 0
		for (i in 1:length(ed.Om2$values)){
			loss <- loss - 1/2 * (1 / (ev.cov[1] + ed.Om2$values[i])) * (t(residuals) %*% (matrix(ev.cov[-1], ncol=1) %x% ed.Om2$vectors[,i, drop=FALSE]))^2
		}
		return(as.numeric(loss))
	}

	# first part log-likelihood
	cov.mat <- .covMatConstruction(object@sigma2s, object@rho, object@corType)
	tau.mat <- matrix(0, nrow=length(as.numeric(object@tau2s)), ncol=length(as.numeric(object@tau2s)))
	diag(tau.mat) <- as.numeric(object@tau2s)
	het.mat <- object@X %*% tau.mat %*% t(object@X)
	ev.cov.mat <- eigen(cov.mat, symmetric=TRUE)$values
	ev.het.mat <- eigen(het.mat, symmetric=TRUE)$values
	ev.het.mat[(length(object@tau2s)+1):length(ev.het.mat)] <- 0
	logdet <- 0
	for (i in 1:dim(Y)[2]){
		for (j in 1:dim(Y)[1]){
			logdet <- logdet - 1/2 * log(ev.cov.mat[j] + ev.het.mat[i])
		}
	}	

	# second part log-likelihood
	X.circ <- matrix(rep(1, dim(Y)[1]), ncol=1) %x% object@X
	Y.res.as.vec <- matrix(as.numeric(t(Y)), ncol=1) - X.circ %*% matrix(object@betas, ncol=1)
	ed.cov <- eigen(cov.mat, symmetric=TRUE)
	ed.cov <- rbind(ed.cov$values, ed.cov$vectors)
	ed.Om2 <- eigen(het.mat, symmetric=TRUE)
	logexp <- sum(apply(ed.cov, 2, loss.part2, ed.Om2, Y.res.as.vec))

	ll <- logexp + logdet - prod(dim(Y)) * log(2 * pi) /2

	return(ll)
})

