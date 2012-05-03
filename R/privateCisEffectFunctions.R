	.tuning <- function(data.in, datacgh.in, allest, alphanewmat, powerunbal, null.dists.tune, test.stat="wmw", seqg, nperm, pi0, fdrcut=0.25, nresamp=100, gridnr=30, minim=10, a=a, nosamp=nosamp, shiftsam=shiftsam, verbose){
		################################################################################################
		# function that does the tuning
		################################################################################################
	
		probs <- seq(0,gridnr/(gridnr+1), 1/(gridnr+1))
		probstruncate <- probs[probs<=((100-minim)/100)]
		allest2 <- allest[seqg]
		alphanewmat2<-alphanewmat[seqg,]
		data.in2 <- data.in[seqg,]
		datacgh.in2 <- datacgh.in[seqg,]
		powerunbal2 <- powerunbal[seqg]
		quant <- quantile(powerunbal2, probs=probstruncate)
		whatthrmat <- c()
		for (j in 1:nresamp){
			if ((j %% 50) == 0){ if (verbose){ cat(paste(j, "of", nresamp, "resamples done, and counting...", sep=" "), "\n") } }
			shiftsampled <- sample(shiftsam, nrow(data.in2), replace=TRUE)
			pi0cutoff <- quantile(shiftsampled, pi0)
			shiftsampled <- sapply(shiftsampled, function(x){if (x <= pi0cutoff){ return(0) } else { return(x)}})
			redrawall <- t(sapply(1:nrow(data.in2), .redrawcol, allest=allest2, alphanmat=alphanewmat2, shiftsampled=shiftsampled, data2=data.in2, a=a, nosamp=nosamp))
			newdata <- cbind(datacgh.in2[,-(1:3)], redrawall) 
			if(test.stat=="wmw"){ test.resamp <- apply(newdata, 1, .wmw.test.stats, nosamp, a) } else { test.resamp <- apply(newdata, 1, .wcvm.test.stats, nosamp, a) }
			rawpvals.test.resamp <- .rawps(test.resamp, null.dists.tune, nperm)
			rawp_powerunbal <- cbind(rawpvals.test.resamp, powerunbal2)
			whatthresh <- sapply(quant, .countdiscoveries, rawp_powerunbal=rawp_powerunbal, fdrcutoff=fdrcut)
			whatthrmat <- rbind(whatthrmat, whatthresh)  
		}
		whatthrmean <- apply(whatthrmat, 2, mean)
		# print(whatthrmean)
		bestthr <- which.max(whatthrmean)    
		return(quant[bestthr])
	}

	.uni.an.wcvm <- function(data.both, nosamp, a, nperm, low.ci.thres=0.10, verbose=TRUE){
		################################################################################################
		# function that performs a univariate analysis with the weighted cvm-statistic,
		# with the efficient p-value calculation.
		################################################################################################

		# Calculate the observed weighted cvm test statistics
		wcvm.obs <- apply(data.both, 1, .wcvm.test.stats, nosamp, a)

		# Actual analysis
		data.perm <- data.both
		# data.perm <- rbind(data.perm, data.perm)
		total.genes.on.chr <- dim(data.both)[1]
		remainders <- c(1:dim(data.perm)[1])
		steps <- sort(unique(c(0,25,50,100,150,200,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,seq(from=10000,to=50000,by=1000),nperm)))
		steps <- steps[steps <= nperm]
		wcvm.mat <- c()
		for(j in 1:(length(steps)-1)){
			for(i in 1:(steps[j+1]-steps[j])){
				if (verbose){  if (((steps[j]+i) %% 100) == 0){ cat(paste(steps[j]+i," of ", steps[length(steps)]," permutations done, and counting...", sep=""), "\n") } }
				x <- sample(1:nosamp,nosamp) + a*nosamp
				data.ran <- cbind(data.perm[, c(1:(a*nosamp))], data.perm[, x])
				wcvm.ran <- apply(data.ran, 1, .wcvm.test.stats, nosamp, a)
				wcvm.mat <- cbind(wcvm.mat, wcvm.ran)
			}

			# compute pval bound and delete row from data set and permutation set when 0.001 < lower bound
			perm.and.obs <- cbind(wcvm.mat, wcvm.obs)
			pvals <- apply(perm.and.obs, 1, .countth)/steps[j+1]
			pbound <- sapply(pvals, .pvalbound, steps[j+1])
			index <- cbind(1:length(pbound), pbound)
			pind <- index[pbound < low.ci.thres, 1]

			if (length(pind) > 2){ remainders <- remainders[pind] } else { pind <- c(1:length(remainders)) }

			wcvm.mat <- wcvm.mat[pind, ]
			wcvm.obs <- wcvm.obs[pind]
		
			if (verbose){   cat(paste(steps[j]+i, "of", steps[length(steps)], " permutations done, and", length(remainders), "of", total.genes.on.chr, "genes remaining...", sep=" "), "\n")   }
			data.perm <- data.both[remainders, , drop=FALSE]
		}

		# Generate list with raw p-values for all genes 
		raw.pvals <- cbind(c(1:dim(data.both)[1]), rep(1,dim(data.both)[1]))
		raw.pvals[remainders,2] <- pvals[pind]

		# Calculate BH-adjusted p-values
		adjpvals <- cbind(raw.pvals, p.adjust(raw.pvals[,2], "BH"))
		adjpvals[,2:3] <- round(adjpvals[,2:3], digits=4)

		# Some editing for presentation purposes
		colnames(adjpvals) <- c("clone.id", "raw.p", "adj.p")

		return(adjpvals)
	}

	.uni.an.prob <- function(data.both, nosamp, a, nperm, low.ci.thres=0.10, verbose=TRUE){
		################################################################################################
		# function that performs a univariate analysis with the normalized weighted wm-statistic,
		# with the efficient p-value calculation.
		################################################################################################

		# Calculate the observed weighted cvm test statistics
		prob.obs <- apply(data.both, 1, .prob.test.stats, nosamp, a)

		# Actual analysis
		data.perm <- data.both
		total.genes.on.chr <- dim(data.both)[1]
		remainders <- c(1:dim(data.perm)[1])
		steps <- sort(unique(c(0,25,50,100,150,200,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,seq(from=10000,to=50000,by=1000),nperm)))
		steps <- steps[steps <= nperm]
		prob.mat <- c()
		for(j in 1:(length(steps)-1)){
			for(i in 1:(steps[j+1]-steps[j])){
				if (verbose){  if (((steps[j]+i) %% 100) == 0){ cat(paste(steps[j]+i," of ", steps[length(steps)]," permutations done, and counting...", sep=""), "\n") }  }
				x <- sample(1:nosamp,nosamp) + a*nosamp
				data.ran <- cbind(data.perm[,c(1:(a*nosamp))],data.perm[,x])
				prob.ran <- apply(data.ran, 1, .prob.test.stats, nosamp, a)
				prob.mat <- cbind(prob.mat, prob.ran)
			}

			# compute pval bound and delete row from data set and permutation set when 0.001 < lower bound
			perm.and.obs <- cbind(prob.mat, prob.obs)
			pvals <- apply(perm.and.obs, 1, .countth)/steps[j+1]
			pbound <- sapply(pvals, .pvalbound, steps[j+1])
			index <- cbind(1:length(pbound), pbound)
			pind <- index[pbound < low.ci.thres,1]

			if (length(pind) > 2){ remainders <- remainders[pind] } else { pind <- c(1:length(remainders)) }

			prob.mat <- prob.mat[pind,]
			prob.obs <- prob.obs[pind]

			# cat(paste(length(remainders), "of", total.genes.on.chr, "genes remaining...", sep=" "), "\n")
			if (verbose){   cat(paste(steps[j]+i, "of", steps[length(steps)], " permutations done, and", length(remainders), "of", total.genes.on.chr, "genes remaining...", sep=" "), "\n")  }
			data.perm <- data.both[remainders, , drop=FALSE]
	    	}

		# Generate list with raw p-values for all genes    
		raw.pvals <- cbind(c(1:dim(data.both)[1]), rep(1,dim(data.both)[1]))
		raw.pvals[remainders, 2] <- pvals[pind]

		# Calculate BH-adjusted p-values
		adjpvals <- cbind(raw.pvals, p.adjust(raw.pvals[,2], "BH"))
		adjpvals[,2:3] <- round(adjpvals[,2:3], digits=4)

		# Some editing for presentation purposes
		colnames(adjpvals) <- c("clone.id", "raw.p", "adj.p")
	
		return(adjpvals)
	}


	.R2.stat <- function(data.both, a, nosamp){
		################################################################################################
		# function that calculates the R^2 value.
		################################################################################################

		# Define function for calculated of numerator of R^2 statistic
		r2.numerator <- function(data.both, a, nosamp){
			alpha.ind <- matrix(data.both[c(1:(a*nosamp))], ncol=a, byrow=TRUE)
			alpha.1.mom <- apply(alpha.ind, 2, mean)
			alpha.2.mom <- t(alpha.ind) %*% alpha.ind / nosamp

			cgh.em <- cbind(matrix(data.both[c(1:(a*nosamp))], ncol=2, byrow=TRUE), data.both[c((a*nosamp+1):((a+1)*nosamp))])

			mu1 <- 1 / nosamp*sum(cgh.em[,3]*(cgh.em[,1]*alpha.1.mom[2] - alpha.2.mom[2,1]) / (alpha.2.mom[1,1]*alpha.1.mom[2]-alpha.1.mom[1]*alpha.2.mom[1,2]))
			mu2 <- 1 / nosamp*sum(cgh.em[,3]*(cgh.em[,2]*alpha.1.mom[1] - alpha.2.mom[2,1]) / (alpha.2.mom[2,2]*alpha.1.mom[1]-alpha.1.mom[2]*alpha.2.mom[1,2]))
			call.means <- matrix(c(mu1, mu2), nrow=1)

        		return(sum((data.both[a*nosamp+c(1:nosamp)] - alpha.ind %*% t(call.means))^2)/(nosamp-1))
		}

		# Calculate R^2 statistic
		r2.denominator <- apply(data.both[, c((a*nosamp+1):((a+1)*nosamp))], 1, var)
		numerator <- apply(data.both, 1, r2.numerator, a=a, nosamp=nosamp)
		R2 <- 1- numerator / r2.denominator
		for (i in 1:length(R2)){ R2[i] <- min(1, max(0, R2[i])) }
		return(R2)
	}


	.shrin.an.prob <- function(data.both, nosamp, a, nperm, low.ci.thres=0.1, datacgh.org, verbose=TRUE){
		################################################################################################
		# function that performs a regional analysis with the normalized weighted mw-statistic, also
		# called the probability in this context, with shrinkage per region and efficient p-value calculation.
		################################################################################################

		shrunken.prob.test.stats <- function(reg.bounds.lambda, data.both, nosamp, a){
			################################################################################################
			# Function that reshuffles the shrunken test statistics to get them in the right format.
			################################################################################################

			shrunken.test.stats.wrong.format <- apply(reg.bounds.lambda, 1, "shrunken.prob.test.stats.per.reg", cgh.em=data.both, nosamp=nosamp, a=a)
			if (is.numeric(shrunken.test.stats.wrong.format)){ 
				shrunken.test.stats <- shrunken.test.stats.wrong.format
			}
			if (is.matrix(shrunken.test.stats.wrong.format)){
				shrunken.test.stats <- as.numeric(shrunken.test.stats.wrong.format)
			}
			if (is.list(shrunken.test.stats.wrong.format)){
				shrunken.test.stats <- NULL
				for (i in 1:dim(reg.bounds.lambda)[1]){
					shrunken.test.stats <- c(shrunken.test.stats, shrunken.test.stats.wrong.format[[i]])
				}
			}
			return(shrunken.test.stats)
		}

		shrunken.prob.test.stats.per.reg <- function(bounds.lambda, cgh.em, nosamp, a){
			################################################################################################
			# Function that calculates the normalized weighted MW test statistics for all clones in a region.
			################################################################################################
	
			lambda <- bounds.lambda[3]
			bounds <- bounds.lambda[c(1:2)]
			if (bounds[1] != bounds[2]){
				cgh.em <- cgh.em[c(bounds[1]:bounds[2]),]
				marg.test.stats <- apply(cgh.em,1, .prob.test.stats, nosamp,a)
				shrunken.test.stats <- lambda*marg.test.stats + (1-lambda)*mean(marg.test.stats)
			} else {
				shrunken.test.stats <- .wcvm.test.stats(cgh.em[bounds[1],],nosamp,a)
			}
			return(shrunken.test.stats)
		}

		if (verbose){   cat("construct regions...", "\n")  }
		reg.bounds <- .find.cgh.reg.data(data.both, a, nosamp, datacgh.org)

		if (verbose){   cat("calculate shrinkage parameters...", "\n")   }
		lambda.of.reg <- apply(reg.bounds, 1, .lambda.per.reg, data.both=data.both, a=a, nosamp=nosamp)
		reg.bounds.lambda <- cbind(reg.bounds, lambda.of.reg)
		colnames(reg.bounds.lambda) <- NULL

		if (verbose){   cat("calculate observed test statistics...", "\n")  }
		shrunken.prob.obs <- shrunken.prob.test.stats(reg.bounds.lambda, data.both, nosamp, a)

		if (verbose){   cat("calculate null distribution...", "\n")    }
		prob.obs <- shrunken.prob.obs
		data.perm <- data.both
		total.genes.on.chr <- dim(data.both)[1]
		remainders <- c(1:dim(data.perm)[1])
		steps <- sort(unique(c(0,25,50,100,150,200,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,seq(from=10000,to=50000,by=1000),nperm)))
		steps <- steps[steps <= nperm]
		prob.mat <- c()
		for(j in 1:(length(steps)-1)){
			for(i in 1:(steps[j+1]-steps[j])){
				if (verbose){   if (((steps[j]+i) %% 100) == 0){ cat(paste(steps[j]+i," of ", steps[length(steps)]," permutations done, and counting...", sep=""), "\n") }  }
				x <- sample(1:nosamp,nosamp) + a*nosamp
				data.ran <- cbind(data.perm[, c(1:(a*nosamp))], data.perm[,x])
				prob.ran <- shrunken.prob.test.stats(reg.bounds.lambda, data.ran, nosamp, a)
				prob.mat <- cbind(prob.mat, prob.ran[remainders])
			}
			# compute pval bound and delete row from data set and permutation set when 0.001 < lower bound
			perm.and.obs <- cbind(prob.mat, prob.obs)
			pvals <- apply(perm.and.obs, 1, .countth)/steps[j+1]
			pbound <- sapply(pvals, .pvalbound, steps[j+1])
			index <- cbind(1:length(pbound), pbound)
			pind <- index[pbound < low.ci.thres, 1]

			if (length(pind) > 2){ remainders <- remainders[pind] } else { pind <- c(1:length(remainders)) }

			prob.mat <- prob.mat[pind, , drop=FALSE]
			prob.obs <- prob.obs[pind]

			if (verbose){   cat(paste(steps[j]+i, "of", steps[length(steps)], " permutations done, and", length(remainders), "of", total.genes.on.chr, "genes remaining...", sep=" "), "\n")  }
		}

		# Generate list with raw p-values for all genes    
		raw.pvals <- cbind(c(1:dim(data.both)[1]), rep(1,dim(data.both)[1]))
		raw.pvals[remainders,2] <- pvals[pind]

		# Calculate BH-adjusted p-values
		adjpvals <- cbind(raw.pvals, p.adjust(raw.pvals[,2], "BH"))
		adjpvals[,2:3] <- round(adjpvals[,2:3], digits=4)

		reg.details <- NULL
		for (i in 1:dim(reg.bounds.lambda)[1]){
			reg.length <- (reg.bounds.lambda[i,2] - reg.bounds.lambda[i,1] +1)
			reg.details <- rbind(reg.details, cbind(rep(i,reg.length), matrix(rep(reg.bounds.lambda[i,], reg.length), nrow=reg.length, byrow=TRUE)))
		}
		adjpvals <- cbind(adjpvals[,1], reg.details, adjpvals[,2:3])
		colnames(adjpvals) <- c("clone.id", "reg.id", "begin.reg", "end.reg", "shrinkage", "raw.p", "adj.p")

		return(adjpvals)
	}

	.shrin.an.wcvm <- function(data.both, nosamp, a, nperm, low.ci.thres=0.1, datacgh.org, verbose=TRUE){
		################################################################################################
		# function that performs a regional analysis with the weighted cvm-statistic,
		# with shrinkage per region and efficient p-value calculation.
		################################################################################################
	
		shrunken.wcvm.test.stats <- function(reg.bounds.lambda, data.both, nosamp, a){
			################################################################################################
			# Function that reshuffles the shrunken test statistics to get them in the right format.
			################################################################################################
	
			shrunken.test.stats.wrong.format <- apply(reg.bounds.lambda, 1, "shrunken.wcvm.test.stats.per.reg", cgh.em=data.both, nosamp=nosamp, a=a)
			if (is.numeric(shrunken.test.stats.wrong.format)){ 
				shrunken.test.stats <- shrunken.test.stats.wrong.format
			}
			if (is.matrix(shrunken.test.stats.wrong.format)){
				shrunken.test.stats <- as.numeric(shrunken.test.stats.wrong.format)
			}
			if (is.list(shrunken.test.stats.wrong.format)){
				shrunken.test.stats <- NULL
				for (i in 1:dim(reg.bounds.lambda)[1]){
					shrunken.test.stats <- c(shrunken.test.stats, shrunken.test.stats.wrong.format[[i]])
				}
			}
			return(shrunken.test.stats)
		}
	
		shrunken.wcvm.test.stats.per.reg <- function(bounds.lambda, cgh.em, nosamp, a){
			################################################################################################
			# Function that calculates the weighted CvMe test statistics for all clones in a region.
			################################################################################################
	
			lambda <- bounds.lambda[3]
			bounds <- bounds.lambda[c(1:2)]
			if (bounds[1] != bounds[2]){
 				cgh.em <- cgh.em[c(bounds[1]:bounds[2]),]
				marg.test.stats <- apply(cgh.em, 1, .wcvm.test.stats, nosamp, a)
				shrunken.test.stats <- lambda * marg.test.stats + (1-lambda) * mean(marg.test.stats)
			} else {
				shrunken.test.stats <- .wcvm.test.stats(cgh.em[bounds[1],], nosamp, a)
			}
			return(shrunken.test.stats)
		}
    
		if (verbose){   cat("construct regions...", "\n")  }
		reg.bounds <- .find.cgh.reg.data(data.both, a, nosamp, datacgh.org)

		if (verbose){   cat("calculate shrinkage parameters...", "\n")   }
		lambda.of.reg <- apply(reg.bounds, 1, .lambda.per.reg, data.both=data.both, a=a, nosamp=nosamp)
		reg.bounds.lambda <- cbind(reg.bounds, lambda.of.reg)
		colnames(reg.bounds.lambda) <- NULL

		if (verbose){   cat("calculate observed test statistics...", "\n")   }
		shrunken.wcvm.obs <- shrunken.wcvm.test.stats(reg.bounds.lambda, data.both, nosamp, a)

		if (verbose){   cat("calculate null distribution...", "\n")  }
		wcvm.obs <- shrunken.wcvm.obs
		data.perm <- data.both
		total.genes.on.chr <- dim(data.both)[1]
		remainders <- c(1:dim(data.perm)[1])
		steps <- sort(unique(c(0,25,50,100,150,200,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,seq(from=10000,to=50000,by=1000),nperm)))
		steps <- steps[steps <= nperm]
		wcvm.mat <- c()
		for(j in 1:(length(steps)-1)){
			for(i in 1:(steps[j+1]-steps[j])){
				if (verbose){   if (((steps[j]+i) %% 100) == 0){ cat(paste(steps[j]+i," of ", steps[length(steps)], " permutations done, and counting...", sep=""), "\n") }   }
				x <- sample(1:nosamp,nosamp) + a * nosamp
				data.ran <- cbind(data.perm[, c(1:(a*nosamp))],data.perm[,x])
				wcvm.ran <- shrunken.wcvm.test.stats(reg.bounds.lambda, data.ran, nosamp, a)
				wcvm.mat <- cbind(wcvm.mat, wcvm.ran[remainders])
			}
			# compute pval bound and delete row from data set and permutation set when 0.001 < lower bound
			perm.and.obs <- cbind(wcvm.mat, wcvm.obs)
			pvals <- apply(perm.and.obs, 1, .countth)/steps[j+1]
			pbound <- sapply(pvals, .pvalbound, steps[j+1])
			index <- cbind(1:length(pbound), pbound)
			pind <- index[pbound < low.ci.thres,1]
	
			if (length(pind) > 2){ remainders <- remainders[pind] } else { pind <- c(1:length(remainders)) }

			wcvm.mat <- wcvm.mat[pind, , drop=FALSE]
			wcvm.obs <- wcvm.obs[pind]

			if (verbose){   cat(paste(steps[j]+i, "of", steps[length(steps)], " permutations done, and", length(remainders), "of", total.genes.on.chr, "genes remaining...", sep=" "), "\n")   }
		}

		# Generate list with raw p-values for all genes    
		raw.pvals <- cbind(c(1:dim(data.both)[1]), rep(1,dim(data.both)[1]))
		raw.pvals[remainders,2] <- pvals[pind]

		# Calculate BH-adjusted p-values
		adjpvals <- cbind(raw.pvals, p.adjust(raw.pvals[,2], "BH"))
		adjpvals[,2:3] <- round(adjpvals[,2:3], digits=4)

		reg.details <- NULL
		for (i in 1:dim(reg.bounds.lambda)[1]){
			reg.length <- (reg.bounds.lambda[i,2] - reg.bounds.lambda[i,1] +1)
			reg.details <- rbind(reg.details, cbind(rep(i, reg.length), matrix(rep(reg.bounds.lambda[i,], reg.length), nrow=reg.length, byrow=TRUE)))
		}
		adjpvals <- cbind(adjpvals[,1], reg.details, adjpvals[,2:3])
		colnames(adjpvals) <- c("clone.id", "reg.id", "begin.reg", "end.reg", "shrinkage", "raw.p", "adj.p")

		return(adjpvals)
	}

	.lambda.per.reg <- function(bounds, data.both, a, nosamp){
		################################################################################################
		# function that calculates the shrinkage factor for a region.
		################################################################################################

		# if the region consists of more than one clone, 
		if (bounds[1] != bounds[2]){
			# select data from region
			data.both.org <- matrix(data.both[c(bounds[1]:bounds[2]),], nrow=(bounds[2]-bounds[1]+1))
			data.both.org <- data.both.org[, c((dim(data.both)[2] - dim(data.both)[2] / 3 + 1):dim(data.both)[2])]
	        
			# calculate mean correlation, and shrinkage parameter
			slh <- cor(t(data.both.org), method="spearman") 
			lambda.final <- 1 - max(0, mean(slh[upper.tri(slh)]))
		}
	
		# if region consists of one clone, set shrinkage parameter equal to 1.
		if (bounds[1] == bounds[2]){ lambda.final <- 1 }
		return(lambda.final)
	}

	.find.cgh.reg.data <- function(data.both, a, nosamp, datacgh.org){
		################################################################################################
		# Compress the data to regions. A new regions starts when the CGH probabilities change.
		################################################################################################

		# every change in the call probabilities implies a new region
		cgh.data <- data.both[,c(1:(a*nosamp))]
    
		splitter <- list()
		splitter[[1]] <- c(1)
		index.temp <- 1
		j <- 1
		for (i in 1:(dim(cgh.data)[1]-1)){
			if (all(cgh.data[i,] == cgh.data[i+1,])){
				index.temp <- c(index.temp,i+1)
				splitter[[j]] <- index.temp
			
			} else {
				index.temp <- i+1
				j <- j + 1
				splitter[[j]] <- index.temp
			}
		}

		region.details <- NULL
		for (i in 1:length(splitter)){
			region.details <- rbind(region.details, c(min(splitter[[i]]), max(splitter[[i]])))
		}
		return(region.details)
	}

	.pvalbound <- function(pval, np){
		################################################################################################
		# p-value lower conf. bound. 3.09 = Z_{0.001}
		################################################################################################
	
		return(pval-sqrt(pval * (1-pval) / np) * 3.09)
	}

	.countth.eff <- function(statlist, threshold){
		################################################################################################
		# counts the number of values in statlist exceed threshold.
		################################################################################################
	
		return(length(statlist[statlist >= threshold]))  
	}

	.countth <- function(statlist){
		################################################################################################
		# Counts the number of values in statlist exceed threshold.
		################################################################################################

		threshold <- as.numeric(statlist[length(statlist)])
		statlist <- as.numeric(statlist[c(1:(length(statlist)-1))])
		return(length(statlist[statlist >= threshold]))
	}

	.pval.perm.marg <- function(observed, permuted, nperm){
		################################################################################################
		# Calculates the marginal p-values.
		# The p-value of gene is calculated using the null-ditribution 
		# resulting from the permutations of that gene.
		################################################################################################

		perm.and.obs <- cbind(permuted, observed)
		return(apply(perm.and.obs, 1, .countth) / nperm)
	}

	.rawps <- function(stats.obs, nulldists, nperm){
		################################################################################################
		# calculate the raw p-values using the observed test statistics
		# and their null distribution.
		################################################################################################

		pval.ln <- .pval.perm.marg(stats.obs, nulldists, nperm)
		return(pval.ln)
	}

	.wcvm.test.stats <- function(cgh.em, nosamp, a){
		################################################################################################
		# function that calculates the wcvm test statistics for one clone.
		################################################################################################
	
		# makes a matrix of call probs
		cgh.2cat <- matrix(cgh.em[c(1:(a*nosamp))], ncol=a, byrow=TRUE)
	
		# calculate contrast coefficients
		alphaas <- t(cgh.2cat) %*% cgh.2cat
		cs <- as.numeric(solve(alphaas) %*% matrix(c(-1,1), ncol=1))

		cgh.em <- cbind(cgh.2cat, cgh.em[c((a*nosamp+1):((a+1)*nosamp))])
		cgh.em <- cbind(cgh.em[order(cgh.em[,3]),], rep(1/dim(cgh.em)[1], dim(cgh.em)[1]))
		cgh.em <- cbind(cgh.em, cumsum(cgh.em[, 4]), cumsum(cgh.em[, 1] * cgh.em[, dim(cgh.em)[2]]) * cs[1], cumsum(cgh.em[, 2] * cgh.em[, dim(cgh.em)[2]]) * cs[2])

		test.stat <- -sum(cgh.em[,7] + cgh.em[,6]) / nosamp
		return(test.stat)
	}

	.wmw.test.stats <- function(cgh.em, nosamp, a){
		################################################################################################
		# function that calculates the wMW-like test statistics for one clone.
		################################################################################################

		# makes a matrix; first 2 call probs, last is expression
		cgh.em <- cbind(matrix(cgh.em[c(1:(a*nosamp))], ncol=a,byrow=TRUE), cgh.em[c((a*nosamp+1):((a+1)*nosamp))]) 

		# sort by expression
		cgh.em <- cgh.em[order(cgh.em[,(a+1)]),]  

		# shift probs 1 position, because I_{Xi < Xj}
		cgh.em <- cbind(cgh.em, rbind(rep(0,a), apply(cgh.em[,1:a], 2, cumsum)[-nosamp,])) 

		# '2' = prob. class 1, a+2 = cum. prob class 2
		test.stat <- sum(cgh.em[,a+2]*cgh.em[,2]) 
		return(test.stat)
	}

	.prob.test.stats <- function(cgh.em, nosamp, a){
		################################################################################################
		# function that calculates the wMW-like test statistics for one clone.
		################################################################################################

		# makes a matrix of call probs
		cgh.2cat <- matrix(cgh.em[c(1:(a*nosamp))], ncol=a, byrow=TRUE)
	
		# calculate contrast coefficients
		alphaas <- t(cgh.2cat) %*% cgh.2cat
		cs <- c(det(alphaas), alphaas[1,2]*sum(alphaas)/2)

		# makes a matrix; first 2 call probs, last is expression
		cgh.em <- cbind(matrix(cgh.em[c(1:(a*nosamp))], ncol=a,byrow=TRUE), cgh.em[c((a*nosamp+1):((a+1)*nosamp))]) 

		# sort by expression
		cgh.em <- cgh.em[order(cgh.em[,(a+1)]),]  

		# shift probs 1 position, because I_{Xi < Xj}
		cgh.em <- cbind(cgh.em,rbind(rep(0,a), apply(cgh.em[,1:a], 2, cumsum)[-nosamp,])) 
	
		# '2' = prob. class 1, a+2 = cum. prob class 2
		test.stat <- (sum(cgh.em[,a+2] * cgh.em[,2]) - cs[2]) / cs[1]
		return(test.stat)
	}

	.pretest <- function(alphascgh){
		################################################################################################
		# function that determines which hypothesis is tested (loss vs. no-loss or no-gain vs. gain).
		# also merges the call probabilities of the aberrated class that is not-dominant.
		################################################################################################

		alphas <- alphascgh[1:3]
		probs <- matrix(alphascgh[-(1:3)], ncol=3, byrow=TRUE)
		if (alphas[1] >= alphas[3]){
			probs2 <- cbind(probs[,1], probs[,2]+probs[,3])
			alphas2 <- c(alphas[1], alphas[2]+alphas[3])
			return(c(1, alphas2, as.vector(t(probs2))))
		} else {
			probs2 <- cbind(probs[,1]+probs[,2],probs[,3])
			alphas2 <- c(alphas[1]+alphas[2], alphas[3])
			return(c(2, alphas2, as.vector(t(probs2))))
		} 
	}

	.alphaest <- function(cgh.em, nosamp, a){
		################################################################################################
		# function that calculated the first order moments of the call probabilities.
		################################################################################################
	
		cgh.em <- matrix(cgh.em[c(1:(a*nosamp))], ncol=a,byrow=TRUE) 
		return(apply(cgh.em,2,mean))
	}

	.alphabivariate <- function(cgh.em, nosamp, a){
		################################################################################################
		# function that calculated the second order moments of the call probabilities.
		################################################################################################
	    
		cgh.em <- matrix(cgh.em[c(1:(a*nosamp))], ncol=a, byrow=TRUE) 
		cgh.em1 <- cgh.em[,1]
		cgh.em2 <- cgh.em[,2]
		return(c(1 / nosamp * (cgh.em1 %*% cgh.em1), 1 / nosamp * (cgh.em2 %*% cgh.em2), 1 / nosamp * (cgh.em1 %*% cgh.em2)))
	}

	.probs2calls <- function(problist){
		################################################################################################
		# function that converts sof calls to hard calls
		################################################################################################

		maxprobposition <- which.max(problist)
		call.list <- rep(0, length(problist))
		call.list[maxprobposition] <- 1
		return(call.list)
	}

	.shift.est <- function(row, nosamp, a, data2, alphabmat, alphanmat, minalphathr){
		################################################################################################
		# function that estimates the effect size.
		# if a gene exceeds the power unbalance threshold NA is returned.
		################################################################################################

		cgh.em <- cbind(matrix(data2[row,c(1:(a*nosamp))], ncol=a, byrow=TRUE), data2[row,c((a*nosamp+1):((a+1)*nosamp))])
		alphasbiv <- alphabmat[row,]
		alphas <- alphanmat[row,]
		c1 <- (alphasbiv[1]/alphasbiv[3] - alphas[1]/alphas[2])^(-1)
		c2 <- (alphasbiv[2]/alphasbiv[3] - alphas[2]/alphas[1])^(-1)
		if (is.na(c1) | is.na(c2)) {
			return(NA)
		} else {
			if (max(c1,c2) >= minalphathr) {
				return(NA)
			} else {
				mu1 <- 1 / nosamp * sum(cgh.em[,3] * (cgh.em[,1] * alphas[2] - alphasbiv[3]) / (alphasbiv[1] * alphas[2] - alphas[1] * alphasbiv[3]))
				mu2 <- 1 / nosamp * sum(cgh.em[,3] * (cgh.em[,2] * alphas[1] - alphasbiv[3]) / (alphasbiv[2] * alphas[1] - alphas[2] * alphasbiv[3]))
				shiftest <- mu2 - mu1
				return(shiftest)
			}
		}
	}
   
	.powerunbalance <- function(row, alphabmat){
		################################################################################################
		# function that calculates the power unbalance
		################################################################################################

		alphasbiv <- alphabmat[row,]
		unbalance <- alphasbiv[1] * alphasbiv[2] - alphasbiv[3]^2
		return(unbalance)
	}	

	.countdiscoveries <- function(powerquant, rawp_powerunbal, fdrcutoff){
		################################################################################################
		# function that converts sof calls to hard calls
		################################################################################################

		rawp_filt <- rawp_powerunbal[rawp_powerunbal[,2]>=powerquant,][,1]
		m <- length(rawp_filt)
		adpjrawp_filt <- cbind(rawp_filt, p.adjust(rawp_filt, "BH"))
		selected <- matrix(adpjrawp_filt[adpjrawp_filt[,2] <= fdrcutoff,],ncol=2)
		count <- nrow(selected) 
		false <- ifelse(count > 0, max(selected[,1])*m, 0)
		countS <- count - false
		return(countS)
	}

	.redraw <- function(col, shiftest, shiftsam, alpharow, dataexp, datacgh){
		################################################################################################
		# function that .....
		################################################################################################
	
		if(is.na(shiftest)){
			return(sample(dataexp,1)) 
		} else {  
			# this returns a re-sample from the data when there is basically only one group
			alpha <- alpharow[1]
			Ig <- sample(c(1,2), 1, prob=c(alpha,1-alpha))
			if (Ig==1){
				xib <- sample(dataexp, 1, prob=datacgh[,1]) + shiftest 
			} else {
				xib <- sample(dataexp, 1, prob=datacgh[,2])
			}
			Ig2 <- sample(c(1,2), 1, prob = c( max(0, min(1,datacgh[col,1])), max(0, min(1, 1-datacgh[col,1]))))
			if (Ig2==1){
				xib2 <- xib
			} else {
				xib2 <- xib + shiftsam
			}
			return(xib2)
		}
	}

	.redrawcol <- function(row, allest, alphanmat, shiftsampled, data2, a, nosamp){
		################################################################################################
		# function that .......
		################################################################################################
	
		shiftest <- allest[row]
		datacgh <- matrix(data2[row,c(1:(a*nosamp))],ncol=a,byrow=TRUE)
		dataexp <- data2[row,c((a*nosamp+1):((a+1)*nosamp))]
		shiftsam <- shiftsampled[row]
		alpharow <- alphanmat[row,]
		result <- sapply(1:nosamp, .redraw, shiftest=shiftest, alpharow=alpharow, shiftsam=shiftsam, dataexp=dataexp, datacgh=datacgh)
		return(result)
	}

	.pi0est <- function(data.in, null.dists.tune, test.stat="wmw", seqg, nperm, a, nosamp){
		################################################################################################
		# function that estimates the proportion of rejected null-hypothesis, i.e., the number
		# of gene whose expression is afffected by copy number changes.
		################################################################################################

		data.in2 <- data.in[seqg,]
		if(test.stat=="wmw"){ 
			test.stats <- apply(data.in2, 1, .wmw.test.stats, nosamp, a)  
		} else { 
			test.stats <- apply(data.in2, 1, .wcvm.test.stats, nosamp, a) 
		}
		rawpvals.test <- .rawps(test.stats, null.dists.tune, nperm)
		pi0 <- convest(rawpvals.test)
		return(pi0)
	}

	.datareduce <- function(powerunbal, unbalthr){
		################################################################################################
		# function that selects rows that pass the power unbalance criterion
		################################################################################################
	
		datarows <- which(powerunbal >= unbalthr)
		return(datarows)
	}
	
	.nulldist.all.wcvm <- function(data.both, nosamp, a, nperm, verbose){
		################################################################################################
		# function that calculates the null distribution of the weighted CvM test statistics for tuning.
		################################################################################################
	    
		# Permutes data and calculates test statistic on permuted data
		cvm.like.mat <- NULL
		for(i in 1:nperm){
			if ((i %% 50) == 0){ if (verbose){ cat(i," of ", nperm, " permutations done, and counting...", "\n") } }
			x <- sample(1:nosamp,nosamp) + a*nosamp
			data.ran <- cbind(data.both[,c(1:(a*nosamp))], data.both[,x])
			cvm.like.ran <- apply(data.ran, 1, .wcvm.test.stats, nosamp, a)
			cvm.like.mat <- cbind(cvm.like.mat,cvm.like.ran)
		}
		return(cvm.like.mat)
	}
	
	.nulldist.all.wmw <- function(data.both, nosamp, a, nperm, verbose){
		################################################################################################
		# function that calculates the null distribution of the weighted CvM test statistics for tuning.
		################################################################################################
	    
		# Permutes data and calculates test statistic on permuted data
		wmw.ln.mat <- NULL
		cghdata <- data.both[,c(1:(a*nosamp))]
		for(i in 1:nperm){
			if ((i %% 50) == 0){ if (verbose){ cat(i," of ", nperm, " permutations done, and counting...", "\n") } }
			x <- sample(1:nosamp, nosamp) + a*nosamp
			data.ran <- cbind(cghdata, data.both[,x])
			wmw.ran <- apply(data.ran, 1, .wmw.test.stats, nosamp, a)
			wmw.ln.mat <- cbind(wmw.ln.mat, wmw.ran)
		}
		return(wmw.ln.mat)
	}

