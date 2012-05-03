# source("")
# cisEffectTuning(CNdata, GEdata, "wcvm", nGenes=150, nPerm=250, minCallProbMass=0.10, verbose=TRUE)
#  testRes <- cisEffectTesting(CNdata, GEdata, genes2test, 1, "univariate", "wcvm", nPerm=100)


cisEffectPlot <- function (geneId, CNdata, GEdata, verbose=FALSE){
	################################################################################################
	# plot expression vs. copy number
	################################################################################################

	# input checks
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	if ( !(as.character(class(GEdata)) == "ExpressionSet" | as.character(class(GEdata)) == "eSet") ){ stop("GEdata not of class ExpressionSet.") }

	# actual start testing	
	if (verbose){ cat("formatting data for testing...", "\n") }
	CNdata <- cghCall2subset(CNdata, rep(geneId, 2), verbose)
	GEdata <- ExpressionSet2subset(GEdata, rep(geneId, 2), verbose)

	# extract data from objects
	nosamp <- ncol(exprs(GEdata))
	cghdata.probs <- numeric()
	for (i in 1:dim(calls(CNdata))[2]){
		if (is.null(probamp(CNdata)[,i])){
			cghdata.probs <- cbind(cghdata.probs, cbind(probloss(CNdata)[,i], probnorm(CNdata)[,i], probgain(CNdata)[,i]))
		} else {
			cghdata.probs <- cbind(cghdata.probs, cbind(probloss(CNdata)[,i], probnorm(CNdata)[,i], probgain(CNdata)[,i] + probamp(CNdata)[,i]))
		}
	}
	nclass <- dim(cghdata.probs)[2] / dim(calls(CNdata))[2]
	annInfo <- fData(GEdata)
	# chrInfo <- annInfo[, GEchr]

	# Estimate alpha parameters per clone, a0 classes
	if (verbose){ cat("pre-test...", "\n") }
	alphamat <- t(apply(cghdata.probs, 1, .alphaest, nosamp=nosamp, a=nclass))

	# Perform pre-test and merge columns
	datacgh2 <- as.matrix(t(apply(cbind(alphamat, cghdata.probs), 1, .pretest)))
	lossorgain <- datacgh2[,1]
	alphas2 <- datacgh2[,2:3]
	datafortest <- as.matrix(cbind(datacgh2, exprs(GEdata))[,-(1:3)])

	# Estimate new alpha and bivariate alpha parameters per clone
	alphabivmat <- t(apply(datafortest, 1, .alphabivariate, nosamp=nosamp, a=2))
	colnames(alphabivmat) <- c("a11", "a22", "a12")

	# gene.id <- which(data.tuned$genestotest==geneId)
	cgh.em <- datafortest[1, ]
	alphasbiv <- .alphabivariate(cgh.em[c(1:(2 * nosamp))], nosamp, 2)
	alphas <- .alphaest(cgh.em[c(1:(2 * nosamp))], nosamp, 2)
	cgh.em <- cbind(matrix(cgh.em[c(1:(2 * nosamp))], ncol = 2, byrow = TRUE), cgh.em[c((2 * nosamp + 1):((2 + 1) * nosamp))])
	c1 <- (alphasbiv[1]/alphasbiv[3] - alphas[1]/alphas[2])^(-1)
	c2 <- (alphasbiv[2]/alphasbiv[3] - alphas[2]/alphas[1])^(-1)
	mu1 <- (1/nosamp) * sum(cgh.em[, 3] * (cgh.em[,1] * alphas[2] - alphasbiv[3])/(alphasbiv[1] * alphas[2] - alphas[1] * alphasbiv[3]))
	mu2 <- (1/nosamp) * sum(cgh.em[, 3] * (cgh.em[,2] * alphas[1] - alphasbiv[3])/(alphasbiv[2] * alphas[1] - alphas[2] * alphasbiv[3]))
	overall.data <- cbind(c(1:2), cbind(alphas, c(mu1, mu2)))
	clone.data <- cbind(rep(c(1:2), nosamp), datafortest[1, c(1:(2 * nosamp))], as.numeric(t(matrix(rep(datafortest[1, c((2 * nosamp + 1):((2 + 1) * nosamp))], 2), nrow = nosamp))))
    
	symbols(overall.data[, 1], overall.data[, 3], circles = overall.data[, 2]/3, fg = "red", bg = "indianred1", ylim = c(min(clone.data[, 
        	3]) - 3/4, max(clone.data[, 3]) + 3/4), xlim = c(1/2, 2 + 1/2), ylab = "gene expression", xlab = "DNA copy number", xaxt = "n", inches = FALSE, lwd = 3, main = rownames(datafortest)[1])
	    if (lossorgain[1] == 1) {
        	axis(1, at = c(1:2), labels = c("loss", "no-loss"))
	} else {
        	axis(1, at = c(1:2), labels = c("no-gain", "gain"))
	}
	symbols(clone.data[, 1], clone.data[, 3], circles = clone.data[,2]/3, fg = "blue", add = TRUE, inches = FALSE)
	return(invisible(NULL))
}




cisEffectTable <- function (testRes, number = 10, sort.by = NULL){
	################################################################################################
	# returns top results
	################################################################################################

	# input checks
	if (as.character(class(testRes)) != "cisTest"){ stop("testRes not of class cisTest.") }

	# determine top
	if (is.null(sort.by)){ 
		top <- 1:number 
	} else {
		if (sort.by == "p.value"){ 
			top <- which(rank(testRes@p.value, ties.method="first") <= number) 
			top <- top[rank(testRes@p.value[top], ties.method="first")]
		}
		if (sort.by == "R2"){ 
			top <- which(rank(testRes@R2, ties.method="first") > length(testRes@R2) - number) 
			top <- top[number - rank(testRes@R2[top], ties.method="first") + 1]
		}
		if (sort.by == "effect"){ 
			top <- which(rank(testRes@effectSize, ties.method="first") > length(testRes@R2) - number) 
			top <- top[number - rank(testRes@effectSize[top], ties.method="first") + 1]
		}
	}
	
	# construct table
	if (testRes@analysisType == "univariate"){
		if (dim(testRes@geneInfo)[2] > 1){
			tab <- data.frame(testRes@geneInfo[top,], geneId = testRes@geneId[top],
				comparison = testRes@comparison[top], av.prob1 = testRes@av.prob1[top],
				av.prob2 = testRes@av.prob2[top], effectSize = testRes@effectSize[top],
				R2 = testRes@R2[top], 	p.value = testRes@p.value[top], adjP.value = testRes@adjP.value[top])
		} else {
			tab <- data.frame(testRes@geneInfo[top], geneId = testRes@geneId[top],
				comparison = testRes@comparison[top], av.prob1 = testRes@av.prob1[top],
				av.prob2 = testRes@av.prob2[top], effectSize = testRes@effectSize[top],
				R2 = testRes@R2[top], 	p.value = testRes@p.value[top], adjP.value = testRes@adjP.value[top])
		}
	}
	if (testRes@analysisType == "regional"){
		if (dim(testRes@geneInfo)[2] > 1){
			tab <- data.frame(testRes@geneInfo[top,], geneId = testRes@geneId[top],
				comparison = testRes@comparison[top], av.prob1 = testRes@av.prob1[top],
				av.prob2 = testRes@av.prob2[top], effectSize = testRes@effectSize[top],
				R2 = testRes@R2[top], regId = testRes@regId[top], beginReg = testRes@beginReg[top],
				endReg = testRes@endReg[top], shrinkage = testRes@shrinkage[top],
				p.value = testRes@p.value[top], adjP.value = testRes@adjP.value[top])
		} else {
			tab <- data.frame(testRes@geneInfo[top,], geneId = testRes@geneId[top],
				comparison = testRes@comparison[top], av.prob1 = testRes@av.prob1[top],
				av.prob2 = testRes@av.prob2[top], effectSize = testRes@effectSize[top],
				R2 = testRes@R2[top], regId = testRes@regId[top], beginReg = testRes@beginReg[top],
				endReg = testRes@endReg[top], shrinkage = testRes@shrinkage[top],
				p.value = testRes@p.value[top], adjP.value = testRes@adjP.value[top])
		}
	}
	return(tab)
}

cisEffectTest <- function(CNdata, GEdata, genes2test=NULL, GEchr, analysisType="univariate", testStatistic="wcvm", nPerm=10000, lowCiThres=0.10, verbose=TRUE){
	################################################################################################
	# wrapper function for testing
	################################################################################################

	# input checks
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	if ( !(as.character(class(GEdata)) == "ExpressionSet" | as.character(class(GEdata)) == "eSet") ){ stop("GEdata not of class ExpressionSet.") }

	# actual start testing	
	if (verbose){ cat("formatting data for testing...", "\n") }
	CNdata <- cghCall2subset(CNdata, genes2test, verbose)
	GEdata <- ExpressionSet2subset(GEdata, genes2test, verbose)

	# extract data from objects
	nosamp <- ncol(exprs(GEdata))
	cghdata.probs <- numeric()
	for (i in 1:dim(calls(CNdata))[2]){
		if (is.null(probamp(CNdata)[,i])){
			cghdata.probs <- cbind(cghdata.probs, cbind(probloss(CNdata)[,i], probnorm(CNdata)[,i], probgain(CNdata)[,i]))
		} else {
			cghdata.probs <- cbind(cghdata.probs, cbind(probloss(CNdata)[,i], probnorm(CNdata)[,i], probgain(CNdata)[,i] + probamp(CNdata)[,i]))
		}
	}
	nclass <- dim(cghdata.probs)[2] / dim(calls(CNdata))[2]
	annInfo <- fData(GEdata)
	chrInfo <- annInfo[, GEchr]

	# Estimate alpha parameters per clone, a0 classes
	if (verbose){ cat("pre-test...", "\n") }
	alphamat <- t(apply(cghdata.probs, 1, .alphaest, nosamp=nosamp, a=nclass))

	# Perform pre-test and merge columns
	datacgh2 <- as.matrix(t(apply(cbind(alphamat, cghdata.probs), 1, .pretest)))
	lossorgain <- datacgh2[,1]
	alphas2 <- datacgh2[,2:3]
	datafortest <- as.matrix(cbind(datacgh2, exprs(GEdata))[,-(1:3)])

	# Estimate new alpha and bivariate alpha parameters per clone
	alphabivmat <- t(apply(datafortest, 1, .alphabivariate, nosamp=nosamp, a=2))
	colnames(alphabivmat) <- c("a11", "a22", "a12")

	rm(GEdata, CNdata)
	gc()

	# actual start testing
	results <- NULL
	for (chr in 1:length(unique(chrInfo))){
		chrs.present <- unique(chrInfo)
		set.seed(7396)
		genes.per.time <- which(chrInfo == chrs.present[chr])
		data.both <- datafortest[genes.per.time, , drop=FALSE]
		if (dim(data.both)[1] == 1){
			no.genes.only.one <- TRUE
			data.both <- rbind(data.both, data.both)
		} else {
			no.genes.only.one <- FALSE
		}

		if (verbose){  cat(paste("chromosome ",  chrs.present[chr], " started...", sep=""), "\n") }
		if (verbose){  cat(paste(length(genes.per.time), " genes to be tested...", sep=""), "\n") }

		# "Robustly" estimate effect sizes. These are used for generating empirical distribution of shift parameters
		alleffects <- sapply(1:nrow(data.both), .shift.est, data2=data.both, alphabmat=alphabivmat[genes.per.time, , drop=FALSE], alphanmat=alphas2[genes.per.time, , drop=FALSE], nosamp=nosamp, a=2, minalphathr=10) 

		R2 <- .R2.stat(data.both, 2, nosamp)
		R2[is.na(alleffects)] <- NA
	
		if (analysisType == "univariate"){
			# Univariate analysis, calculates the p-values, raw and adjusted.
			if (testStatistic == "wcvm"){
				adjpvals <- .uni.an.wcvm(data.both, nosamp, 2, nPerm, low.ci.thres=lowCiThres, verbose=verbose)
			}
			if (testStatistic == "wmw"){
				adjpvals <- .uni.an.prob(data.both, nosamp, 2, nPerm, low.ci.thres=lowCiThres, verbose=verbose)
			}
			# REFORMATTING FOR RESULTS PRESENTATION   
			results.temp <- data.frame()

			if (no.genes.only.one){
				results.temp <- c(genes2test[genes.per.time], 
						lossorgain[genes.per.time], 
						round(alphas2[genes.per.time,], digits=4), 
						round(alleffects, digits=4), 
						round(R2[1], digits=4), round(adjpvals[1,2:3], digits=4))
				names(results.temp)[1:6] <- c("gene.id", "comparison", "av.probs.1", "av.probs.2", "effect.size", "R2")
				results <- rbind(results, c(annInfo[genes.per.time,], results.temp))
				rownames(results)[nrow(results)] <- rownames(annInfo)[genes.per.time]
			} else {
				results.temp <- cbind(genes2test[genes.per.time], 
						lossorgain[genes.per.time], 
						round(alphas2[genes.per.time,], digits=4), 
						round(alleffects, digits=4), 
						round(R2, digits=4), round(adjpvals[,2:3], digits=4))
				colnames(results.temp)[1:6] <- c("gene.id", "comparison", "av.probs.1", "av.probs.2", "effect.size", "R2")
				results <- rbind(results, cbind(annInfo[genes.per.time,], results.temp))
			}
		}

		if (analysisType == "regional"){
			# Regional analysis, calculates the p-values, raw and adjusted.
			if (testStatistic == "wcvm"){
				adjpvals <- .shrin.an.wcvm(data.both, nosamp, 2, nPerm, low.ci.thres=lowCiThres, datacgh.org=data.both[,c(1:(2*nosamp))], verbose=verbose)
			}
			if (testStatistic == "wmw"){
				adjpvals <- .shrin.an.prob(data.both, nosamp, 2, nPerm, low.ci.thres=lowCiThres, datacgh.org=data.both[,c(1:(2*nosamp))], verbose=verbose)
			}

			# REFORMATTING FOR RESULTS PRESENTATION
			results.temp <- data.frame()   

			# "Robustly" estimate effect sizes. These are used for generating empirical distribution of shift parameters
			alleffects <- sapply(1:nrow(data.both), .shift.est, data2=data.both, alphabmat=alphabivmat[genes.per.time, , drop=FALSE], alphanmat=alphas2[genes.per.time, , drop=FALSE], nosamp=nosamp, a=2, minalphathr=10) 

			R2 <- .R2.stat(data.both, 2, nosamp)
			R2[is.na(alleffects)] <- NA
	
			if (no.genes.only.one){
			results.temp <- c(genes2test[genes.per.time], 
					lossorgain[genes.per.time],
					round(alphas2[genes.per.time,], digits=4), 
					round(alleffects, digits=4),
					round(R2[1], digits=4), adjpvals[1,-1])
				names(results.temp)[1:6] <- c("gene.id", "comparison", "av.probs.1", "av.probs.2", "effect.size", "R2")
				results.temp[12] <- round(results.temp[12], digits=4)
				results.temp[11] <- round(results.temp[11], digits=4)
				results.temp[9] <- results.temp[8]
				results <- rbind(results, c(annInfo[genes.per.time,], results.temp))
				rownames(results)[nrow(results)] <- rownames(annInfo)[genes.per.time]
			} else {
				results.temp <- cbind(genes2test[genes.per.time], 
						lossorgain[genes.per.time],
						round(alphas2[genes.per.time,], digits=4), 
						round(alleffects, digits=4),
						round(R2, digits=4), adjpvals[,-1])
				colnames(results.temp)[1:6] <- c("gene.id", "comparison", "av.probs.1", "av.probs.2", "effect.size", "R2")
				results.temp[,12] <- round(results.temp[,12], digits=4)
				results.temp[,11] <- round(results.temp[,11], digits=4)
				results <- rbind(results, cbind(annInfo[genes.per.time,], results.temp))
			}
	
		}
		if (verbose){  cat(paste("chromosome ",  chrs.present[chr], " done...", sep=""), "\n") }

	}
	if (is.list(results[,dim(results)[2]-1])){
		results[,dim(results)[2]] <- p.adjust(unlist(results[,dim(results)[2]-1]), "BH")
	} else {
		results[,dim(results)[2]] <- p.adjust(results[,dim(results)[2]-1], "BH")
	}

	if (verbose){  cat("ready: testing done", "\n") }

	# format output
	if (analysisType == "univariate"){
		testRes <- new("cisTest", geneInfo = results[,1:(ncol(results)-8)], geneId = results[,ncol(results)-7],
				comparison = results[,ncol(results)-6], av.prob1 = results[,ncol(results)-5],
				av.prob2 = results[,ncol(results)-4], effectSize = results[,ncol(results)-3],
				R2 = results[,ncol(results)-2], p.value = results[,ncol(results)-1], adjP.value = results[,ncol(results)],
				testStatistic = testStatistic, analysisType = analysisType, nPerm = nPerm)
	}
	if (analysisType == "regional"){
		testRes <- new("cisTest", geneInfo = results[,1:(ncol(results)-12)], geneId = results[,ncol(results)-11],
				comparison = results[,ncol(results)-10], av.prob1 = results[,ncol(results)-9],
				av.prob2 = results[,ncol(results)-8], effectSize = results[,ncol(results)-7],
				R2 = results[,ncol(results)-6], 
				regId = results[,ncol(results)-5], beginReg = results[,ncol(results)-4],
				endReg = results[,ncol(results)-3], shrinkage = results[,ncol(results)-2],
				p.value = results[,ncol(results)-1], adjP.value = results[,ncol(results)],
				testStatistic = testStatistic, analysisType = analysisType, nPerm = nPerm)
	}

	return(testRes)
}



cisEffectTune <- function(CNdata, GEdata, testStatistic, nGenes=250, nPerm=250, minCallProbMass=0.10, verbose=TRUE){
	######################################################################################################
	# function that filters genes that have a low probability of becoming significant
	# and determines which test is performed: loss-vs-no.loss or no.gain-vs-gain.
	######################################################################################################

	# check for presence missing values.
	if (sum(is.na(exprs(GEdata))) > 0){ 
		stop("Gene expression matrix contains missing values: not allowed.")
	}

	# check for equal number of samples.
	if (dim(exprs(GEdata))[2] != dim(copynumber(CNdata))[2]){ 
		stop("Gene expression and copy number matrices contain unequal number of samples.")
	}

	# check for equal number of probes.
	if (dim(exprs(GEdata))[1] != dim(copynumber(CNdata))[1]){ 
		stop("Gene expression and copy number matrices contain unequal number of samples: impossible after matching.")
	}

	set.seed(7396)

	# extract data from objects
	data <- list()
	data$ann <- fData(GEdata)
	data$em <- exprs(GEdata)
	nosamp <- ncol(exprs(GEdata))
	cghdata.probs <- numeric()
	for (i in 1:dim(calls(CNdata))[2]){
		if (is.null(probamp(CNdata)[,i])){
			cghdata.probs <- cbind(cghdata.probs, cbind(probloss(CNdata)[,i], probnorm(CNdata)[,i], probgain(CNdata)[,i]))
		} else {
			cghdata.probs <- cbind(cghdata.probs, cbind(probloss(CNdata)[,i], probnorm(CNdata)[,i], probgain(CNdata)[,i] + probamp(CNdata)[,i]))
		}
	}
	data$cgh <- cghdata.probs
	nclass <- dim(cghdata.probs)[2] / dim(calls(CNdata))[2]

	# Estimate alpha parameters per clone, a0 classes
	if (verbose){ cat("pre-test...", "\n") }
	alphamat <- t(apply(data$cgh, 1, .alphaest, nosamp=nosamp, a=nclass))

	# Perform pre-test and merge columns
	datacgh2 <- as.matrix(t(apply(cbind(alphamat, data$cgh), 1, .pretest)))
	lossorgain <- datacgh2[,1]
	alphas2 <- datacgh2[,2:3]
	data.both <- as.matrix(cbind(datacgh2, data$em)[,-(1:3)])

	# remove features with too little spread in call probs mass
	suffCPM <- which(apply(alphas2, 1, min) > minCallProbMass)
	data$ann <- data$ann[suffCPM,]
	data$cgh <- data$cgh[suffCPM,]
	data$em <- data$em[suffCPM,]
	alphamat <- alphamat[suffCPM,]
	datacgh2 <- datacgh2[suffCPM,]
	lossorgain <- lossorgain[suffCPM]
	alphas2 <- alphas2[suffCPM,]
	data.both <- data.both[suffCPM,]
	init.prop.kept <- length(suffCPM)/dim(exprs(GEdata))[1]
	
	# Estimate new alpha and bivariate alpha parameters per clone
	# alphanewmat <- alphas2
	alphabivmat <- t(apply(data.both, 1, .alphabivariate, nosamp=nosamp, a=2))
	colnames(alphabivmat) <- c("a11", "a22", "a12")

	# start of tuning
	if (verbose){ cat("tuning started...", "\n") }

	# Calculate the measure of unbalance and estimate the effect sizes
	# UNBALANCE, SHIFT ESTIMATES AND TUNING, needs to be done AFTER null distr computation, BEFORE final p-value and fdr computation.
	powerunbal <- sapply(1:nrow(data.both), .powerunbalance, alphabmat=alphabivmat)
	allest <- sapply(1:nrow(data.both), .shift.est, data2=data.both, alphabmat=alphabivmat, alphanmat=alphas2, nosamp=nosamp, a=2, minalphathr=10)

	# "Robustly" estimate effect sizes. These are used for generating empirical distribution of shift parameters
	estrobust <- sapply(1:nrow(data.both), .shift.est, data2=data.both, alphabmat=alphabivmat, alphanmat=alphas2, nosamp=nosamp, a=2, minalphathr=3) 

	# Draw histogram of robust effect sizes
	shiftsam <- estrobust[!is.na(estrobust)]

	# Determine null distribution for tuning genes.
	seqgenes <- floor(seq(1,nrow(data.both), length.out=nGenes))
	if (testStatistic == "wcvm") {
		null.dists.tune <- .nulldist.all.wcvm(data.both[seqgenes,], nosamp, 2, nperm=nPerm, verbose)
	} else {
		null.dists.tune <- .nulldist.all.wmw(data.both[seqgenes,], nosamp, 2, nperm=nPerm, verbose)
	}   

	# Determine the proportion of diff exp genes on small data set
	pi0 <- .pi0est(data.in=data.both, null.dists.tune=null.dists.tune, seqg=seqgenes, test.stat=testStatistic, nperm=nPerm, a=2, nosamp=nosamp) 

	# Do actual tuning, returns list of genes that are propagated into the test
	data.tuned <- .tuning(data.in=data.both, datacgh.in=datacgh2, allest=allest, alphanewmat=alphas2, powerunbal=powerunbal, null.dists.tune=null.dists.tune, seqg=seqgenes, test.stat=testStatistic, nperm=nPerm, pi0=pi0, fdrcut=0.05, nresamp=100, gridnr=200, minim=10, a=2, nosamp=nosamp, shiftsam=shiftsam, verbose)
	genes2test <- .datareduce(powerunbal, as.real(data.tuned))
	if (verbose){ cat("ready: tuning done", "\n") }

	# return(data.tuned)
	return(as.numeric(suffCPM[genes2test]))
}
