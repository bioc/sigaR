profilesPlot <- function(CNdata, GEdata, sampleNo, chr=0, verbose=TRUE){
	################################################################################################
	# function for plotting the copy number and gene expression profiles of an individual sample
	################################################################################################

	profilePlot <- function (x, segment, z){

		chrom           <- CGHbase::chromosomes(x)
		chrom.labels    <- unique(chrom)
		nclone          <- length(chrom)
	
		genomdat        <- CGHbase::copynumber(x)[,1]
		probsdraw       <- cbind(CGHbase::probloss(x)[,1], CGHbase::probnorm(x)[,1], CGHbase::probgain(x))
		if (!is.null(probamp(x))){ probsdraw <- cbind(probsdraw, probamp(x)[,1]) }
        
		widths          <- segment[,3] - segment[,2] + 1
		plot.data       <- probsdraw[segment[,2],,drop=FALSE]
    
		par(mar=c(5, 4, 4, 4) + 0.2)
        
		### Plot the probability bars
		if (z){ barplot(t(plot.data), width=widths, border=FALSE, space=0, col=c("#FF000033", "white", "#00FF0033"), las=1, cex.axis=1, cex.lab=1, xaxt="n")
		} else { barplot(t(plot.data), width=widths, border=FALSE, space=0, col=c("red","white","green"), las=1, cex.axis=1, cex.lab=1, xaxt="n") 	}
        
		lim <- par("usr")
		if (z) { YaxisLimits <- c(floor(min(copynumber(x)[,1], na.rm=TRUE)), ceiling(max(copynumber(x)[,1], na.rm=TRUE))); lim[3:4] <- YaxisLimits
		} else { lim[3:4] <- c(-5, 5)	}
		par(usr=lim)
		if (z) { dticks <- seq(YaxisLimits[1], YaxisLimits[2], by=1)
		} else { dticks <- seq(-5, 5, by=1) 	}
		axis(4, at=dticks, labels=dticks, srt=270, las=1, cex.axis=1, cex.lab=1)
		box()
        
		### Add axis labels
		if (z) { mtext("norm. GE", side=4, line=3, srt=270) } else { mtext("norm. CN", side=4, line=3, srt=270) }
		mtext("probability", side=2, line=3, srt=270)
        
		#### add vert lines at chromosome ends
		abline(h=0) 
		if (z) { for (iii in 1:length(cumsum(table(chrom)))) { segments(cumsum(table(chrom))[[iii]], YaxisLimits[1], cumsum(table(chrom))[[iii]], YaxisLimits[2], lty=2) } 
		} else { for (iii in 1:length(cumsum(table(chrom)))) { segments(cumsum(table(chrom))[[iii]], -5, cumsum(table(chrom))[[iii]], 5, lty=2) }	}

		if (z) { title(paste("GE profile:", sampleNames(x)[1], sep=" "))
		} else { title(paste("CN profile:", sampleNames(x)[1], sep=" ")) }

        	### Add log2ratios
		points((1:nclone)-.5,genomdat,cex=.1)
        
		### X-axis with chromosome labels
		ax <- (cumsum(table(chrom))+c(0, cumsum(table(chrom))[-length(cumsum(table(chrom)))]))/2
		axis(side=1, at=ax, labels=chrom.labels, cex=.2, lwd=.5, las=1, cex.axis=1, cex.lab=1) # bottom axis
            
		### Blue lines for segment means
		if (z) { for (jjj in (1:nrow(segment))){ segments(segment[jjj,2], segment[jjj,1], segment[jjj,3], segment[jjj,1], col="blue", lwd=3) }
		} else { for (jjj in (1:nrow(segment))){ segments(segment[jjj,2], segment[jjj,1], segment[jjj,3], segment[jjj,1], col="blue", lwd=3) } }		
	}

	# input checks
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	if ( !(as.character(class(GEdata)) == "ExpressionSet" | as.character(class(GEdata)) == "eSet") ){ stop("GEdata not of class ExpressionSet.") }
	if (dim(fData(CNdata))[1] != dim(fData(GEdata))[1]){ stop("CN and GE data have different number of rows.") }
	if (!all(fData(CNdata)[,1] == fData(GEdata)[,1])){ stop("chrosome annotation between CN and GE does not match.") }
	if (!(as.character(class(sampleNo))=="numeric" | as.character(class(sampleNo))=="integer" )){ stop("sampleNo of wrong class.") }
	if (!(as.character(class(chr))=="numeric" | as.character(class(chr))=="integer" )){ stop("sampleNo of wrong class.") }
	if ( !(sampleNo >= 1 | sampleNo <= ncol(exprs(GEdata)) ) ){ stop("sampleNo of wrong length.") }
	if ( !(chr >= 0 | chr <= 23) ){ stop("sampleNo of wrong length.") }
	if ( !(is.logical(verbose)) ){ stop("verbose wrongly specified.") }

	# filter data corresponding to chromosome 
	if (chr != 0){
		ids <- which(fData(CNdata)[,1] == chr)
		if (length(ids) == 0){ stop("supplied chromosome does not exist.") }
		CNdata <- cghCall2subset(CNdata, ids, verbose=verbose)
		GEdata <- ExpressionSet2subset(GEdata, ids, verbose=verbose)
	}

	# calculate median expression per segment
	SegExpr <- numeric()
	# SegData <- segmented(CNdata[,sampleNo])
	# segmentsCN <- sigaR:::.makeSegments(segmented(CNdata[,sampleNo]))	
	segmentsCN <- .makeSegments(CGHbase::segmented(CNdata[,sampleNo]))	
	segmentsGE <- segmentsCN
	for (j in 1:dim(segmentsCN)[1]){
		ids <- c(segmentsCN[j,2]:segmentsCN[j,3])
		medSegExpr <- median(Biobase::exprs(GEdata)[ids, sampleNo])
		segmentsGE[j, 1] <- medSegExpr
		SegExpr <- c(SegExpr, rep(medSegExpr, length(ids)))
	}

	# modify CGHcall object
	CNdataGE <- CNdata
	segmented(CNdataGE)[,sampleNo] <- SegExpr
	copynumber(CNdataGE)[,sampleNo] <- Biobase::exprs(GEdata)[,sampleNo]

	# plot profiles
	op <- par(mfrow = c(2, 1), pty = "m")
	profilePlot(CNdataGE[,sampleNo], segment=segmentsGE, z=TRUE)	
	profilePlot(CNdata[,sampleNo], segment=segmentsCN, z=FALSE)
	par(op)
	return(invisible(NULL))
}


CNGEheatmaps <- function(CNdata, GEdata, location="mode", colorbreaks="equiquantiles"){

	makeRegionsSigaR <- function(CNsegData){
		# determine regions
		splitter <- list()
		splitter[[1]] <- c(1)
		index.temp <- 1
		j <- 1
		for (i in 1:(dim(CNsegData)[1]-1)){
			if (all(CNsegData[i,] == CNsegData[i+1,])){
    				index.temp <- c(index.temp,i+1)
	           		splitter[[j]] <- index.temp
			} else {
				index.temp <- i+1
	     			j <- j + 1
      	     		splitter[[j]] <- index.temp
			}
		}
		regDetails <- NULL
		for (i in 1:length(splitter)){
			regDetails <- rbind(regDetails, c(min(splitter[[i]]),max(splitter[[i]])))
		}
		return(regDetails)
	}

	makeColorScheme <- function(xTemp, colorbreaks){
		# make color palette for heatmap from matrix X
		# xTemp <- as.numeric(X)
		histres <- hist(xTemp, plot=FALSE, n=511)
		if (location == "median"){ xMode <- median(xTemp) }
		if (location == "mean"){ xMode <- mean(xTemp) }
		if (location == "mode"){ xMode <- histres$mids[which.max(histres$density)] }
		# xMode <- histres$mids[which.max(histres$density)]
		xTempBelowMode <- xTemp[xTemp < xMode]
		xTempAboveMode <- xTemp[xTemp >= xMode]
		xTempBelowMode <- cbind(xTempBelowMode, ecdf(xTempBelowMode)(xTempBelowMode))[order(xTempBelowMode),]
		xTempAboveMode <- cbind(xTempAboveMode, ecdf(xTempAboveMode)(xTempAboveMode))[order(xTempAboveMode),]
		if(any(xTempBelowMode[,2]==0)){ xTempBelowMode[any(xTempBelowMode[,2]==0),2] <- 10^(-10) }
		if(any(xTempBelowMode[,2]==1)){ xTempBelowMode[any(xTempBelowMode[,2]==1),2] <- 1-10^(-10) }
		if(any(xTempAboveMode[,2]==0)){ xTempAboveMode[any(xTempAboveMode[,2]==0),2] <- 10^(-10) }
		if(any(xTempAboveMode[,2]==1)){ xTempAboveMode[any(xTempAboveMode[,2]==1),2] <- 1-10^(-10) }
		if (colorbreaks == "equiquantiles"){
			histresBM <- hist(xTempBelowMode[,2], plot=FALSE, n=100)
			histresAM <- hist(xTempAboveMode[,2], plot=FALSE, n=101)
			breaksX <- c(quantile(xTempBelowMode[,1], probs=histresBM$breaks, na.rm=TRUE), xMode, quantile(xTempAboveMode[,1], probs=histresAM$breaks, na.rm=TRUE))
			collistX <- c(maPalette(low = "red", high="black", k=length(histresBM$breaks)), maPalette(low="black", high="green", k=length(histresAM$breaks)))
		}
		if (colorbreaks == "equidistant"){
			collistBelowMode <- unique(maPalette(low = "red", high="black", k=100))
			collistAboveMode <- unique(maPalette(low = "black", high="green", k=100))
			breaksX <- c(seq(min(xTemp), xMode, length.out=length(collistBelowMode)+1), seq(xMode, max(xTemp), length.out=length(collistAboveMode))[-1])
			collistX <- unique(c(collistBelowMode, collistAboveMode))
		}
		return(list(collist=collistX, breaks=breaksX))
	}

	# input checks
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	if ( !(as.character(class(GEdata)) == "ExpressionSet" | as.character(class(GEdata)) == "eSet") ){ stop("GEdata not of class ExpressionSet.") }
	if (dim(fData(CNdata))[1] != dim(fData(GEdata))[1]){ stop("CN and GE data have different number of rows.") }
	if (!all(fData(CNdata)[,1] == fData(GEdata)[,1])){ stop("chrosome annotation between CN and GE does not match.") }
	if (!(location %in% c("mode", "median", "mean"))){ stop("location parameter ill-specified.") }
	if (!(colorbreaks %in% c("equidistant", "equiquantiles"))){ stop("colorbreaks parameter ill-specified.") }

	# determine regions
	regDetails <- makeRegionsSigaR(CGHbase::segmented(CNdata))
	regChr <- fData(CNdata)[regDetails[ ,1], 1]

	# calculate and extract segment level data for expression and copy number (respectively)
	SegData <- list()
	for (sampleNo in 1:ncol(GEdata)){
		# calculate median expression per segment
		SegExpr <- numeric()
		# segments <- sigaR:::.makeSegments(segmented(CNdata[,sampleNo]))	
		segments <- .makeSegments(segmented(CNdata[,sampleNo]))	
		segments <- segments[,c(2,3,1), drop=FALSE]
		for (j in 1:dim(segments)[1]){
			ids <- c(segments[j,1]:segments[j,2])
			medSegExpr <- median(exprs(GEdata)[ids, sampleNo])
			SegExpr <- c(SegExpr, medSegExpr)
		}
		SegData[[sampleNo]] <- cbind(segments, SegExpr)
	}
	rm(GEdata, CNdata)
	gc(verbose=FALSE)	

	# expand segment level data to region level data
	regSegLog2 <- matrix(NA, nrow=nrow(regDetails), ncol=length(SegData))
	regSegExprs <- matrix(NA, nrow=nrow(regDetails), ncol=length(SegData))
	for (i in 1:length(SegData)){
		reps <- apply(SegData[[i]][,1:2, drop=FALSE], 1, function(Z, Y){ length(intersect(which(Y[,1] >= Z[1]), which(Y[,2] <= Z[2]))) }, Y=regDetails)
		regSegLog2[,i] <- unlist(apply(cbind(SegData[[i]][,3], reps), 1, function(Z){ rep(Z[1], Z[2]) }))
		regSegExprs[,i] <- unlist(apply(cbind(SegData[[i]][,4], reps), 1, function(Z){ rep(Z[1], Z[2]) }))
	}

	# make color palette for heatmaps
	colorSchemeExprs <- makeColorScheme(unique(as.numeric(regSegExprs)), colorbreaks)
	colorSchemeCopynum <- makeColorScheme(unique(as.numeric(regSegLog2)), colorbreaks)

	# preparation of the plotting of the heatmap
	# generate alternating colors for chromosomes.
	chrInd <- rep(0, length(regChr))
	chrInd[(regChr %% 2 == 0)] <- 1
	chrColor <- rep("blue", length(regChr))
	chrColor[(regChr %% 2 == 0)] <- c("yellow")

	# generate labels for begin points of chromosomes.
	Y <- rep(FALSE, length(regChr))
	for (i in 2:length(regChr)){
		if ((regChr[i-1] != regChr[i])){ Y[i] <- TRUE }
	}
	Y[1] <- TRUE
	beginChr <- rep("", length(regChr))
	beginChr[Y] <- regChr[Y]

	# start plotting
	def.par <- par
	fl <- layout(matrix(c(1,2,3,1,2,3,1,2,3), 3, 3, byrow = TRUE), widths=c(1,9,9))
	# layout.show(fl)
	par(mar=c(3,2,4,0))
	image(z=matrix(chrInd, nrow=1), xaxt="n", yaxt="n", col=c("blue", "yellow"))
	axis(2, at=(which(Y)-1)/(length(Y)-1), labels=regChr[Y], tick=FALSE, las=1)
	par(mar=c(3,1,4,1))
	image(z=t(regSegLog2), xaxt="n", yaxt="n", col=colorSchemeCopynum$collist, breaks=colorSchemeCopynum$breaks, main="copy number data")
	par(mar=c(3,1,4,1))
	image(z=t(regSegExprs), xaxt="n", yaxt="n", col=colorSchemeExprs$collist, breaks=colorSchemeExprs$breaks, main="gene expression data")
	par(def.par)
	return(invisible(NULL))
}

