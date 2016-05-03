

cghCall2maximumSubset <- function(CNdata, featuresAndWeights, chr, bpstart, bpend, ncpus=1, verbose=TRUE){
   	############################################################################
	# function shrinks cghCall-object to a subset of the features
	############################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	# if (!(as.character(class(featureSubset))=="numeric" | as.character(class(featureSubset))=="integer" )){ stop("featureSubset of wrong class.") }
	# if (!(as.character(class(weights))=="numeric" | as.character(class(weights))=="integer" )){ stop("weights of wrong class.") }
	# if ( !(length(featureSubset) >= 1) ){ stop("featureSubset contains no row numbers.") }
	# if ( !(min(featureSubset) >= 1 | max(featureSubset) <= nrow(fData(CNdata))) ){ stop("featureSubset contains illegal row numbers.") }
	# if ( length(featureSubset) != length(weights) ){ stop("featureSubset and weights of unequal length.") }

	if (!(as.character(class(chr))=="numeric" | as.character(class(chr))=="integer" )){ stop("chr of wrong class.") }
	if (length(chr) != 1){ stop("chr of wrong length.") }
	if ( !(chr >= 1 | chr <= ncol(fData(CNdata))) ){ stop("chr exceeds number of columns in featureData.") }
	if (!(as.character(class(bpstart))=="numeric" | as.character(class(bpstart))=="integer" )){ stop("bpstart of wrong class.") }
	if (length(bpstart) != 1){ stop("bpstart of wrong length.") }
	if ( !(bpstart >= 1 | bpstart <= ncol(fData(CNdata))) ){ stop("bpstart exceeds number of columns in featureData.") }
	if (!(as.character(class(bpend))=="numeric" | as.character(class(bpend))=="integer" )){ stop("bpend of wrong class.") }
	if (length(bpend) != 1){ stop("bpend of wrong length.") }
	if ( !(bpend >= 1 | chr <= ncol(fData(CNdata))) ){ stop("bpend exceeds number of columns in featureData.") }

	# disecting the cghCall-object
	amplifications <- FALSE
	copynumberMat <- CGHbase::copynumber(CNdata)
	segmentedMat <- CGHbase::segmented(CNdata)
	problossMat <- CGHbase::probloss(CNdata)
	probnormMat <- CGHbase::probnorm(CNdata)
	probgainMat <- CGHbase::probgain(CNdata)
	if (!is.null(CGHbase::probamp(CNdata))){ probampMat <- CGHbase::probamp(CNdata); amplifications <- TRUE }
	if (!is.null(CGHbase::probdloss(CNdata))){ probdlossMat <- CGHbase::probdloss(CNdata); doubleLoss <- TRUE }
	CNann <- fData(CNdata)
	rm(CNdata)
	gc()

	# identify maximum deviating features
	if (verbose){ cat("identify maximum deviating features ...", "\n") }
	if (ncpus == 1){ idMat <- t(sapply(featuresAndWeights, .slhFuncMaxCGHcallData, segmentedMat, simplify=TRUE)) }
	if (ncpus > 1){ 
		sfInit(parallel=TRUE, cpus=ncpus)
		sfLibrary(sigaR, verbose=FALSE)
		idMat <- t(sfSapply(featuresAndWeights, .slhFuncMaxCGHcallData, segmentedMat, simplify=TRUE)) 
	}

	# subsetting of segmented data
	if (verbose){ cat("generate subsetted matrices ...", "\n") }
	if (ncpus == 1){ copynumberMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, copynumberMat, simplify=TRUE)) }
	if (ncpus > 1){ copynumberMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, copynumberMat, simplify=TRUE)) }
	if (ncpus == 1){ segmentedMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, segmentedMat, simplify=TRUE)) }
	if (ncpus > 1){ segmentedMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, segmentedMat, simplify=TRUE)) }
	if (ncpus == 1){ problossMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, problossMat, simplify=TRUE)) }
	if (ncpus > 1){ problossMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, problossMat, simplify=TRUE)) }
	if (ncpus == 1){ probnormMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probnormMat, simplify=TRUE)) }
	if (ncpus > 1){ probnormMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probnormMat, simplify=TRUE)) }
	if (ncpus == 1){ probgainMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probgainMat, simplify=TRUE)) }
	if (ncpus > 1){ probgainMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probgainMat, simplify=TRUE)) }
	if (amplifications){ 
		if (ncpus == 1){ probampMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probampMat, simplify=TRUE)) }
		if (ncpus > 1){ probampMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probampMat, simplify=TRUE)) }
	}
	if (doubleLoss){ 
		if (ncpus == 1){ probdlossMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probdlossMat, simplify=TRUE)) }
		if (ncpus > 1){ probdlossMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probdlossMat, simplify=TRUE)) }
	}

	# subsetting of call data
	if (verbose){ cat("generate subsetted call matrices ...", "\n")	}
	callsMat <- numeric()
	for (i in 1:ncol(problossMat)){
		if (!amplifications & !doubleLoss){
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 2) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 2) }
		} 
		if (amplifications & !doubleLoss){
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 2) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 2) }
		}
		if (!amplifications & doubleLoss){
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(probdlossMat[,i], problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 3) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(probdlossMat[,i], problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 3) }
		}
		if (amplifications & doubleLoss){
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(probdlossMat[,i], problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 3) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(probdlossMat[,i], problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 3) }
		}
	}	

	# weighted subsetting of annotation data
	if (verbose){ cat("reorganize annotation ...", "\n") }
	if (ncpus == 1){ annotation <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, CNann, chr, bpstart, bpend, simplify=TRUE)) }
	if (ncpus > 1){	annotation <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, CNann, chr, bpstart, bpend, simplify=TRUE)) }
	colnames(annotation) <- c("Chromosome", "Start", "End")
	if (ncpus == 1){ rownames(annotation) <- as.character(sapply(featuresAndWeights, .slhFuncCombineFeatureNames, rownames(CNann), simplify=TRUE)) }
	if (ncpus > 1){	rownames(annotation) <- as.character(sfSapply(featuresAndWeights, .slhFuncCombineFeatureNames, rownames(CNann), simplify=TRUE)) }
	if (length(rownames(annotation)) != length(unique(rownames(annotation)))){
		warning("modified feature names are returned")
		newRowNames <- rownames(annotation)
		slh <- rep(1, length(newRowNames))
		for (i in 1:length(unique(newRowNames))){
			ids <- which(newRowNames == unique(newRowNames)[i])
			slh[ids] <- 1:length(ids)
		}
		newRowNames <- paste(newRowNames, slh, sep="__")		
		rownames(annotation) <- newRowNames
	}
	fInfo <- new("AnnotatedDataFrame", data=data.frame(annotation), varMetadata=data.frame(labelDescription=c("Chromosome", "Start", "End")))
	if (ncpus > 1){ sfStop() }

	rownames(copynumberMat) <- rownames(annotation)
	rownames(segmentedMat) <- rownames(annotation)
	rownames(problossMat) <- rownames(annotation)
	rownames(probnormMat) <- rownames(annotation)
	rownames(probgainMat) <- rownames(annotation)
	rownames(callsMat) <- rownames(annotation)
	if (amplifications){ rownames(probampMat) <- rownames(annotation) }
	if (doubleLoss){ rownames(probdlossMat) <- rownames(annotation) }

	# make cghCall object
	if (verbose){ cat("merge into cghCall object ...", "\n") }
	if (!amplifications & !doubleLoss){ 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	}
	if (amplifications & !doubleLoss){ 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			probamp = data.matrix(data.frame(probampMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	}
	if (!amplifications & doubleLoss){ 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probdloss = data.matrix(data.frame(probdlossMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	}
	if (amplifications & doubleLoss){ 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probdloss = data.matrix(data.frame(probdlossMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			probamp = data.matrix(data.frame(probampMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	}
	return(CNdata)
}






cghCall2weightedSubset <- function(CNdata, featuresAndWeights, chr, bpstart, bpend, ncpus=1, verbose=TRUE){
   	############################################################################
	# function shrinks cghCall-object to a subset of the features
	############################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	# if (!(as.character(class(featureSubset))=="numeric" | as.character(class(featureSubset))=="integer" )){ stop("featureSubset of wrong class.") }
	# if (!(as.character(class(weights))=="numeric" | as.character(class(weights))=="integer" )){ stop("weights of wrong class.") }
	# if ( !(length(featureSubset) >= 1) ){ stop("featureSubset contains no row numbers.") }
	# if ( !(min(featureSubset) >= 1 | max(featureSubset) <= nrow(fData(CNdata))) ){ stop("featureSubset contains illegal row numbers.") }
	# if ( length(featureSubset) != length(weights) ){ stop("featureSubset and weights of unequal length.") }

	if (!(as.character(class(chr))=="numeric" | as.character(class(chr))=="integer" )){ stop("chr of wrong class.") }
	if (length(chr) != 1){ stop("chr of wrong length.") }
	if ( !(chr >= 1 | chr <= ncol(fData(CNdata))) ){ stop("chr exceeds number of columns in featureData.") }
	if (!(as.character(class(bpstart))=="numeric" | as.character(class(bpstart))=="integer" )){ stop("bpstart of wrong class.") }
	if (length(bpstart) != 1){ stop("bpstart of wrong length.") }
	if ( !(bpstart >= 1 | bpstart <= ncol(fData(CNdata))) ){ stop("bpstart exceeds number of columns in featureData.") }
	if (!(as.character(class(bpend))=="numeric" | as.character(class(bpend))=="integer" )){ stop("bpend of wrong class.") }
	if (length(bpend) != 1){ stop("bpend of wrong length.") }
	if ( !(bpend >= 1 | chr <= ncol(fData(CNdata))) ){ stop("bpend exceeds number of columns in featureData.") }

	# disecting the cghCall-object
	amplifications <- FALSE
	doubleLoss <- FALSE
	copynumberMat <- copynumber(CNdata)
	segmentedMat <- segmented(CNdata)
	problossMat <- probloss(CNdata)
	probnormMat <- probnorm(CNdata)
	probgainMat <- probgain(CNdata)
	if (!is.null(probamp(CNdata))){ probampMat <- probamp(CNdata); amplifications <- TRUE }
	if (!is.null(CGHbase::probdloss(CNdata))){ probdlossMat <- CGHbase::probdloss(CNdata); doubleLoss <- TRUE }
	CNann <- fData(CNdata)
	rm(CNdata)
	gc()

	# weighted subsetting of raw copy number data
	if (verbose){ cat("generate subsetted raw copy number matrix ...", "\n") }
	if (ncpus == 1){ copynumberMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, copynumberMat, simplify=TRUE))	}
	if (ncpus > 1){
		sfInit(parallel=TRUE, cpus=ncpus)
		sfLibrary(sigaR, verbose=FALSE)
		copynumberMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, copynumberMat, simplify=TRUE))	
	}

	# weighted subsetting of segmented data
	if (verbose){ cat("generate subsetted segmented matrix ...", "\n") }
	if (ncpus == 1){ segmentedMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, segmentedMat, simplify=TRUE)) }
	if (ncpus > 1){ segmentedMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, segmentedMat, simplify=TRUE)) }

	# weighted subsetting of call probability data
	if (verbose){ cat("generate subsetted call probability matrices ...", "\n") }
	if (ncpus == 1){ problossMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, problossMat, simplify=TRUE)) }
	if (ncpus > 1){ problossMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, problossMat, simplify=TRUE)) }
	if (ncpus == 1){ probnormMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, probnormMat, simplify=TRUE)) }
	if (ncpus > 1){ probnormMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, probnormMat, simplify=TRUE)) }
	if (ncpus == 1){ probgainMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, probgainMat, simplify=TRUE)) }
	if (ncpus > 1){ probgainMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, probgainMat, simplify=TRUE)) }
	if (amplifications){ 
		if (ncpus == 1){ probampMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, probampMat, simplify=TRUE)) }
		if (ncpus > 1){ probampMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, probampMat, simplify=TRUE)) }
	}
	if (doubleLoss){ 
		if (ncpus == 1){ probdlossMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, probdlossMat, simplify=TRUE)) }
		if (ncpus > 1){ probdlossMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, probdlossMat, simplify=TRUE)) }
	}

	# subsetting of call data
	if (verbose){ cat("generate subsetted call matrices ...", "\n")	}
	callsMat <- numeric()
	for (i in 1:ncol(problossMat)){
		if (!amplifications & !doubleLoss){
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 2) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 2) }
		} 
		if (amplifications & !doubleLoss){
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 2) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 2) }
		}
		if (!amplifications & doubleLoss){
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(probdlossMat[,i], problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 3) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(probdlossMat[,i], problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 3) }
		}
		if (amplifications & doubleLoss){
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(probdlossMat[,i], problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 3) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(probdlossMat[,i], problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 3) }
		}
	}	

	# weighted subsetting of annotation data
	if (verbose){ cat("reorganize annotation ...", "\n") }
	if (ncpus == 1){ annotation <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, CNann, chr, bpstart, bpend, simplify=TRUE)) }
	if (ncpus > 1){	annotation <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, CNann, chr, bpstart, bpend, simplify=TRUE)) }
	colnames(annotation) <- c("Chromosome", "Start", "End")
	if (ncpus == 1){ rownames(annotation) <- as.character(sapply(featuresAndWeights, .slhFuncCombineFeatureNames, rownames(CNann), simplify=TRUE)) }
	if (ncpus > 1){	rownames(annotation) <- as.character(sfSapply(featuresAndWeights, .slhFuncCombineFeatureNames, rownames(CNann), simplify=TRUE)) }
	if (length(rownames(annotation)) != length(unique(rownames(annotation)))){
		warning("modified feature names are returned")
		newRowNames <- rownames(annotation)
		slh <- rep(1, length(newRowNames))
		for (i in 1:length(unique(newRowNames))){
			ids <- which(newRowNames == unique(newRowNames)[i])
			slh[ids] <- 1:length(ids)
		}
		newRowNames <- paste(newRowNames, slh, sep="__")		
		rownames(annotation) <- newRowNames
	}
	fInfo <- new("AnnotatedDataFrame", data=data.frame(annotation), varMetadata=data.frame(labelDescription=c("Chromosome", "Start", "End")))
	if (ncpus > 1){ sfStop() }

	rownames(copynumberMat) <- rownames(annotation)
	rownames(segmentedMat) <- rownames(annotation)
	rownames(problossMat) <- rownames(annotation)
	rownames(probnormMat) <- rownames(annotation)
	rownames(probgainMat) <- rownames(annotation)
	rownames(callsMat) <- rownames(annotation)
	if (amplifications){ rownames(probampMat) <- rownames(annotation) }
	if (doubleLoss){ rownames(probdlossMat) <- rownames(annotation) }

	# make cghCall object
	if (verbose){ cat("merge into cghCall object ...", "\n") }
	if (!amplifications & !doubleLoss){ 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	}
	if (amplifications & !doubleLoss){ 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			probamp = data.matrix(data.frame(probampMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	}
	if (!amplifications & doubleLoss){ 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probdloss = data.matrix(data.frame(probdlossMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	}
	if (amplifications & doubleLoss){ 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probdloss = data.matrix(data.frame(probdlossMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			probamp = data.matrix(data.frame(probampMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	}
	return(CNdata)
}




cghCall2maximumSubset <- function(CNdata, featuresAndWeights, chr, bpstart, bpend, ncpus=1, verbose=TRUE){
   	############################################################################
	# function shrinks cghCall-object to a subset of the features
	############################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	# if (!(as.character(class(featureSubset))=="numeric" | as.character(class(featureSubset))=="integer" )){ stop("featureSubset of wrong class.") }
	# if (!(as.character(class(weights))=="numeric" | as.character(class(weights))=="integer" )){ stop("weights of wrong class.") }
	# if ( !(length(featureSubset) >= 1) ){ stop("featureSubset contains no row numbers.") }
	# if ( !(min(featureSubset) >= 1 | max(featureSubset) <= nrow(fData(CNdata))) ){ stop("featureSubset contains illegal row numbers.") }
	# if ( length(featureSubset) != length(weights) ){ stop("featureSubset and weights of unequal length.") }

	if (!(as.character(class(chr))=="numeric" | as.character(class(chr))=="integer" )){ stop("chr of wrong class.") }
	if (length(chr) != 1){ stop("chr of wrong length.") }
	if ( !(chr >= 1 | chr <= ncol(fData(CNdata))) ){ stop("chr exceeds number of columns in featureData.") }
	if (!(as.character(class(bpstart))=="numeric" | as.character(class(bpstart))=="integer" )){ stop("bpstart of wrong class.") }
	if (length(bpstart) != 1){ stop("bpstart of wrong length.") }
	if ( !(bpstart >= 1 | bpstart <= ncol(fData(CNdata))) ){ stop("bpstart exceeds number of columns in featureData.") }
	if (!(as.character(class(bpend))=="numeric" | as.character(class(bpend))=="integer" )){ stop("bpend of wrong class.") }
	if (length(bpend) != 1){ stop("bpend of wrong length.") }
	if ( !(bpend >= 1 | chr <= ncol(fData(CNdata))) ){ stop("bpend exceeds number of columns in featureData.") }

	# disecting the cghCall-object
	amplifications <- FALSE
	copynumberMat <- copynumber(CNdata)
	segmentedMat <- segmented(CNdata)
	problossMat <- probloss(CNdata)
	probnormMat <- probnorm(CNdata)
	probgainMat <- probgain(CNdata)
	if (!is.null(probamp(CNdata))){ probampMat <- probamp(CNdata); amplifications <- TRUE }
	CNann <- fData(CNdata)
	rm(CNdata)
	gc()

	# identify maximum deviating features
	if (verbose){ cat("generate subsetted segmented matrix ...", "\n") }
	if (ncpus == 1){ idMat <- t(sapply(featuresAndWeights, .slhFuncMaxCGHcallData, segmentedMat, simplify=TRUE)) }
	if (ncpus > 1){ 
		sfInit(parallel=TRUE, cpus=ncpus)
		sfLibrary(sigaR, verbose=FALSE)
		idMat <- t(sfSapply(featuresAndWeights, .slhFuncMaxCGHcallData, segmentedMat, simplify=TRUE)) 
	}

	# subsetting of segmented data
	if (ncpus == 1){ copynumberMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, copynumberMat, simplify=TRUE)) }
	if (ncpus > 1){ copynumberMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, copynumberMat, simplify=TRUE)) }
	if (ncpus == 1){ segmentedMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, segmentedMat, simplify=TRUE)) }
	if (ncpus > 1){ segmentedMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, segmentedMat, simplify=TRUE)) }
	if (ncpus == 1){ problossMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, problossMat, simplify=TRUE)) }
	if (ncpus > 1){ problossMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, problossMat, simplify=TRUE)) }
	if (ncpus == 1){ probnormMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probnormMat, simplify=TRUE)) }
	if (ncpus > 1){ probnormMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probnormMat, simplify=TRUE)) }
	if (ncpus == 1){ probgainMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probgainMat, simplify=TRUE)) }
	if (ncpus > 1){ probgainMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probgainMat, simplify=TRUE)) }
	if (amplifications){ 
		if (ncpus == 1){ probampMat <- t(sapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probampMat, simplify=TRUE)) }
		if (ncpus > 1){ probampMat <- t(sfSapply(1:nrow(idMat), .slhFunctMaxFeatures, featuresAndWeights, idMat, probampMat, simplify=TRUE)) }
	}

	# subsetting of call data
	if (verbose){ cat("generate subsetted call matrices ...", "\n")	}
	callsMat <- numeric()
	for (i in 1:ncol(problossMat)){
		if (!amplifications){
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 2) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i]), 1, which.max) - 2) }
		} else {
			if (ncpus == 1){ callsMat <- cbind(callsMat, apply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 2) }
			if (ncpus > 1){ callsMat <- cbind(callsMat, sfApply(cbind(problossMat[,i], probnormMat[,i], probgainMat[,i], probampMat[,i]), 1, which.max) - 2) }
		}
	}	

	# weighted subsetting of annotation data
	if (verbose){ cat("reorganize annotation ...", "\n") }
	if (ncpus == 1){ annotation <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, CNann, chr, bpstart, bpend, simplify=TRUE)) }
	if (ncpus > 1){	annotation <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, CNann, chr, bpstart, bpend, simplify=TRUE)) }
	colnames(annotation) <- c("Chromosome", "Start", "End")
	if (ncpus == 1){ rownames(annotation) <- as.character(sapply(featuresAndWeights, .slhFuncCombineFeatureNames, rownames(CNann), simplify=TRUE)) }
	if (ncpus > 1){	rownames(annotation) <- as.character(sfSapply(featuresAndWeights, .slhFuncCombineFeatureNames, rownames(CNann), simplify=TRUE)) }
	if (length(rownames(annotation)) != length(unique(rownames(annotation)))){
		warning("modified feature names are returned")
		newRowNames <- rownames(annotation)
		slh <- rep(1, length(newRowNames))
		for (i in 1:length(unique(newRowNames))){
			ids <- which(newRowNames == unique(newRowNames)[i])
			slh[ids] <- 1:length(ids)
		}
		newRowNames <- paste(newRowNames, slh, sep="__")		
		rownames(annotation) <- newRowNames
	}
	fInfo <- new("AnnotatedDataFrame", data=data.frame(annotation), varMetadata=data.frame(labelDescription=c("Chromosome", "Start", "End")))
	if (ncpus > 1){ sfStop() }

	rownames(copynumberMat) <- rownames(annotation)
	rownames(segmentedMat) <- rownames(annotation)
	rownames(problossMat) <- rownames(annotation)
	rownames(probnormMat) <- rownames(annotation)
	rownames(probgainMat) <- rownames(annotation)
	rownames(callsMat) <- rownames(annotation)
	if (amplifications){ rownames(probampMat) <- rownames(annotation) }

	# make cghCall object
	if (verbose){ cat("merge into cghCall object ...", "\n") }
	if (!amplifications){ 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	} else { 
		CNdata <- new('cghCall', 
			copynumber = data.matrix(data.frame(copynumberMat, row.names=rownames(annotation))), 
			segmented = data.matrix(data.frame(segmentedMat, row.names=rownames(annotation))), 
			calls = data.matrix(data.frame(callsMat, row.names=rownames(annotation))), 
			probloss = data.matrix(data.frame(problossMat, row.names=rownames(annotation))), 
			probnorm = data.matrix(data.frame(probnormMat, row.names=rownames(annotation))), 
			probgain = data.matrix(data.frame(probgainMat, row.names=rownames(annotation))), 
			probamp = data.matrix(data.frame(probampMat, row.names=rownames(annotation))), 
			featureData = fInfo) 
	}
	return(CNdata)
}


cghCall2subset <- function(CNdata, featureSubset, verbose=TRUE){
   	############################################################################
	# function shrinks cghCall-object to a subset of the features
	############################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	if (!(as.character(class(featureSubset))=="numeric" | as.character(class(featureSubset))=="integer" )){ stop("featureSubset of wrong class.") }
	if ( !(length(featureSubset) >= 1) ){ stop("featureSubset contains no row numbers.") }
	if ( !(min(featureSubset) >= 1 | max(featureSubset) <= nrow(fData(CNdata))) ){ stop("featureSubset contains illegal row numbers.") }

	# sort out new annotation info
	fd <- fData(CNdata)[featureSubset, , drop=FALSE]
	newRowNames <- rownames(copynumber(CNdata))[featureSubset]
	if (length(newRowNames) != length(unique(newRowNames))){
		warning("modified feature names are returned")
		slh <- rep(1, length(newRowNames))
		for (i in 1:length(unique(newRowNames))){
			ids <- which(newRowNames == unique(newRowNames)[i])
			slh[ids] <- 1:length(ids)
		}
		newRowNames <- paste(newRowNames, slh, sep="__")		
	}
	rownames(fd) <- newRowNames
	# sNames <- sampleNames(CNdata)
	metaData <- data.frame(labelDescription=colnames(fd)) 
	fd <- new("AnnotatedDataFrame", data=data.frame(fd), varMetadata=metaData)

	if (is.null(CGHbase::probdloss(CNdata)[featureSubset, , drop=FALSE])){ 
		if (is.null(probamp(CNdata)[featureSubset, , drop=FALSE])){ 
			CNdata <- new('cghCall', 
				copynumber = data.matrix(data.frame(copynumber(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				segmented = data.matrix(data.frame(segmented(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				calls = data.matrix(data.frame(calls(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probloss = data.matrix(data.frame(probloss(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probnorm = data.matrix(data.frame(probnorm(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probgain = data.matrix(data.frame(probgain(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				featureData = fd) 
		} else { 
			CNdata <- new('cghCall', 
				copynumber = data.matrix(data.frame(copynumber(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				segmented = data.matrix(data.frame(segmented(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				calls = data.matrix(data.frame(calls(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probloss = data.matrix(data.frame(probloss(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probnorm = data.matrix(data.frame(probnorm(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probgain = data.matrix(data.frame(probgain(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probamp = data.matrix(data.frame(probamp(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				featureData = fd) 
		}
	} else {
		if (is.null(probamp(CNdata)[featureSubset, , drop=FALSE])){ 
			CNdata <- new('cghCall', 
				copynumber = data.matrix(data.frame(copynumber(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				segmented = data.matrix(data.frame(segmented(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				calls = data.matrix(data.frame(calls(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probdloss = data.matrix(data.frame(CGHbase::probdloss(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)),  
				probloss = data.matrix(data.frame(probloss(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)),
				probnorm = data.matrix(data.frame(probnorm(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probgain = data.matrix(data.frame(probgain(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				featureData = fd) 
		} else { 
			CNdata <- new('cghCall', 
				copynumber = data.matrix(data.frame(copynumber(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				segmented = data.matrix(data.frame(segmented(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				calls = data.matrix(data.frame(calls(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probdloss = data.matrix(data.frame(CGHbase::probdloss(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)),  
				probloss = data.matrix(data.frame(probloss(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probnorm = data.matrix(data.frame(probnorm(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probgain = data.matrix(data.frame(probgain(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				probamp = data.matrix(data.frame(probamp(CNdata)[featureSubset, , drop=FALSE], row.names=newRowNames)), 
				featureData = fd) 
		}		
	}
	return(CNdata)
}





merge2cghCalls <- function(CNdata1, CNdata2, verbose=TRUE){
   	############################################################################
	# function merges cghCall-objects into one cghCall-object
	############################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata1)) != "cghCall"){ stop("CNdata1 not of class cghCall.") }
	if (as.character(class(CNdata2)) != "cghCall"){ stop("CNdata2 not of class cghCall.") }
	if ( dim(CNdata1)[2] != dim(CNdata2)[2] ){ stop("CNdata1 and CNdata2 have an unequal number of samples.") }
	if ( dim(fData(CNdata1))[2] != dim(fData(CNdata2))[2] ){ stop("CNdata1 and CNdata2 have an unequal annotation columns.") }
	amp1 <- is.null(probamp(CNdata1))
	amp2 <- is.null(probamp(CNdata2))
	if (!(amp1 == amp2)){ stop("One of the cghCall-objects has amplification probabilities whereas the other has not.") }

	# prepare new annotation
	annotation <- rbind(fData(CNdata1), fData(CNdata2))
	newRowNames <- c(rownames(fData(CNdata1)), rownames(fData(CNdata2)))
	if (length(newRowNames) != length(unique(newRowNames))){
		warning("modified feature names are returned")
		slh <- rep(1, length(newRowNames))
		for (i in 1:length(unique(newRowNames))){
			ids <- which(newRowNames == unique(newRowNames)[i])
			slh[ids] <- 1:length(ids)
		}
		newRowNames <- paste(newRowNames, slh, sep="__")		
	}
	rownames(annotation) <- newRowNames
	metaData <- data.frame(labelDescription=colnames(annotation)) 
	annotation <- new("AnnotatedDataFrame", data=data.frame(annotation), varMetadata=metaData)

	if (is.null(CGHbase::probdloss(CNdata1))){ 
		if (is.null(probamp(CNdata1))){ 
			CNdata <- new('cghCall', 
				copynumber = data.matrix(data.frame(rbind(copynumber(CNdata1), copynumber(CNdata2)), row.names=newRowNames)), 
				segmented = data.matrix(data.frame(rbind(segmented(CNdata1), segmented(CNdata2)), row.names=newRowNames)), 
				calls = data.matrix(data.frame(rbind(calls(CNdata1), calls(CNdata2)), row.names=newRowNames)), 
				probloss = data.matrix(data.frame(rbind(probloss(CNdata1), probloss(CNdata2)), row.names=newRowNames)), 
				probnorm = data.matrix(data.frame(rbind(probnorm(CNdata1), probnorm(CNdata2)), row.names=newRowNames)), 
				probgain = data.matrix(data.frame(rbind(probgain(CNdata1), probgain(CNdata2)), row.names=newRowNames)), 
				featureData = annotation) 
		} else { 
			CNdata <- new('cghCall', 
				copynumber = data.matrix(data.frame(rbind(copynumber(CNdata1), copynumber(CNdata2)), row.names=newRowNames)), 
				segmented = data.matrix(data.frame(rbind(segmented(CNdata1), segmented(CNdata2)), row.names=newRowNames)), 
				calls = data.matrix(data.frame(rbind(calls(CNdata1), calls(CNdata2)), row.names=newRowNames)), 
				probloss = data.matrix(data.frame(rbind(probloss(CNdata1), probloss(CNdata2)), row.names=newRowNames)), 
				probnorm = data.matrix(data.frame(rbind(probnorm(CNdata1), probnorm(CNdata2)), row.names=newRowNames)), 
				probgain = data.matrix(data.frame(rbind(probgain(CNdata1), probgain(CNdata2)), row.names=newRowNames)), 
				probamp = data.matrix(data.frame(rbind(probamp(CNdata1), probamp(CNdata2)), row.names=newRowNames)), 
				featureData = annotation) 
		}
	} else {
		if (is.null(probamp(CNdata1))){ 
			CNdata <- new('cghCall', 
				copynumber = data.matrix(data.frame(rbind(copynumber(CNdata1), copynumber(CNdata2)), row.names=newRowNames)), 
				segmented = data.matrix(data.frame(rbind(segmented(CNdata1), segmented(CNdata2)), row.names=newRowNames)), 
				calls = data.matrix(data.frame(rbind(calls(CNdata1), calls(CNdata2)), row.names=newRowNames)), 
				probdloss = data.matrix(data.frame(rbind(CGHbase::probdloss(CNdata1), probloss(CNdata2)), row.names=newRowNames)), 
				probloss = data.matrix(data.frame(rbind(probloss(CNdata1), probloss(CNdata2)), row.names=newRowNames)), 
				probnorm = data.matrix(data.frame(rbind(probnorm(CNdata1), probnorm(CNdata2)), row.names=newRowNames)), 
				probgain = data.matrix(data.frame(rbind(probgain(CNdata1), probgain(CNdata2)), row.names=newRowNames)), 
				featureData = annotation) 
		} else { 
			CNdata <- new('cghCall', 
				copynumber = data.matrix(data.frame(rbind(copynumber(CNdata1), copynumber(CNdata2)), row.names=newRowNames)), 
				segmented = data.matrix(data.frame(rbind(segmented(CNdata1), segmented(CNdata2)), row.names=newRowNames)), 
				calls = data.matrix(data.frame(rbind(calls(CNdata1), calls(CNdata2)), row.names=newRowNames)), 
				probdloss = data.matrix(data.frame(rbind(CGHbase::probdloss(CNdata1), probloss(CNdata2)), row.names=newRowNames)), 
				probloss = data.matrix(data.frame(rbind(probloss(CNdata1), probloss(CNdata2)), row.names=newRowNames)), 
				probnorm = data.matrix(data.frame(rbind(probnorm(CNdata1), probnorm(CNdata2)), row.names=newRowNames)), 
				probgain = data.matrix(data.frame(rbind(probgain(CNdata1), probgain(CNdata2)), row.names=newRowNames)), 
				probamp = data.matrix(data.frame(rbind(probamp(CNdata1), probamp(CNdata2)), row.names=newRowNames)), 
				featureData = annotation) 
		}
	}

	return(CNdata)
}




cghCall2order <- function(CNdata, chr, bpstart, verbose=TRUE){
   	###########################################################################################
	# function that orders cghCall-object genomically
	###########################################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	if (!(as.character(class(chr))=="numeric" | as.character(class(chr))=="integer" )){ stop("chr of wrong class.") }
	if (length(chr) != 1){ stop("chr of wrong length.") }
	if ( !(chr >= 1 | chr <= ncol(fData(CNdata))) ){ stop("chr exceeds number of columns in featureData.") }
	if (!(as.character(class(bpstart))=="numeric" | as.character(class(bpstart))=="integer" )){ stop("bpstart of wrong class.") }
	if (length(bpstart) != 1){ stop("bpstart of wrong length.") }
	if ( !(bpstart >= 1 | bpstart <= ncol(fData(CNdata))) ){ stop("GEbpstart exceeds number of columns in featureData.") }

	# extract annotation data for ordering
	CNann <- fData(CNdata)[, c(chr, bpstart), drop=FALSE]

	# in case columns are of wrong class ....
	for (i in 1:2){
		if (is.factor(CNann[,i])){ CNann[,i] <- as.numeric(levels(CNann[,i]))[CNann[,i]] }
		if (is.character(CNann[,i])){ CNann[,i] <- as.numeric(CNann[,i]) }
	}

	# order ExpressionSet-object genomically	
	fData(CNdata) <- fData(CNdata)[order(CNann[,1], CNann[,2]),]
	copynumber(CNdata) <- copynumber(CNdata)[order(CNann[,1], CNann[,2]),]
	segmented(CNdata) <- segmented(CNdata)[order(CNann[,1], CNann[,2]),]
	calls(CNdata) <- calls(CNdata)[order(CNann[,1], CNann[,2]),]
	probloss(CNdata) <- probloss(CNdata)[order(CNann[,1], CNann[,2]),]
	probnorm(CNdata) <- probnorm(CNdata)[order(CNann[,1], CNann[,2]),]
	probgain(CNdata) <- probgain(CNdata)[order(CNann[,1], CNann[,2]),]
	if (is.null(probamp(CNdata))){ } else { probamp(CNdata) <- probamp(CNdata)[order(CNann[,1], CNann[,2]),] }
	if (is.null(CGHbase::probdloss(CNdata))){ } else { CGHbase::probdloss(CNdata) <- CGHbase::probdloss(CNdata)[order(CNann[,1], CNann[,2]),] }
	return(CNdata)
}



cghSeg2order <- function(CNdata, chr, bpstart, verbose=TRUE){
   	###########################################################################################
	# function that orders cghCall-object genomically
	###########################################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata)) != "cghSeg"){ stop("CNdata not of class cghSeg.") }
	if (!(as.character(class(chr))=="numeric" | as.character(class(chr))=="integer" )){ stop("chr of wrong class.") }
	if (length(chr) != 1){ stop("chr of wrong length.") }
	if ( !(chr >= 1 | chr <= ncol(fData(CNdata))) ){ stop("chr exceeds number of columns in featureData.") }
	if (!(as.character(class(bpstart))=="numeric" | as.character(class(bpstart))=="integer" )){ stop("bpstart of wrong class.") }
	if (length(bpstart) != 1){ stop("bpstart of wrong length.") }
	if ( !(bpstart >= 1 | bpstart <= ncol(fData(CNdata))) ){ stop("bpstart exceeds number of columns in featureData.") }

	# extract annotation data for ordering
	CNann <- fData(CNdata)[, c(chr, bpstart), drop=FALSE]

	# in case columns are of wrong class ....
	for (i in 1:2){
		if (is.factor(CNann[,i])){ CNann[,i] <- as.numeric(levels(CNann[,i]))[CNann[,i]] }
		if (is.character(CNann[,i])){ CNann[,i] <- as.numeric(CNann[,i]) }
	}

	# order ExpressionSet-object genomically	
	fData(CNdata) <- fData(CNdata)[order(CNann[,1], CNann[,2]),]
	copynumber(CNdata) <- copynumber(CNdata)[order(CNann[,1], CNann[,2]),]
	segmented(CNdata) <- segmented(CNdata)[order(CNann[,1], CNann[,2]),]
	return(CNdata)
}



cghCall2cghSeg <- function (CNdata, verbose = TRUE){
    if (verbose) {
        cat("perform input checks...", "\n")
    }
    if (as.character(class(CNdata)) != "cghCall") {
        stop("CNdata not of class cghCall.")
    }
   

    fd <- fData(CNdata)
    newRowNames <- rownames(copynumber(CNdata))
    if (length(newRowNames) != length(unique(newRowNames))) {
        warning("modified feature names are returned")
        slh <- rep(1, length(newRowNames))
        for (i in 1:length(unique(newRowNames))) {
            ids <- which(newRowNames == unique(newRowNames)[i])
            slh[ids] <- 1:length(ids)
        }
        newRowNames <- paste(newRowNames, slh, sep = "__")
    }
    rownames(fd) <- newRowNames
    metaData <- data.frame(labelDescription = colnames(fd))
    fd <- new("AnnotatedDataFrame", data = data.frame(fd), varMetadata = metaData)
    CNdata <- new("cghSeg", copynumber = data.matrix(data.frame(copynumber(CNdata), row.names = newRowNames)), 
		segmented = data.matrix(data.frame(segmented(CNdata), row.names = newRowNames)), featureData = fd)

    return(CNdata)
}



cghSeg2subset <- function (CNdata, featureSubset, verbose = TRUE){
    if (verbose) {
        cat("perform input checks...", "\n")
    }
    if (as.character(class(CNdata)) != "cghSeg") {
        stop("CNdata not of class cghSeg.")
    }
    if (!(as.character(class(featureSubset)) == "numeric" | as.character(class(featureSubset)) == 
        "integer")) {
        stop("featureSubset of wrong class.")
    }
    if (!(length(featureSubset) >= 1)) {
        stop("featureSubset contains no row numbers.")
    }
    if (!(min(featureSubset) >= 1 | max(featureSubset) <= nrow(fData(CNdata)))) {
        stop("featureSubset contains illegal row numbers.")
    }
    fd <- fData(CNdata)[featureSubset, , drop = FALSE]
    newRowNames <- rownames(copynumber(CNdata))[featureSubset]
    if (length(newRowNames) != length(unique(newRowNames))) {
        warning("modified feature names are returned")
        slh <- rep(1, length(newRowNames))
        for (i in 1:length(unique(newRowNames))) {
            ids <- which(newRowNames == unique(newRowNames)[i])
            slh[ids] <- 1:length(ids)
        }
        newRowNames <- paste(newRowNames, slh, sep = "__")
    }
    rownames(fd) <- newRowNames
    metaData <- data.frame(labelDescription = colnames(fd))
    fd <- new("AnnotatedDataFrame", data = data.frame(fd), varMetadata = metaData)
    CNdata <- new("cghSeg", copynumber = data.matrix(data.frame(copynumber(CNdata)[featureSubset, 
            , drop = FALSE], row.names = newRowNames)), segmented = data.matrix(data.frame(segmented(CNdata)[featureSubset, 
            , drop = FALSE], row.names = newRowNames)), featureData = fd)

    return(CNdata)
}




cghSeg2weightedSubset <- function (CNdata, featuresAndWeights, chr, bpstart, bpend, ncpus = 1, verbose = TRUE){
    if (verbose) {
        cat("perform input checks...", "\n")
    }
    if (as.character(class(CNdata)) != "cghSeg") {
        stop("CNdata not of class cghSeg.")
    }
    if (as.character(class(CNdata)) != "cghSeg") {
        stop("CNdata not of class cghSeg.")
    }
    if (!(as.character(class(chr)) == "numeric" | as.character(class(chr)) == 
        "integer")) {
        stop("chr of wrong class.")
    }
    if (length(chr) != 1) {
        stop("chr of wrong length.")
    }
    if (!(chr >= 1 | chr <= ncol(fData(CNdata)))) {
        stop("chr exceeds number of columns in featureData.")
    }
    if (!(as.character(class(bpstart)) == "numeric" | as.character(class(bpstart)) == 
        "integer")) {
        stop("bpstart of wrong class.")
    }
    if (length(bpstart) != 1) {
        stop("bpstart of wrong length.")
    }
    if (!(bpstart >= 1 | bpstart <= ncol(fData(CNdata)))) {
        stop("bpstart exceeds number of columns in featureData.")
    }
    if (!(as.character(class(bpend)) == "numeric" | as.character(class(bpend)) == 
        "integer")) {
        stop("bpend of wrong class.")
    }
    if (length(bpend) != 1) {
        stop("bpend of wrong length.")
    }
    if (!(bpend >= 1 | chr <= ncol(fData(CNdata)))) {
        stop("bpend exceeds number of columns in featureData.")
    }
    copynumberMat <- copynumber(CNdata)
    segmentedMat <- segmented(CNdata)
    CNann <- fData(CNdata)
    rm(CNdata)
    gc()
    if (verbose) {
        cat("generate subsetted raw copy number matrix ...", 
            "\n")
    }
    if (ncpus == 1) {
        copynumberMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, 
            copynumberMat, simplify = TRUE))
    }
    if (ncpus > 1) {
        sfInit(parallel = TRUE, cpus = ncpus)
        sfLibrary(sigaR, verbose = FALSE)
        copynumberMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, 
            copynumberMat, simplify = TRUE))
    }
    if (verbose) {
        cat("generate subsetted segmented matrix ...", "\n")
    }
    if (ncpus == 1) {
        segmentedMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, 
            segmentedMat, simplify = TRUE))
    }
    if (ncpus > 1) {
        segmentedMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, 
            segmentedMat, simplify = TRUE))
    }
    if (verbose) {
        cat("generate subsetted call matrices ...", "\n")
    }
    if (verbose) {
        cat("reorganize annotation ...", "\n")
    }
    if (ncpus == 1) {
        annotation <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, 
            CNann, chr, bpstart, bpend, simplify = TRUE))
    }
    if (ncpus > 1) {
        annotation <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, 
            CNann, chr, bpstart, bpend, simplify = TRUE))
    }
    colnames(annotation) <- c("Chromosome", "Start", "End")
    if (ncpus == 1) {
        rownames(annotation) <- as.character(sapply(featuresAndWeights, 
            .slhFuncCombineFeatureNames, rownames(CNann), simplify = TRUE))
    }
    if (ncpus > 1) {
        rownames(annotation) <- as.character(sfSapply(featuresAndWeights, 
            .slhFuncCombineFeatureNames, rownames(CNann), simplify = TRUE))
    }
    if (length(rownames(annotation)) != length(unique(rownames(annotation)))) {
        warning("modified feature names are returned")
        newRowNames <- rownames(annotation)
        slh <- rep(1, length(newRowNames))
        for (i in 1:length(unique(newRowNames))) {
            ids <- which(newRowNames == unique(newRowNames)[i])
            slh[ids] <- 1:length(ids)
        }
        newRowNames <- paste(newRowNames, slh, sep = "__")
        rownames(annotation) <- newRowNames
    }
    fInfo <- new("AnnotatedDataFrame", data = data.frame(annotation), 
        varMetadata = data.frame(labelDescription = c("Chromosome", 
            "Start", "End")))
    if (ncpus > 1) {
        sfStop()
    }
    rownames(copynumberMat) <- rownames(annotation)
    rownames(segmentedMat) <- rownames(annotation)
    if (verbose) {
        cat("merge into cghSeg object ...", "\n")
    }
    CNdata <- new("cghSeg", copynumber = data.matrix(data.frame(copynumberMat, 
		row.names = rownames(annotation))), segmented = data.matrix(data.frame(segmentedMat, 
		row.names = rownames(annotation))), featureData = fInfo)
    return(CNdata)
}





splitMatchingAtBreakpoints <- function(matchedIDs, CNdata){
   	############################################################################
	# function splits matching list-object at breakpoints
	############################################################################

	nSamplesWithBreakpoint <- function(X, CNdata){ 
		sum(apply(segmented(CNdata)[X[,2], , drop=FALSE], 2, function(Z){ length(unique(Z)) }) > 1) 
	}

	matrixRows2listwithReg <- function(X, regInfo){
		subList <- sapply(1:nrow(regInfo), function(u, Z, regInfo){ Z[regInfo[u,1]:regInfo[u,2], , drop=FALSE] }, Z=X, regInfo, simplify=FALSE)
		names(subList) <- rep(X[1, 1], nrow(X))
		return(subList)
	}
	
	makeRegions <- function(CNsegData){
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

	expandMatches <- function(id, matchedIDs, CNdata){
		nEntries <- length(matchedIDs)
		mFeatures <- nrow(matchedIDs[[id]])
		regionsWithinGene <- makeRegions(segmented(CNdata)[as.numeric(matchedIDs[[id]][,2]),,drop=FALSE])
		mIDsInsert <- matrixRows2listwithReg(matchedIDs[[id]], regionsWithinGene)
		if (id < nEntries){ mIDsEnd <- matchedIDs[(id+1):length(matchedIDs)] } else { mIDsEnd <- NULL }
		if (id > 1){ mIDsBegin <- matchedIDs[1:(id-1)] } else { mIDsBegin <- NULL }
		if (id > 1 & id < nEntries){ 
			matchedIDs <- mIDsBegin
			matchedIDs[(length(matchedIDs)+1):(length(matchedIDs)+mFeatures)] <- mIDsInsert
			names(matchedIDs)[(length(matchedIDs)-mFeatures+1):length(matchedIDs)] <- names(mIDsInsert)
			matchedIDs[(length(matchedIDs)+1):(length(matchedIDs)+length(mIDsEnd))] <- mIDsEnd
			names(matchedIDs)[(length(matchedIDs)-length(mIDsEnd)+1):length(matchedIDs)] <- names(mIDsEnd)
		}
		if (id == 1 & id < nEntries){ 
			matchedIDs <- mIDsInsert
			matchedIDs[(length(matchedIDs)+1):(length(matchedIDs)+length(mIDsEnd))] <- mIDsEnd
			names(matchedIDs)[(length(matchedIDs)-length(mIDsEnd)+1):length(matchedIDs)] <- names(mIDsEnd)
		}
		if (id > 1 & id == nEntries){ 
			matchedIDs <- mIDsBegin
			matchedIDs[(length(matchedIDs)+1):(length(matchedIDs)+mFeatures)] <- mIDsInsert
			names(matchedIDs)[(length(matchedIDs)-mFeatures+1):length(matchedIDs)] <- names(mIDsInsert)
		}
		return(matchedIDs)
	}

	# determine matches with breakpoints
	IDs <- sort(which(sapply(matchedIDs, nSamplesWithBreakpoint, CNdata, simplify=TRUE) > 1), decreasing=TRUE)

	# split at breakpoints
	if (length(IDs) > 0){
		for (j in 1:length(IDs)){
			matchedIDs <- expandMatches(IDs[j], matchedIDs, CNdata)
		}
	}
	return(matchedIDs)
}






merge2ExpressionSets <- function(GEdata1, GEdata2, verbose=TRUE){
   	############################################################################
	# function merges ExpressionSet-objects into one ExpressionSet-object
	############################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if ( !(as.character(class(GEdata1)) == "ExpressionSet" | as.character(class(GEdata1)) == "eSet") ){ stop("GEdata1 not of class ExpressionSet.") }
	if ( !(as.character(class(GEdata2)) == "ExpressionSet" | as.character(class(GEdata2)) == "eSet") ){ stop("GEdata2 not of class ExpressionSet.") }
	if ( dim(GEdata1)[2] != dim(GEdata2)[2] ){ stop("GEdata1 and GEdata2 have an unequal number of samples.") }
	if ( dim(fData(GEdata1))[2] != dim(fData(GEdata2))[2] ){ stop("GEdata1 and GEdata2 have an unequal annotation columns.") }

	# prepare new annotation
	annotation <- rbind(fData(GEdata1), fData(GEdata2))
	newRowNames <- c(rownames(fData(GEdata1)), rownames(fData(GEdata2)))
	if (length(newRowNames) != length(unique(newRowNames))){
		warning("modified feature names are returned")
		slh <- rep(1, length(newRowNames))
		for (i in 1:length(unique(newRowNames))){
			ids <- which(newRowNames == unique(newRowNames)[i])
			slh[ids] <- 1:length(ids)
		}
		newRowNames <- paste(newRowNames, slh, sep="__")		
	}
	rownames(annotation) <- newRowNames
	metaData <- data.frame(labelDescription=colnames(annotation)) 
	probeinfo <- new("AnnotatedDataFrame", data=data.frame(annotation), varMetadata=metaData)
	expMat <- rbind(exprs(GEdata1), exprs(GEdata2))
	rownames(expMat) <- newRowNames
	return(new("ExpressionSet", exprs=expMat, featureData = probeinfo))
}





nBreakpoints <- function(featuresAndWeights, CNdata){
   	############################################################################
	# function calculates the number of samples with at least one breakpoint 
	# within the transcript to which the DNA copy number data have been matched
	############################################################################

	nSamplesWithBreakpoint <- function(X, CNdata){ 
		sum(apply(segmented(CNdata)[X[,1], , drop=FALSE], 2, function(Z){ length(unique(Z)) }) > 1) 
	}
	return(sapply(featuresAndWeights, nSamplesWithBreakpoint, CNdata, simplify=TRUE))
}


expandMatching2singleIDs <- function(matchedIDs){
   	############################################################################
	# function unfolds matching list-object 
	############################################################################

	expandMatches <- function(id, matchedIDs){
		nEntries <- length(matchedIDs)
		mFeatures <- nrow(matchedIDs[[id]])
		mIDsInsert <- .matrixRows2list(matchedIDs[[id]])
		if (id < nEntries){ mIDsEnd <- matchedIDs[(id+1):length(matchedIDs)] } else { mIDsEnd <- NULL }
		if (id > 1){ mIDsBegin <- matchedIDs[1:(id-1)] } else { mIDsBegin <- NULL }
		if (id > 1 & id < nEntries){ 
			matchedIDs <- mIDsBegin
			matchedIDs[(length(matchedIDs)+1):(length(matchedIDs)+mFeatures)] <- mIDsInsert
			names(matchedIDs)[(length(matchedIDs)-mFeatures+1):length(matchedIDs)] <- names(mIDsInsert)
			matchedIDs[(length(matchedIDs)+1):(length(matchedIDs)+length(mIDsEnd))] <- mIDsEnd
			names(matchedIDs)[(length(matchedIDs)-length(mIDsEnd)+1):length(matchedIDs)] <- names(mIDsEnd)
		}
		if (id == 1 & id < nEntries){ 
			matchedIDs <- mIDsInsert
			matchedIDs[(length(matchedIDs)+1):(length(matchedIDs)+length(mIDsEnd))] <- mIDsEnd
			names(matchedIDs)[(length(matchedIDs)-length(mIDsEnd)+1):length(matchedIDs)] <- names(mIDsEnd)
		}
		if (id > 1 & id == nEntries){ 
			matchedIDs <- mIDsBegin
			matchedIDs[(length(matchedIDs)+1):(length(matchedIDs)+mFeatures)] <- mIDsInsert
			names(matchedIDs)[(length(matchedIDs)-mFeatures+1):length(matchedIDs)] <- names(mIDsInsert)
		}
		return(matchedIDs)
	}

	# actual expansion
	IDs <- sort(as.numeric(which(sapply(matchedIDs, function(X){ nrow(X) }, simplify=TRUE) > 1)), decreasing=TRUE)
	if (length(IDs) > 1){
		for (j in 1:length(IDs)){
			matchedIDs <- expandMatches(IDs[j], matchedIDs)
		}
	}
	return(matchedIDs)
}





ExpressionSet2weightedSubset <- function(GEdata, featuresAndWeights, chr, bpstart, bpend, ncpus=1, verbose=TRUE){
   	############################################################################
	# function shrinks ExpressionSet-object to a weighted subset of the features
	############################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if ( !(as.character(class(GEdata)) == "ExpressionSet" | as.character(class(GEdata)) == "eSet") ){ stop("GEdata not of class ExpressionSet.") }
	if (!(as.character(class(featuresAndWeights))=="list")){ stop("featuresAndWeights of wrong class.") }

	if (!(as.character(class(chr))=="numeric" | as.character(class(chr))=="integer" )){ stop("chr of wrong class.") }
	if (length(chr) != 1){ stop("chr of wrong length.") }
	if ( !(chr >= 1 | chr <= ncol(fData(GEdata))) ){ stop("chr exceeds number of columns in featureData.") }
	if (!(as.character(class(bpstart))=="numeric" | as.character(class(bpstart))=="integer" )){ stop("bpstart of wrong class.") }
	if (length(bpstart) != 1){ stop("bpstart of wrong length.") }
	if ( !(bpstart >= 1 | bpstart <= ncol(fData(GEdata))) ){ stop("bpstart exceeds number of columns in featureData.") }
	if (!(as.character(class(bpend))=="numeric" | as.character(class(bpend))=="integer" )){ stop("bpend of wrong class.") }
	if (length(bpend) != 1){ stop("bpend of wrong length.") }
	if ( !(bpend >= 1 | chr <= ncol(fData(GEdata))) ){ stop("bpend exceeds number of columns in featureData.") }

	# disecting the ExpressionSet-object
	exprMat <- exprs(GEdata)
	GEann <- fData(GEdata)
	rm(GEdata)
	gc()

	# weighted subsetting of expression data
	if (verbose){ cat("generate matched expression matrix ...", "\n") }
	if (ncpus == 1){ exprMat <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, exprMat, simplify=TRUE)) }
	if (ncpus > 1){ 
		sfInit(parallel=TRUE, cpus=ncpus)
		sfLibrary(sigaR, verbose=FALSE)
		exprMat <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetData, exprMat, simplify=TRUE))
	}

	# weighted subsetting of annotation data
	if (verbose){ cat("reorganize annotation ...", "\n") }
	if (ncpus == 1){ annotation <- t(sapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, GEann, chr, bpstart, bpend, simplify=TRUE)) }
	if (ncpus > 1){ annotation <- t(sfSapply(featuresAndWeights, .slhFuncWeightedCGHcallExpressionSetAnnotation, GEann, chr, bpstart, bpend, simplify=TRUE)) }
	colnames(annotation) <- c("Chromosome", "Start", "End")
	if (ncpus == 1){ rownames(annotation) <- as.character(sapply(featuresAndWeights, .slhFuncCombineFeatureNames, rownames(GEann), simplify=TRUE)) }
	if (ncpus > 1){ rownames(annotation) <- as.character(sfSapply(featuresAndWeights, .slhFuncCombineFeatureNames, rownames(GEann), simplify=TRUE)) }
	if (length(rownames(annotation)) != length(unique(rownames(annotation)))){
		warning("modified feature names are returned")
		newRowNames <- rownames(annotation)
		slh <- rep(1, length(newRowNames))
		for (i in 1:length(unique(newRowNames))){
			ids <- which(newRowNames == unique(newRowNames)[i])
			slh[ids] <- 1:length(ids)
		}
		newRowNames <- paste(newRowNames, slh, sep="__")		
		rownames(annotation) <- newRowNames
	}
	fInfo <- new("AnnotatedDataFrame", data=data.frame(annotation), varMetadata=data.frame(labelDescription=c("Chromosome", "Start", "End")))
	if (ncpus > 1){ sfStop() }
	rownames(exprMat) <- rownames(annotation)

	# make ExpressionSet object
	if (verbose){ cat("merge into ExpressionSet object ...", "\n") }
	return(new("ExpressionSet", exprs=exprMat, featureData=fInfo))
}






matchCGHcall2ExpressionSet <- function(CNdata, GEdata, CNchr, CNbpstart, CNbpend, GEchr, GEbpstart, GEbpend, method="distance", reference=1, ncpus=1, verbose=TRUE){
   	###########################################################################################
	# function that matches features from cghCall and ExpressionSet-objects on genomic location
	###########################################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	if ( !(as.character(class(GEdata)) == "ExpressionSet" | as.character(class(GEdata)) == "eSet") ){ stop("GEdata not of class ExpressionSet.") }
	if ( !(method %in% c("distance", "overlap", "overlapPlus")) ){ stop("method wrongly specified.") }

	if (!(as.character(class(GEchr))=="numeric" | as.character(class(GEchr))=="integer" )){ stop("GEchr of wrong class.") }
	if (length(GEchr) != 1){ stop("GEchr of wrong length.") }
	if ( !(GEchr >= 1 | GEchr <= ncol(fData(GEdata))) ){ stop("GEchr exceeds number of columns in featureData.") }

	if (!(as.character(class(GEbpstart))=="numeric" | as.character(class(GEbpstart))=="integer" )){ stop("GEbpstart of wrong class.") }
	if (length(GEbpstart) != 1){ stop("GEbpstart of wrong length.") }
	if ( !(GEbpstart >= 1 | GEbpstart <= ncol(fData(GEdata))) ){ stop("GEbpstart exceeds number of columns in featureData.") }

	if (!(as.character(class(GEbpend))=="numeric" | as.character(class(GEbpend))=="integer" )){ stop("GEbpend of wrong class.") }
	if (length(GEbpend) != 1){ stop("GEbpend of wrong length.") }
	if ( !(GEbpend >= 1 | GEbpend <= ncol(fData(GEdata))) ){ stop("GEbpend exceeds number of columns in featureData.") }

	if (!(as.character(class(CNchr))=="numeric" | as.character(class(CNchr))=="integer" )){ stop("CNchr of wrong class.") }
	if (length(CNchr) != 1){ stop("CNchr of wrong length.") }
	if ( !(CNchr >= 1 | CNchr <= ncol(fData(CNdata))) ){ stop("CNchr exceeds number of columns in featureData.") }

	if (!(as.character(class(CNbpstart))=="numeric" | as.character(class(CNbpstart))=="integer" )){ stop("CNbpstart of wrong class.") }
	if (length(CNbpstart) != 1){ stop("CNbpstart of wrong length.") }
	if ( !(CNbpstart >= 1 | CNbpstart <= ncol(fData(CNdata))) ){ stop("GEbpstart exceeds number of columns in featureData.") }

	if (!(as.character(class(CNbpend))=="numeric" | as.character(class(CNbpend))=="integer" )){ stop("CNbpend of wrong class.") }
	if (length(CNbpend) != 1){ stop("CNbpend of wrong length.") }
	if ( !(CNbpend >= 1 | CNbpend <= ncol(fData(CNdata))) ){ stop("CNbpend exceeds number of columns in featureData.") }

	# extract annotation info
	CNann <- cbind(fData(CNdata)[, c(CNchr, CNbpstart, CNbpend), drop=FALSE], c(1:nrow(fData(CNdata))))
	GEann <- cbind(fData(GEdata)[, c(GEchr, GEbpstart, GEbpend), drop=FALSE], c(1:nrow(fData(GEdata))))

	# in case columns are of wrong class ....
	for (i in 1:4){
		if (is.factor(CNann[,i])){ CNann[,i] <- as.numeric(levels(CNann[,i]))[CNann[,i]] }
		if (is.factor(GEann[,i])){ GEann[,i] <- as.numeric(levels(GEann[,i]))[GEann[,i]] }			
		if (is.character(CNann[,i])){ CNann[,i] <- as.numeric(CNann[,i]) }
		if (is.character(GEann[,i])){ GEann[,i] <- as.numeric(GEann[,i]) }			
	}

	# ensure data from gene expression platform is genomicly ordered.
	if (verbose){ cat("order both data sets genomically...", "\n") }
	GEdata <- ExpressionSet2order(GEdata, GEchr, GEbpstart, verbose=verbose)
	GEann <- GEann[order(GEann[,1], GEann[,2]), ]

	# ensure data from copy number platform is genomicly ordered.	
	CNdata <- cghCall2order(CNdata, CNchr, CNbpstart, verbose=verbose)
	CNann <- CNann[order(CNann[,1], CNann[,2]), ]

	# filter genes with complete annotation
	if (verbose){ cat("select probes with complete annotation...", "\n") }
	GEfeaturesKept <- which(apply(GEann, 1, function(x){ all(is.na(x) == FALSE) }))
	if (length(GEfeaturesKept) < 1){ stop("no annotation info for GE-array features provided") }
	GEann <- GEann[GEfeaturesKept, , drop=FALSE]
	CNfeaturesKept <- which(apply(CNann, 1, function(x){ all(is.na(x) == FALSE) }))
	if (length(CNfeaturesKept) < 1){ stop("no annotation info for CN-array features provided") }
	CNann <- CNann[CNfeaturesKept, , drop=FALSE]

	# chromosomes with genes on the expression array
	if (verbose){ cat("select probes on chromosomes present in both data sets...", "\n") }
	GEchrs <- sort(unique(GEann[,1]))
	CNfeaturesKept.2 <- which(CNann[,1] %in% GEchrs)
	CNann <- CNann[CNfeaturesKept.2, , drop=FALSE]
	CNfeaturesKept <- CNfeaturesKept[CNfeaturesKept.2]

	# chromosomes with features on the CN array
	CNchrs <- sort(unique(CNann[,1]))
	GEfeaturesKept.2 <- which(GEann[,1] %in% CNchrs)
	GEann <- GEann[GEfeaturesKept.2, , drop=FALSE]
	GEfeaturesKept <- GEfeaturesKept[GEfeaturesKept.2]
	
	# match by percentage of overlap
	if (verbose){ cat("start actual matching procedure...", "\n") }
	if (ncpus > 1){
		sfInit(parallel=TRUE, cpus=ncpus)
		sfLibrary(sigaR, verbose=FALSE)
	}
	if (method == "overlap"){
		# find features with biggest overlap and remove features without any overlap
		if (ncpus == 1){ CNGEmatched.features <- apply(GEann[,-4], 1, .overlapPercentage_CGHcall2ExpressionSet, CNann[,-4], reference) }
		if (ncpus > 1){ CNGEmatched.features <- sfApply(GEann[,-4], 1, .overlapPercentage_CGHcall2ExpressionSet, CNann[,-4], reference) }
		no.match.found <- (is.na(CNGEmatched.features) == FALSE)
		GEfeaturesKept <- GEfeaturesKept[no.match.found]		
		CNGEmatched.features <- CNGEmatched.features[no.match.found]
		CNfeaturesKept <- CNfeaturesKept[CNGEmatched.features]
	}

	# match by percentage of overlap with some interpolation
	if (method == "overlapPlus"){
			# find features with biggest overlap and remove features without any overlap
			if (ncpus == 1){ CNGEmatched.features <- apply(GEann[,-4], 1, .overlapPlusInterpolation_CGHcall2ExpressionSet, CNann[,-4], segmented(CNdata)[CNfeaturesKept, , drop=FALSE], reference) }
			if (ncpus > 1){ CNGEmatched.features <- sfApply(GEann[,-4], 1, .overlapPlusInterpolation_CGHcall2ExpressionSet, CNann[,-4], segmented(CNdata)[CNfeaturesKept, , drop=FALSE], reference) }
			no.match.found <- (is.na(CNGEmatched.features) == FALSE)
			GEfeaturesKept <- GEfeaturesKept[no.match.found]		
			CNGEmatched.features <- CNGEmatched.features[no.match.found]
			CNfeaturesKept <- CNfeaturesKept[CNGEmatched.features]
	}

	# match by distance of mid base pair
	if (method == "distance"){
		# determine mid base pair
		CNannTemp <- cbind(CNann[,1], apply(CNann[,2:3,drop=FALSE], 1, mean, na.rm=TRUE))
		GEannTemp <- cbind(GEann[,1], apply(GEann[,2:3,drop=FALSE], 1, mean, na.rm=TRUE))
		
		# per chromosome find CN-features closed to GE-features
		if (ncpus == 1){ CNGEmatched.features <- apply(GEannTemp, 1, .minDistCghFeature_CGHcall2ExpressionSet, CNannTemp) }
		if (ncpus > 1){ CNGEmatched.features <- sfApply(GEannTemp, 1, .minDistCghFeature_CGHcall2ExpressionSet, CNannTemp) }
		CNfeaturesKept <- CNfeaturesKept[CNGEmatched.features]
	}
	if (ncpus > 1){ sfStop() }

	# in the unlikely event that some genes could not be matched: remove
	if (sum(is.na(CNfeaturesKept)) > 0){
		GEfeaturesKept <- GEfeaturesKept[is.na(CNfeaturesKept) == FALSE]		
		CNfeaturesKept <- CNfeaturesKept[is.na(CNfeaturesKept) == FALSE]
	}

	# in the unlikely event that nothing could be matched: stop
	if ((length(CNfeaturesKept) == 0) | (length(GEfeaturesKept) == 0)){
		stop("No features could be matched. Check format of annotation information.")
	}
	matchedFeatures <- cbind(CNann[CNfeaturesKept,4], GEann[GEfeaturesKept,4])
	return(matchedFeatures)
}





uniqGenomicInfo <- function(chr, bpstart, bpend, verbose=FALSE){
   	############################################################################
	# function shrinks ExpressionSet-object to a weighted subset of the features
	############################################################################

	# check input
	if (verbose){ cat("perform input checks ...", "\n") }
	if (!(as.character(class(chr))=="numeric" | as.character(class(chr))=="integer" )){ stop("chr of wrong class.") }
	if (!(as.character(class(bpstart))=="numeric" | as.character(class(bpstart))=="integer" )){ stop("bpstart of wrong class.") }
	if (!(as.character(class(bpend))=="numeric" | as.character(class(bpend))=="integer" )){ stop("bpend of wrong class.") }
	if (length(chr) != length(bpstart) | length(chr) != length(bpend) | length(bpstart) != length(bpend)){ stop("chr, bpstart and bpend do not have identical length.") }

	# create identifier of unique genomic info's
	if (verbose){ cat("find features with identical genomic info ...", "\n") }
	ann <- cbind(chr, bpstart, bpend)[order(chr, bpstart, bpend),]
	uniques <- apply(cbind(ann[-nrow(ann),], ann[-1,]), 1, function(Z){ !all(Z[1:3] == Z[4:6]) })
	uniques <- c(1, cumsum(uniques) + 1)

	# put features with identical genomic info in list object
	ann <- cbind(1:nrow(ann), ann)
	colnames(ann)[1] <- "featureNo"
	slhFunc <- function(Z, uniques, ann){ ann[which(uniques==Z), , drop=FALSE] }
	uniqueGenomicInfo <- lapply(1:max(uniques), slhFunc, uniques, ann)
	return(uniqueGenomicInfo)
}







matchAnn2Ann <- function(chr1, bpstart1, bpend1, chr2, bpstart2, bpend2, method="distance", maxDist=10000, minPerc=0, reference=1, ncpus=1, verbose=TRUE){
   	###########################################################################################
	# function that matches annotation of to sets of features on genomic location
	###########################################################################################

	# check input
	if (verbose){ cat("perform input checks ...", "\n") }
	if (as.character(class(chr1)) != "numeric" & as.character(class(chr1)) != "integer"){ stop("chr1 of wrong class.") }
	if (as.character(class(bpstart1)) != "numeric" & as.character(class(bpstart1)) != "integer"){ stop("bpstart1 of wrong class.") }
	if (as.character(class(bpend1)) != "numeric" & as.character(class(bpend1)) != "integer"){ stop("bpend1 of wrong class.") }
	if ((length(chr1) != length(bpstart1)) & (length(chr1) != length(bpend1)) & (length(bpstart1) != length(bpend1))){ stop("length of chr1, bpstart1, and bpend1 do not match") }
	if (any(bpstart1 > bpend1)){ stop("inconsistency in bpstart1 and bpend1: for some feature bpstart1 > bpend1") }

	if (as.character(class(chr2)) != "numeric" & as.character(class(chr2)) != "integer"){ stop("chr2 of wrong class.") }
	if (as.character(class(bpstart2)) != "numeric" & as.character(class(bpstart2)) != "integer"){ stop("bpstart2 of wrong class.") }
	if (as.character(class(bpend2)) != "numeric" & as.character(class(bpend2)) != "integer"){ stop("bpend2 of wrong class.") }
	if ((length(chr2) != length(bpstart2)) & (length(chr2) != length(bpend2)) & (length(bpstart2) != length(bpend2))){ stop("length of chr2, bpstart2, and bpend2 do not match") }
	if (any(bpstart2 > bpend2)){ stop("inconsistency in bpstart2 and bpend2: for some feature bpstart2 > bpend2") }

	if ( !(method %in% c("distance", "overlap")) ){ stop("method wrongly specified.") }
	if (as.character(class(maxDist)) != "numeric" & as.character(class(maxDist)) != "integer"){ stop("maxDist of wrong class.") }
	if (as.character(class(minPerc)) != "numeric" & as.character(class(minPerc)) != "integer"){ stop("minPerc of wrong class.") }
	if (as.character(class(ncpus)) != "numeric" & as.character(class(ncpus)) != "integer"){ stop("ncpus of wrong class.") }
	if (maxDist < 1){ stop("maxDist wrongly specified.") }
	if (minPerc < 0 & minPerc > 1){ stop("minPerc wrongly specified.") }

	# extract annotation info
	ann1 <- cbind(chr1, bpstart1, bpend1, 1:length(chr1))
	colnames(ann1)[4] <- "probeNo1"
	ann2 <- cbind(chr2, bpstart2, bpend2, 1:length(chr2))
	colnames(ann2)[4] <- "probeNo2"
	# rm(chr1, chr2, bpstart1, bpstart2, bpend1, bpend2)

	# in case columns are of wrong class ....
	for (i in 1:3){
		if (is.factor(ann1[,i])){ ann1[,i] <- as.numeric(levels(ann1[,i]))[ann1[,i]] }
		if (is.factor(ann2[,i])){ ann2[,i] <- as.numeric(levels(ann2[,i]))[ann2[,i]] }			
		if (is.character(ann1[,i])){ ann1[,i] <- as.numeric(ann1[,i]) }
		if (is.character(ann2[,i])){ ann2[,i] <- as.numeric(ann2[,i]) }			
	}

	# ensure data from gene expression platform is genomicly ordered.
	if (verbose){ cat("order both data sets genomically ...", "\n")	}
	ann1 <- ann1[order(ann1[,1], ann1[,2]), ]
	ann2 <- ann2[order(ann2[,1], ann2[,2]), ]

	# filter features with complete annotation for platform 1
	if (verbose){ cat("select features with complete annotation ...", "\n")	}
	featuresKept1 <- which(apply(ann1, 1, function(x){ all(is.na(x) == FALSE) }))
	if (length(featuresKept1) < 1){ stop("no annotation info for features of platform 1 provided") }
	ann1 <- ann1[featuresKept1, , drop=FALSE]
	ann1removed <- ann1[-featuresKept1, , drop=FALSE]

	# filter features with complete annotation for platform 2
	featuresKept2 <- which(apply(ann2, 1, function(x){ all(is.na(x) == FALSE) }))
	if (length(featuresKept2) < 1){ stop("no annotation info for features of platform 2 provided") }
	ann2removed <- ann2[-featuresKept2, , drop=FALSE]
	ann2 <- ann2[featuresKept2, , drop=FALSE]

	# chromosomes with features on platform 2
	if (verbose){ cat("select features on chromosomes present in both data sets ...", "\n")	}
	chrs2 <- sort(unique(ann2[,1]))
	featuresKept1 <- which(ann1[,1] %in% chrs2)
	ann1removed <- rbind(ann1removed, ann1[-featuresKept1, , drop=FALSE])
	ann1 <- ann1[featuresKept1, , drop=FALSE]

	# chromosomes with features on platform 1
	chrs1 <- sort(unique(ann1[,1]))
	featuresKept2 <- which(ann2[,1] %in% chrs1)
	ann2removed <- rbind(ann2removed, ann2[-featuresKept2, , drop=FALSE])
	ann2 <- ann2[featuresKept2, , drop=FALSE]

	# match by distance of mid base pair
	if (verbose){ cat("start actual matching procedure ...", "\n") }
	if (ncpus > 1){	
		sfInit(parallel=TRUE, cpus=ncpus)
		sfLibrary(sigaR, verbose=FALSE)
	}
	if (method == "distance"){
		# determine mid base pair
		ann1temp <- cbind(ann1[,1], apply(ann1[,2:3,drop=FALSE], 1, mean, na.rm=TRUE), ann1[,4])
		ann2temp <- cbind(ann2[,1], apply(ann2[,2:3,drop=FALSE], 1, mean, na.rm=TRUE), ann2[,4])
		
		# per chromosome find CN-features closed to GE-features
		if (ncpus == 1){ matchedFeatures <- apply(ann2temp, 1, .minDistFeature, ann1temp, maxDist) }
		if (ncpus > 1){	matchedFeatures <- sfApply(ann2temp, 1, .minDistFeature, ann1temp, maxDist) }
	}

	# match by percentage of overlap
	if (method == "overlap"){
		# find features with biggest overlap and remove features without any overlap
		if (ncpus == 1){ matchedFeatures <- apply(ann2, 1, .overlapPercentage, ann1, minPerc, reference) }
		if (ncpus > 1){ matchedFeatures <- sfApply(ann2, 1, .overlapPercentage, ann1, minPerc, reference) }
	}
	if (ncpus > 1){ sfStop() }

	if (!is.list(matchedFeatures)){ matchedFeatures <- .matrixRows2list(matchedFeatures) }

	# reformatting
	if (verbose){ cat("reformat matching results ...", "\n") }
	names(matchedFeatures) <- as.character(ann2[,4])
	idsNull <- as.numeric(which(sapply(matchedFeatures, is.null, simplify=TRUE)))
	if (length(idsNull) > 1){ matchedFeatures <- matchedFeatures[-idsNull] }

	return(matchedFeatures)
}



getSegFeatures <- function(featureNo, CNdata, verbose=TRUE){
	###########################################################
	# function that determine features with similar CN signature as featureNo 
	###########################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if (as.character(class(CNdata)) != "cghCall"){ stop("CNdata not of class cghCall.") }
	if (!(as.character(class(featureNo))=="numeric" | as.character(class(featureNo))=="integer" )){ stop("featureSubset of wrong class.") }
	if ( length(featureNo) != 1 ){ stop("featureNo wrongly specified.") }
	if ( !(featureNo >= 1 | featureNo <= nrow(fData(CNdata))) ){ stop("featureNo out of range.") }

	# get features
	CNfeaturesKept <- which(chromosomes(CNdata) == chromosomes(CNdata)[featureNo])
	segSignature <- segmented(CNdata)[featureNo, , drop=FALSE]
	CNfeaturesKept <- CNfeaturesKept[apply(segmented(CNdata)[CNfeaturesKept, , drop=FALSE], 1, function(Y, Z){ all(as.numeric(Y) == as.numeric(Z)) }, Z=segSignature)]

	# in case possible gap between features with same signature
	idTemp <- which(CNfeaturesKept == featureNo)
	if (idTemp == 1){
		checkSeq <- c(0:(length(CNfeaturesKept)-idTemp))
	} else {
		checkSeq <- c(-c((idTemp-1):1), c(0:(length(CNfeaturesKept)-idTemp)))
	}
	CNfeaturesKept <- CNfeaturesKept[which(checkSeq == (CNfeaturesKept - featureNo))]
	return(CNfeaturesKept)
}



ExpressionSet2subset <- function(GEdata, featureSubset, verbose=TRUE){
   	############################################################################
	# function shrinks ExpressionSet-object to a subset of the features
	############################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if ( !(as.character(class(GEdata)) == "ExpressionSet" | as.character(class(GEdata)) == "eSet") ){ stop("GEdata not of class ExpressionSet.") }
	if (!(as.character(class(featureSubset))=="numeric" | as.character(class(featureSubset))=="integer" )){ stop("featureSubset of wrong class.") }
	if ( !(length(featureSubset) >= 1) ){ stop("featureSubset contains no row numbers.") }
	if ( !(min(featureSubset) >= 1 | max(featureSubset) <= nrow(fData(GEdata))) ){ stop("featureSubset contains illegal row numbers.") }

	# subsetting
	fData(GEdata) <- fData(GEdata)[featureSubset, , drop=FALSE]	
	exprs(GEdata) <- exprs(GEdata)[featureSubset, , drop=FALSE]	
	return(GEdata)
}



ExpressionSet2order <- function(GEdata, chr, bpstart, verbose=TRUE){
   	###########################################################################################
	# function that orders ExpressionSet-object genomically
	###########################################################################################

	# check input
	if (verbose){ cat("perform input checks...", "\n") }
	if ( !(class(GEdata) == "ExpressionSet" | class(GEdata) == "eSet") ){ stop("GEdata not of class ExpressionSet.") }
	if (!(class(chr)=="numeric" | class(chr)=="integer" )){ stop("chr of wrong class.") }
	if (length(chr) != 1){ stop("chr of wrong length.") }
	if ( !(chr >= 1 | chr <= ncol(fData(GEdata))) ){ stop("chr exceeds number of columns in featureData.") }
	if (!(class(bpstart)=="numeric" | class(bpstart)=="integer" )){ stop("bpstart of wrong class.") }
	if (length(bpstart) != 1){ stop("bpstart of wrong length.") }
	if ( !(bpstart >= 1 | bpstart <= ncol(fData(GEdata))) ){ stop("bpstart exceeds number of columns in featureData.") }

	# extract annotation data for ordering
	GEann <- fData(GEdata)[, c(chr, bpstart), drop=FALSE]

	# in case columns are of wrong class ....
	for (i in 1:2){
		if (is.factor(GEann[,i])){ GEann[,i] <- as.numeric(levels(GEann[,i]))[GEann[,i]] }			
		if (is.character(GEann[,i])){ GEann[,i] <- as.numeric(GEann[,i]) }			
	}

	# order ExpressionSet-object genomically
	exprs(GEdata) <- exprs(GEdata)[order(GEann[,1], GEann[,2]), ]
	fData(GEdata) <- fData(GEdata)[order(GEann[,1], GEann[,2]), ]
	return(GEdata)
}


