\name{ExpressionSet2weightedSubset}
\alias{ExpressionSet2weightedSubset}
\title{
Weighted subsetting ExpressionSet-objects.
}
\description{
Limit an \code{\link{ExpressionSet}} object to a subset of its features, using weighted averaging of the expression signal.
}
\usage{
ExpressionSet2weightedSubset(GEdata, featuresAndWeights, chr, bpstart, bpend, ncpus = 1, verbose=TRUE)
}
\arguments{
  \item{GEdata}{ Object of class \code{\link{ExpressionSet}}. }
  \item{featuresAndWeights}{ Object of class \code{list}. Each list item is a \code{matrix}. The first column of this \code{matrix} contains the row numbers of features to be maintained in the \code{\link{ExpressionSet}}-object. The second column contains the weights of each features, to be used in the calculation of the weighted average gene expression signal. }
  \item{chr}{ Column in the slot \code{featureData} of the \code{\link{ExpressionSet}}-object specifying the chromosome information of the features. }
  \item{bpstart}{ Column in the slot \code{featureData} of the \code{\link{ExpressionSet}}-object specifying the start basepair information of the features. }
  \item{bpend}{ Column in the slot \code{featureData} of the \code{\link{ExpressionSet}}-object specifying the end basepair information of the features. }
  \item{ncpus}{ Number of cpus to be used in computations. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\details{
Annotation information of features with multiplicity larger than one is compressed as follows. It is assumed that all features map to the same chromosome, leaving no ambiguity. 
The start base pair of the "new" feature is the smallest start base pair of features from which it has been formed.
The end base pair of the "new" feature is the largest end base pair of features from which it has been formed.
}
\value{ Object of class \code{\link{ExpressionSet}}, restricted to the specified subset of features. }
\references{	
Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2011), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{BMC Bioinformatics}, 13:80. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\note{
This is a more intricate version of the \code{ExpressionSet2subset} function. They exists parallel because this function is much slower than its counterpart.
}
\section{Warning}{The \code{phenoData}, \code{experimentData}, and other slots of the \code{ExpressionSet}-object are currently not passed on to the subsetted object. }
\seealso{\code{ExpressionSet2subset}}
\examples{
# load data
data(pollackGE16)

# extract genomic information from ExpressionSet-object
chr <- fData(pollackGE16)[,1]
bpstart <- fData(pollackGE16)[,2]
bpend <- fData(pollackGE16)[,3]

# find unique genomic locations
uniqInfo <- uniqGenomicInfo(chr, bpstart, bpend, verbose = FALSE) 

# subset ExpressionSet-object to features with unique genomic locations
pollackGE16 <- ExpressionSet2weightedSubset(pollackGE16, uniqInfo, 1, 2, 3) 
}
