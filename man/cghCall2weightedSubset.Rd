\name{cghCall2weightedSubset}
\alias{cghCall2weightedSubset}
\title{
Weighted subsetting cghCall-objects.
}
\description{
Limit an \code{\link{cghCall}} object to a subset of its features, using weighted averaging of the copy number signal.
}
\usage{
cghCall2weightedSubset(CNdata, featuresAndWeights, chr, bpstart, bpend, ncpus = 1, verbose=TRUE)
}
\arguments{
  \item{CNdata}{ Object of class \code{\link{cghCall}}. }
  \item{featuresAndWeights}{ Object of class \code{list}. Each list item is a \code{matrix}. The first column of this \code{matrix} contains the row numbers of features to be maintained in the \code{\link{cghCall}}-object. The second column contains the weights of each features, to be used in the calculation of the weighted average copy number signal. }
  \item{chr}{ Column in the slot \code{featureData} of the \code{\link{cghCall}}-object specifying the chromosome information of the features. }
  \item{bpstart}{ Column in the slot \code{featureData} of the \code{\link{cghCall}}-object specifying the start basepair information of the features. }
  \item{bpend}{ Column in the slot \code{featureData} of the \code{\link{cghCall}}-object specifying the end basepair information of the features. }
  \item{ncpus}{ Number of cpus to be used in computations. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{ Object of class \code{\link{cghCall}}, restricted to the specified subset of features. }
\references{	
Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2011), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{submitted for publication}. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\note{
This is a more intricate version of the \code{cghCall2subset} function. They exists parallel because this function is (much) slower than its counterpart.
}

\section{Warning}{The \code{phenoData}, \code{experimentData}, and other slots of the \code{cghCall}-object are currently not passed on to the subsetted object. }

\seealso{\code{cghCall2subset}}
\examples{
# load data
data(pollackCN16)

# extract genomic information from ExpressionSet-object
chr <- fData(pollackCN16)[,1]
bpstart <- fData(pollackCN16)[,2]
bpend <- fData(pollackCN16)[,3]

# find unique genomic locations
uniqInfo <- uniqGenomicInfo(chr, bpstart, bpend, verbose = FALSE) 

# subset cghCall-object to features with unique genomic locations
pollackCN16 <- cghCall2weightedSubset(pollackCN16, uniqInfo, 1, 2, 3) 
}
