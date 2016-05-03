\name{cghCall2maximumSubset}
\alias{cghCall2maximumSubset}
\title{
Maximum subsetting cghCall-objects.
}
\description{
Limit an \code{\link{cghCall}} object to a subset of its features, selecting those features with the most deviating copy number signal.
}
\usage{
cghCall2maximumSubset(CNdata, featuresAndWeights, chr, bpstart, 
bpend, ncpus = 1, verbose=TRUE)
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

\details{
Per entry of the \code{featuresAndWeights}-object and per sample the feature with the maximum absolute segmented DNA copy number signal is selected.
}
\value{ Object of class \code{\link{cghCall}}, restricted to the specified subset of features. }
\references{	
Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2012), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{BMC Bioinformatics}, 13:80. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
\code{\link{matchAnn2Ann}}
}
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
pollackCN16 <- cghCall2maximumSubset(pollackCN16, uniqInfo, 1, 2, 3) 
}
