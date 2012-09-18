\name{nBreakpoints}
\alias{nBreakpoints}
\title{
Number of breakpoints
}
\description{
The number of samples with at least one breakpoint is calculated for each transcipt.
}
\usage{
nBreakpoints(featuresAndWeights, CNdata)
}
\arguments{
  \item{featuresAndWeights}{ Object of class \code{list}. Each list item is a \code{matrix}. The first column of this \code{matrix} contains the row numbers of features to be maintained in subsetting of the \code{\link{cghCall}}-object. 
	The second column contains the weights of each features, to be used in the calculation of the weighted average copy number signal. }
  \item{CNdata}{ Object of class \code{\link{cghCall}} }
}
\details{
For each item of the object \code{featuresAndWeights} the segmented data from the \code{cghCall}-object is used to determine whether a sample exhibits a breakpoint for this transcript.
}
\value{
Object of class \code{numeric} containing the number of samples with at least one breakpoint. It is of the same length as the \code{featuresAndWeights}-object.
}
\references{
Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2012), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{submitted for publication}. 
}
\author{
Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} 
}
\seealso{
\code{\link{matchAnn2Ann}}.
}
\examples{
# load data
data(pollackCN16)
data(pollackGE16)

# extract genomic information from cghCall-object
chr1 <- fData(pollackCN16)[,1]
bpstart1 <- fData(pollackCN16)[,2]
bpend1 <- fData(pollackCN16)[,3]

# extract genomic information from ExpressionSet-object
chr2 <- fData(pollackGE16)[,1]
bpstart2 <- fData(pollackGE16)[,2]
bpend2 <- fData(pollackGE16)[,3]

# match features from both platforms
matchedIDs <- matchAnn2Ann(chr1, bpstart1, bpend1, chr2, bpstart2, bpend2, method = "distance", maxDist = 10000) 

# extract ids for object subsetting
matchedIDsCN <- lapply(matchedIDs, function(Z){ return(Z[, -1, drop=FALSE]) }) 

# calculate the number of breakpoints
nBreakpoints(matchedIDsCN, pollackCN16) 
}
