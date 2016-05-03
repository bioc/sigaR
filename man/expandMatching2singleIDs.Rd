\name{expandMatching2singleIDs}
\alias{expandMatching2singleIDs}
\title{
Expand matching to single entries
}
\description{
In case a feature of platform 1 has been matched to multiple features of another platform, 
instead of averaging the data from these features, one may consider maintaining all features, each matched individually the feature of platform 1.
This function modifies the results from the matching function \code{matchAnn2Ann} to facilitate this. The result can than directly be used in 
the subsetting functions \code{cghCall2weightedSubset} and \code{ExpressionSet2weightedSubset}.
}
\usage{
expandMatching2singleIDs(matchedIDs)
}
\arguments{
  \item{matchedIDs}{ An object of class \code{list}, as returned by the \code{matchAnn2Ann}-function. }
}
\value{
An object of class \code{list}, similar to that returned by the \code{matchAnn2Ann}-function.
}
\references{
Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2012), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{BMC Bioinformatics}, 13:80. 
}
\author{
Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} 
}

\seealso{
\code{\link{matchAnn2Ann}}, \code{\link{cghCall2weightedSubset}}, \code{\link{ExpressionSet2weightedSubset}}.
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
matchedFeatures <- matchAnn2Ann(chr1, bpstart1, bpend1, chr2, 
	bpstart2, bpend2, method = "distance", maxDist = 10000) 

# expand
matchedFeatures <- expandMatching2singleIDs(matchedFeatures) 
}
