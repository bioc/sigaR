\name{uniqGenomicInfo}
\alias{uniqGenomicInfo}
\title{
Unique genomic location information
}
\description{
Finds unique genomic location information.
}
\usage{
uniqGenomicInfo(chr, bpstart, bpend, verbose = FALSE)
}
\arguments{
  \item{chr}{ Object of class \code{numeric} containing chromosome information of features. }
  \item{bpstart}{ Object of class \code{numeric} containing start base pair information of features. Of same length as \code{chr}. }
  \item{bpend}{ Object of class \code{numeric} containing end base pair information of features. Of same length as \code{chr}. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
An object of class \code{list}. Each list item is a four-column matrix with the matched features information. The first column contains feature numbers of features with identical genomic location.
The second, third and fourth column contain the chromosome, start and end base pair information of the features (should be the same for each feature).
}
\references{	
Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2011), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{submitted for publication}. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{\code{ExpressionSet2weightedSubset}, \code{cghCall2weightedSubset}
}
\examples{
# load data
data(pollackGE16)

# extract genomic information from ExpressionSet-object
chr <- fData(pollackGE16)[,1]
bpstart <- fData(pollackGE16)[,2]
bpend <- fData(pollackGE16)[,3]

# find unique genomic locations
uniqInfo <- uniqGenomicInfo(chr, bpstart, bpend, verbose = FALSE) 
}
