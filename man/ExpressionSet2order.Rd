\name{ExpressionSet2order}
\alias{ExpressionSet2order}
\title{
Genomic ordering of ExpressionSet-objects.
}
\description{
Orders the features within a \code{\link{ExpressionSet}}-object in accordance with their genomic order.
}
\usage{
ExpressionSet2order(GEdata, chr, bpstart, verbose=TRUE)
}
\arguments{
  \item{GEdata}{ Object of class \code{\link{ExpressionSet}}. }
  \item{chr}{ Column in the slot \code{featureData} of the \code{\link{ExpressionSet}}-object specifying the chromosome information of the features. }
  \item{bpstart}{ Column in the slot \code{featureData} of the \code{\link{ExpressionSet}}-object specifying the start basepair information of the features. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{\link{ExpressionSet}}, now genomically ordered.
}

\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
\code{\link{ExpressionSet}}.
}
\examples{
# load data
data(pollackGE16) 

# order the copy number data genomically
pollackGE16 <- ExpressionSet2order(pollackGE16, 1, 2) 
}
