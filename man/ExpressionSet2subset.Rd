\name{ExpressionSet2subset}
\alias{ExpressionSet2subset}
\title{
Subsetting ExpressionSet-objects.
}
\description{
Limit an \code{\link{ExpressionSet}} object to a subset of its features.
}
\usage{
ExpressionSet2subset(GEdata, featureSubset, verbose=TRUE)
}
\arguments{
  \item{GEdata}{ Object of class \code{\link{ExpressionSet}}. }
  \item{featureSubset}{ Object of class \code{numeric}, containing the row numbers of features to be maintained in the \code{\link{ExpressionSet}}-object.}
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{\link{ExpressionSet}}, restricted to the specified subset of features.
}

\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
\code{\link{ExpressionSet}}.
}
\examples{
# load data
data(pollackGE16) 

# order the copy number data genomically
pollackGE16 <- ExpressionSet2subset(pollackGE16, c(1:50)) 
}
