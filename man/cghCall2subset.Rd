\name{cghCall2subset}
\alias{cghCall2subset}
\title{
Subsetting cghCall-objects.
}
\description{
Limit an \code{\link{cghCall}} object to a subset of its features.
}
\usage{
cghCall2subset(CNdata, featureSubset, verbose=TRUE)
}
\arguments{
  \item{CNdata}{ Object of class \code{\link{cghCall}}. }
  \item{featureSubset}{ Object of class \code{numeric}, containing the row numbers of features to be maintained in the \code{\link{cghCall}}-object. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{\link{cghCall}}, restricted to the specified subset of features.
}
\references{
Van de Wiel, M.A., Kim, K.I., Vosse, S.J., Van Wieringen, W.N., Wilting, S.M. , Ylstra, B. (2007), "CGHcall: an algorithm for calling aberrations for multiple array CGH tumor profiles", Bioinformatics, 23, 892-894.
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }


\seealso{
\code{\link{cghCall}}.
}
\examples{
# load data
data(pollackCN16) 

# order the copy number data genomically
pollackCN16 <- cghCall2subset(pollackCN16, c(1:50)) 
}

