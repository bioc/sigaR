\name{cghSeg2subset}
\alias{cghSeg2subset}
\title{
Subsetting cghSeg-objects.
}
\description{
Limit an \code{\link{cghSeg}} object to a subset of its features.
}
\usage{
cghSeg2subset(CNdata, featureSubset, verbose=TRUE)
}
\arguments{
  \item{CNdata}{ Object of class \code{\link{cghSeg}}. }
  \item{featureSubset}{ Object of class \code{numeric}, containing the row numbers of features to be maintained in the \code{\link{cghSeg}}-object. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{\link{cghSeg}}, restricted to the specified subset of features.
}
\references{
Van de Wiel, M.A., Kim, K.I., Vosse, S.J., Van Wieringen, W.N., Wilting, S.M. , Ylstra, B. (2007), "CGHcall: an algorithm for calling aberrations for multiple array CGH tumor profiles", Bioinformatics, 23, 892-894.
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }


\seealso{
\code{\link{cghSeg}}.
}
\examples{
# load data
data(pollackCN16) 

# transform the cghCall-object to a cghSeg-object
pollackCN16 <- cghCall2cghSeg(pollackCN16)

# subset the copy number data
pollackCN16 <- cghSeg2subset(pollackCN16, c(1:50)) 
}

