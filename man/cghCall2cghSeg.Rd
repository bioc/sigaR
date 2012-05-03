\name{cghCall2cghSeg}
\alias{cghCall2cghSeg}
\title{
Genomic ordering of cghSeg-objects.
}
\description{
Transforms a \code{\link{cghCall}}-object to a \code{\link{cghSeg}}-object, by removing the slots present in the former but not in the latter.
}
\usage{
cghCall2cghSeg(CNdata, verbose=TRUE)
}
\arguments{
  \item{CNdata}{ Object of class \code{\link{cghCall}}. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{\link{cghSeg}}.
}
\references{
Van de Wiel, M.A., Kim, K.I., Vosse, S.J., Van Wieringen, W.N., Wilting, S.M. , Ylstra, B. (2007), "CGHcall: an algorithm for calling aberrations for multiple array CGH tumor profiles", Bioinformatics, 23, 892-894.
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
\code{\link{cghCall}}, \code{\link{cghSeg}}.
}
\examples{
# load data
data(pollackCN16) 

# reduce the cghCall-object to a cghSeg-object
pollackCN16seg <- cghCall2cghSeg(pollackCN16) 
}
