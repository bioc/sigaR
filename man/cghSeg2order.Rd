\name{cghSeg2order}
\alias{cghSeg2order}
\title{
Genomic ordering of cghSeg-objects.
}
\description{
Orders the features within a \code{\link{cghSeg}}-object in accordance with their genomic order.
}
\usage{
cghSeg2order(CNdata, chr, bpstart, verbose=TRUE)
}
\arguments{
  \item{CNdata}{ Object of class \code{\link{cghSeg}}. }
  \item{chr}{ Column in the slot \code{featureData} of the \code{\link{cghSeg}}-object specifying the chromosome information of the features. }
  \item{bpstart}{ Column in the slot \code{featureData} of the \code{\link{cghSeg}}-object specifying the start basepair information of the features. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{\link{cghSeg}}, now genomically ordered.
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

# order the copy number data genomically
pollackCN16 <- cghSeg2order(pollackCN16, 1, 2) 
}
