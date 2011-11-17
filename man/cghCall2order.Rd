\name{cghCall2order}
\alias{cghCall2order}
\title{
Genomic ordering of cghCall-objects.
}
\description{
Orders the features within a \code{\link{cghCall}}-object in accordance with their genomic order.
}
\usage{
cghCall2order(CNdata, chr, bpstart, verbose=TRUE)
}
\arguments{
  \item{CNdata}{ Object of class \code{\link{cghCall}}. }
  \item{chr}{ Column in the slot \code{featureData} of the \code{\link{cghCall}}-object specifying the chromosome information of the features. }
  \item{bpstart}{ Column in the slot \code{featureData} of the \code{\link{cghCall}}-object specifying the start basepair information of the features. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{\link{cghCall}}, now genomically ordered.
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
pollackCN16 <- cghCall2order(pollackCN16, 1, 2) 
}
