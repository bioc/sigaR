\name{merge2cghCalls}
\alias{merge2cghCalls}
\title{
Merge two cghCall-objects into one cghCall-object
}
\description{
Merge two cghCall-objects into one cghCall-object.
}
\usage{
merge2cghCalls(CNdata1, CNdata2, verbose=TRUE)
}
\arguments{
  \item{CNdata1}{ Object of class \code{\link{cghCall}}. }
  \item{CNdata2}{ Object of class \code{\link{cghCall}}. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{\link{cghCall}}, restricted to the specified subset of features.
}
\details{
Data of the two objects is assumed to originate from the same samples, and are presented in the same order.

Only the experimental data and annotation information is inherited by the merged object. 
}
\references{
Van de Wiel, M.A., Kim, K.I., Vosse, S.J., Van Wieringen, W.N., Wilting, S.M. , Ylstra, B. (2007), "CGHcall: an algorithm for calling aberrations for multiple array CGH tumor profiles", Bioinformatics, 23, 892-894.

Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2012), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{BMC Bioinformatics}, 13:80. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
\code{\link{cghCall}}.
}
\examples{
# load data
data(pollackCN16) 

# create two cghCall-objects
ids1 <- sample(1:dim(pollackCN16)[1], 10) 
CNdata1 <- pollackCN16[ids1,] 
CNdata2 <- pollackCN16[-ids1,] 

# order the copy number data genomically
pollackCN16 <- merge2cghCalls(CNdata1, CNdata2) 
}

