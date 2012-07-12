\name{merge2ExpressionSets}
\alias{merge2ExpressionSets}
\title{
Merge two ExpressionSet-objects into one ExpressionSet-object
}
\description{
Merge two ExpressionSet-objects into one ExpressionSet-object
}
\usage{
merge2ExpressionSets(GEdata1, GEdata2, verbose=TRUE)
}
\arguments{
  \item{GEdata1}{ Object of class \code{\link{ExpressionSet}}. }
  \item{GEdata2}{ Object of class \code{\link{ExpressionSet}}. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{\link{ExpressionSet}}, restricted to the specified subset of features.
}
\details{
Data of the two objects is assumed to originate from the same samples, and are presented in the same order.

Only the experimental data and annotation information is inherited by the merged object. 
}
\references{
Van de Wiel, M.A., Kim, K.I., Vosse, S.J., Van Wieringen, W.N., Wilting, S.M. , Ylstra, B. (2007), "CGHcall: an algorithm for calling aberrations for multiple array CGH tumor profiles", Bioinformatics, 23, 892-894.

Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2011), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{BMC Bioinformatics}, 13:80. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
\code{\link{ExpressionSet}}.
}
\examples{
# load data
data(pollackGE16) 

# create two cghCall-objects
ids1 <- sample(1:dim(pollackGE16)[1], 10) 
GEdata1 <- pollackGE16[ids1,] 
GEdata2 <- pollackGE16[-ids1,] 

# order the copy number data genomically
pollackGE16 <- merge2ExpressionSets(GEdata1, GEdata2) 
}

