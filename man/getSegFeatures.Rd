\name{getSegFeatures}
\alias{getSegFeatures}
\title{
Identical signature features selection from cghCall-object.
}
\description{
Given an example, selects features (contiguous to the example) with the same signature (as the example) across samples from an \code{\link{cghCall}}-object.
}
\usage{
getSegFeatures(featureNo, CNdata, verbose=TRUE)
}
\arguments{
  \item{featureNo}{ Row number of example feature. }
  \item{CNdata}{ Object of class \code{\link{cghCall}}. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\value{
Object of class \code{numeric}, containing the row numbers of those contiguous features with the same segmented log2-ratio signatures as \code{featureNo} across samples.
}

\references{
Van Wieringen, W.N., Berkhof, J., Van de Wiel, M.A. (2010), "A random coefficients model for regional co-expression associated with DNA copy number", \emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Issue1, Article 25, 1-28. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
\code{\link{cghCall}}, \code{\link{RCMestimation}}.
}
\examples{
# load data
data(pollackCN16) 

# feature of interest
featureNo <- 7

# extract all features with identical copy number signature (over the samples)
getSegFeatures(featureNo, pollackCN16) 
}
