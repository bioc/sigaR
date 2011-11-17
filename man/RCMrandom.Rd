\name{RCMrandom}
\alias{RCMrandom}
\title{
Random data from the random coefficients model.
}
\description{
The significance of hypotheses regarding parameters of the random coefficients model is assessed by means of the
parametric bootstrap. Hereto random data from the fitted model under the null hypothesis of interest are drawn.
This function provides.
}
\usage{
RCMrandom(object)
}
\arguments{
  \item{object}{ Object of class \code{rcmFit}.
}
}
\details{
Details on the type of random coefficients model from which data are drawn are specified in the reference below.
}
\value{
  A \code{matrix} of dimension (number of genes) times (number of samples).
}
\references{
Van Wieringen, W.N., Berkhof, J., Van de Wiel, M.A. (2010), "A random coefficients model for regional co-expression associated with DNA copy number", \emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Issue1, Article 25, 1-28. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{  \code{\link{RCMestimation}}, \code{rcmFit}. }
\examples{
# load data
data(pollackCN16)
data(pollackGE16)

# select features belonging to a region
ids <- getSegFeatures(20, pollackCN16)

# extract segmented log2 ratios of the region
X <- t(segmented(pollackCN16)[ids[1], , drop=FALSE])

# extract segmented log2 ratios of the region
Y <- exprs(pollackGE16)[ids,]

# center the expression data (row-wise)
Y <- t(Y - apply(Y, 1, mean))

# specify the linear constraint matrix
R <- matrix(1, nrow=1)

# fit the random coefficients model to the random data
RCMresults <- RCMestimation(Y, X, R)

# draw random data
Yrandom <- RCMrandom(RCMresults)
}
