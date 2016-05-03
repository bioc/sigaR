\name{RCMestimation}
\alias{RCMestimation}
\title{ Fitting of the random coefficients model. }
\description{
The parameters of the random coefficients model are estimated by means of the maximum likelihood method.
The implemented maximum likelihood procedure has been optimized with respect to computational efficiency and
memory usage.
}
\usage{
RCMestimation(Y, X, R, hypothesis = "H2", shrinkType = "none", 
estType = "normal", corType = "unif", maxNoIt = 100, 
minSuccDist = 0.005, verbose = FALSE)
}
\arguments{
  \item{Y}{ The \code{matrix} containing the (e.g., expression) data (number of columns equal to number of features, number of rows equal to number of samples). }
  \item{X}{ The design matrix (number of rows equal to number of samples, number of columns equal to number of covariates). }
  \item{R}{ The linear constraint matrix (number of columns equal to the number of covariates). }
  \item{hypothesis}{ The hypothesis under which the model is fitted: \code{H0} (H0  : R beta = 0 & tau2 = 0), \code{H1} (H1 : R beta >= 0 & tau2 = 0), \code{H2} (H2  : R beta >= 0 & tau2 >= 0). }
  \item{shrinkType}{ The type of shrinkage to be applied to the error variances: \code{none} (shrinkage parameter is set equal to zero: no shrinkage), \code{opt} (shrinkage parameter is chosen to minimize the mean squared error criterion) or \code{full} (shrinkage parameter is set equal to one).}
  \item{estType}{ Type of estimation, either \code{normal} (non-robust) or \code{robust}. }
  \item{corType}{ Correlation structure to be used, either \code{unif} or \code{ar1}. }
  \item{maxNoIt}{ Maximum number of iterations in the ML procedure. }
  \item{minSuccDist}{ Minimum distance between estimates of two successive iterations to be achieved.}
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\details{
Details on the type of random coefficients model that is actually fitted are specified in the reference below.
}
\value{
Object of class \code{rcmFit}.
}
\references{
Van Wieringen, W.N., Berkhof, J., Van de Wiel, M.A. (2010), "A random coefficients model for regional co-expression associated with DNA copy number", \emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Issue1, Article 25, 1-28. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\section{Note}{
In case a covariate for the intercept is included in the design matrix \code{X} we strongly recommend the center, per feature, the data around zero.
}

\seealso{
\code{\link{RCMrandom}}, \code{\link{RCMtest}}, \code{rcmTest}.
}
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
}
