\name{RCMtest}
\alias{RCMtest}
\title{
Hypothesis testing within the random coefficient model.
}
\description{
Function that evaluates various hypothesis within the random coefficients model via bootstrap resampling.
}
\usage{
RCMtest(Y, X, R, testType = "I", nBoot = 100, lowCiThres = 0.1, shrinkType = "none", estType = "normal", corType = "unif", maxNoIt = 100, minSuccDist = 0.005, returnNullDist = FALSE, ncpus=1, verbose = FALSE)
}
\arguments{
  \item{Y}{ The \code{matrix} containing the (e.g., expression) data (number of columns equal to number of features, number of rows equal to number of samples). }
  \item{X}{ The design matrix (number of rows equal to number of samples, number of columns equal to number of covariates). }
  \item{R}{ The linear constraint matrix (number of columns equal to the number of covariates). }
  \item{testType}{ The hypothesis to be tested: \code{I} (H0  : R beta = 0 & tau2 = 0) vs. (H2  : R beta >= 0 V tau2 >= 0), \code{II} (H0  : R beta = 0 & tau2 = 0) vs. (H1 : R beta >= 0 & tau2 = 0), \code{III} (H1 : R beta >= 0 & tau2 = 0) vs. (H2  : R beta >= 0 & tau2 >= 0). }
  \item{nBoot}{ Number of bootstraps. }
  \item{lowCiThres}{ A value between 0 and 1. Determines speed of efficient p-value calculation. If the probability of a p-value being below \code{lowCiThres} is smaller than 0.001 (read: the test is unlikely to become significant), bootstrapping is terminated and a p-value of 1.00 is reported. }
  \item{shrinkType}{ The type of shrinkage to be applied to the error variances: \code{none} (shrinkage parameter is set equal to zero: no shrinkage), \code{opt} (shrinkage parameter is chosen to minimize the mean squared error criterion) or \code{full} (shrinkage parameter is set equal to one).}
  \item{estType}{ Type of estimation, either \code{normal} (non-robust) or \code{robust}. }
  \item{corType}{ Correlation structure to be used, either \code{unif} or \code{ar1}. }
  \item{maxNoIt}{ Maximum number of iterations in the ML procedure. }
  \item{minSuccDist}{ Minimum distance between estimates of two successive iterations to be achieved.}
  \item{returnNullDist}{ Logical indicator: should the null distribution be returned? }
  \item{ncpus}{ Number of cpus used for the bootstrap. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\details{
Details on the type of random coefficients model that is actually fitted are specified in the reference below.
}
\value{
Object of class \code{rcmTest}.}
\references{
Van Wieringen, W.N., Berkhof, J., Van de Wiel, M.A. (2010), "A random coefficients model for regional co-expression associated with DNA copy number", \emph{Statistical Applications in Genetics and Molecular Biology}, Volume 9, Issue1, Article 25, 1-28. 

Van Wieringen, W.N., Van de Wiel, M.A., Van der Vaart, A.W. (2008), "A test for partial differential expression", \emph{Journal of the American Statistical Association}, 103(483), 1039-1049.
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\section{Warning}{
In case a covariate for the intercept is included in the design matrix \code{X} we strongly recommend the center, per feature, the data around zero.
}


\seealso{
\code{\link{RCMestimation}}, \code{\link{RCMrandom}}, \code{rcmTest}.
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

# test for significance of effect of X on Y
RCMtestResults <- RCMtest(Y, X, R, nBoot=2) 
summary(RCMtestResults) 
}
