\name{rcmTest-class}
\Rdversion{1.1}
\docType{class}
\alias{rcmTest-class}

\title{Class "rcmTest" for storing the results of the function RCMtest. }
\description{	
The class \code{rcmTest} is the output of a call to \code{\link{RCMtest}}.
It stores results from a hypothesis test.
}

\section{Slots}{
  \describe{
    \item{\code{statistic}:}{Object of class \code{"numeric"}. Observed test statistic (i.e., estimated mutual information). }
    \item{\code{p.value}:}{Object of class \code{"numeric"}. P-value for the mutual information test. }
    \item{\code{betas}:}{Object of class \code{"numeric"}. Vector of estimated global regression coefficients for each of the covariates in the design matrix. }
    \item{\code{tau2s}:}{Object of class \code{"numeric"}. Vector of estimated regression coefficient variances for each of the covariates in the design matrix. }
    \item{\code{sigma2s}:}{Object of class \code{"numeric"}. Vector of estimated error variances for all features. }
    \item{\code{rho}:}{Object of class \code{"numeric"}. Estimated correlation parameter between the error of two contiguous features. }
    \item{\code{av.sigma2s}:}{Object of class \code{"numeric"}. Average of the unshrunken estimated error variances. }
    \item{\code{shrinkage}:}{Object of class \code{"numeric"}. Type of shrinkage applied in the estimation. }
    \item{\code{loglik}:}{Object of class \code{"numeric"}. The log-likelihood of the fitted model. }
    \item{\code{nBoot}:}{Object of class \code{"numeric"}. Number of bootstraps used for p-value calculation. }
    \item{\code{corType}:}{Object of class \code{"character"}. Correlation structure used in the fitted model. }
    \item{\code{null.dist}:}{Object of class \code{"numeric"}. The permutation null distribution for the test statistic.}
    \item{\code{remark}:}{Object of class \code{"character"}. Tells whether the bootstrapping was terminated prematurely or not. }
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "rcmTest")}: Prints the test results. }
	 }
}

\author{Wessel van Wieringen: \email{w.vanwieringen@vumc.nl}}

\seealso{
\code{\link{RCMtest}}}
\examples{
showClass("rcmTest")
}
\keyword{classes}
