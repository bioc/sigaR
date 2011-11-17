\name{rcmFit-class}
\Rdversion{1.1}
\docType{class}
\alias{rcmFit-class}

\title{Class "rcmFit" for storing the results of the function RCMestimation. }
\description{	
The class \code{rcmFit} is the output of a call to \code{\link{RCMestimation}}.
It stores results from fitting a random coefficients model.
}

\section{Slots}{
  \describe{
    \item{\code{betas}:}{ Object of class \code{"numeric"}. Vector of estimated global regression coefficients for each of the covariates in the design matrix.}
    \item{\code{tau2s}:}{ Object of class \code{"numeric"}. Vector of estimated regression coefficient variances for each of the covariates in the design matrix \code{X}.}
    \item{\code{sigma2s}:}{ Object of class \code{"numeric"}. Vector of estimated error variances for all genes. }
    \item{\code{rho}:}{Object of class \code{"numeric"}. Estimated correlation parameter between the error of two contiguous features. }
    \item{\code{av.sigma2s}:}{Object of class \code{"numeric"}. Average of the unshrunken estimated error variances. }
    \item{\code{shrinkage}:}{Object of class \code{"numeric"}. Applied shrinkage parameters in fitting the model. }
    \item{\code{loglik}:}{Object of class \code{"numeric"}. The log-likelihood of the fitted model. }
    \item{\code{corType}:}{Object of class \code{"character"}. Correlation structure of the error used. }
    \item{\code{X}:}{Object of class \code{"matrix"}. The design matrix. }
  }
}
\section{Methods}{
  \describe{
    \item{.RCMloss}{\code{signature(object = "rcmFit")}: Calculates the log-likelihood associated with the fitted model. }
    \item{RCMrandom}{\code{signature(object = "rcmFit")}: Samples from the distribution induced by the fitted model. }
    \item{summary}{\code{signature(object = "rcmFit")}: Prints the estimation result. }
 }
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
\code{\link{RCMestimation}}, \code{\link{RCMrandom}}.
}
\examples{
showClass("rcmFit")
}
\keyword{classes}
