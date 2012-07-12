\name{pathwayFit-class}
\Rdversion{1.1}
\docType{class}
\alias{pathwayFit-class}

\title{Class "pathwayFit" for storing the results of the function pathway1sample or pathway2sample. }
\description{	
The class \code{pathwayFit} is the output of a call to \code{pathway1sample} and \code{pathway2sample}.
It stores results from fitting a simultaneous-equations model from DNA copy number and gene expression data.
}

\section{Slots}{
  \describe{
    \item{\code{Cis}:}{ Object of class \code{"numeric"}. Vector of estimated \emph{cis}-effect.}
    \item{\code{Trans}:}{ Object of class \code{"matrix"}. Matrix containing the \emph{trans}-effects (one-sample only). }
    \item{\code{Trans1}:}{ Object of class \code{"matrix"}. Matrix containing the \emph{trans}-effects of the first sample (two-sample only). }
    \item{\code{Trans2}:}{Object of class \code{"matrix"}. Matrix containing the \emph{trans}-effects of the second sample (two-sample only). }
    \item{\code{Sigma}:}{Object of class \code{"numeric"}. Vector of estimated residual variances. }
    \item{\code{lambda1}:}{Object of class \code{"matrix"}. Lasso parameter(s) employed. }
    \item{\code{lambdaF}:}{Object of class \code{"matrix"}. Fused lasso parameter(s) employed. }
    \item{\code{constr}:}{Object of class \code{"logical"}. Indicator for parameter constraints on \emph{cis}-effect. }
    \item{\code{epsilon}:}{Object of class \code{"numeric"}. Constant used for the stabilization of estimation in a high-dimensional context. }
    \item{\code{method}:}{Object of class \code{"character"}. Indicator for method used in model fitting. }
  }
}
\section{Methods}{
  \describe{
    \item{pathwayPlot}{\code{signature(object = "pathwayFit")}: Plots the \code{pathwayFit}-object. }
 }
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
See also \code{pathway1sample} and \code{pathway2sample}.
}
\examples{
showClass("pathwayFit")
}
\keyword{classes}
