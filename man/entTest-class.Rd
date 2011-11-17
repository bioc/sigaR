\name{entTest-class}
\Rdversion{1.1}
\docType{class}
\alias{entTest-class}


\title{Class "entTest" for storing the results of the function entropyTest. }


\description{	
The class \code{entTest} is the output of a call to \code{\link{entropyTest}}.
It stores results from a hypothesis test.
}


\section{Slots}{
  \describe{
    \item{\code{statistic}:}{Object of class \code{"numeric"}. Observed test statistic (i.e., estimated mutual information). }
    \item{\code{p.value}:}{Object of class \code{"numeric"}. P-value for the mutual information test. }
    \item{\code{null.dist}:}{Object of class \code{"numeric"}. The permutation null distribution for the test statistic.}
    \item{\code{nperm}:}{Object of class \code{"numeric"}. Number of permutation used for p-value calculation. }
    \item{\code{remark}:}{Object of class \code{"character"}. Tells whether the permutation algorithm was terminated prematurely or not. }
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "entTest")}: Prints the test results. }
	 }
}

\author{Wessel van Wieringen: \email{w.vanwieringen@vumc.nl}}

\seealso{
\code{entTest}}
\examples{
showClass("entTest")
}
\keyword{classes}
