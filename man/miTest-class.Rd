\name{miTest-class}
\Rdversion{1.1}
\docType{class}
\alias{miTest-class}


\title{Class "miTest" for storing the results of the function mutInfTest. }
\description{	
The class \code{miTest} is the output of a call to \code{\link{mutInfTest}}.
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
    \item{summary}{\code{signature(object = "miTest")}: Prints the test results. }
	 }
}

\author{Wessel van Wieringen: \email{w.vanwieringen@vumc.nl}}

\seealso{
\code{\link{mutInfTest}}
}
\examples{
showClass("miTest")
}
\keyword{classes}
