\name{cisTest-class}
\Rdversion{1.1}
\docType{class}
\alias{cisTest-class}

\title{Class "cisTest" for storing the results of the function cisEffectTest. }
\description{	
The class \code{cisTest} is the output of a call to \code{\link{cisEffectTest}}.
It stores results from a hypothesis test.
}

\section{Slots}{
  \describe{
    \item{\code{geneInfo}:}{Object of class \code{"data.frame"}. E.g., annotation information of genes. }
    \item{\code{geneId}:}{Object of class \code{"numeric"}. Row number in \code{\link{ExpressionSet}}-object used in \code{{cisEffectTes}}, corresponding to a gene. }
    \item{\code{comparison}:}{Object of class \code{"numeric"}. Indicator of test performed, either \code{"1"} (loss vs. no-loss) or \code{"2"} (no-gain vs. gain). }
    \item{\code{av.prob1}:}{Object of class \code{"numeric"}. The estimated marginal call probability. }
    \item{\code{av.prob2}:}{Object of class \code{"numeric"}. The estimated marginal call probability. }
    \item{\code{effectSize}:}{Object of class \code{"numeric"}. Estimated genomic cis-effect on gene expression. }
    \item{\code{R2}:}{Object of class \code{"numeric"}. Percentage of explained variation in expression levels by . }
    \item{\code{regId}:}{Object of class \code{"numeric"}. Indicator for the region (\code{NULL} in the \code{regional}-analysis). }
    \item{\code{beginReg}:}{Object of class \code{"numeric"}. Row number in \code{ExpressionSet}-object corresponding to the first gene of the region (\code{NULL} in the \code{regional}-analysis). }
    \item{\code{endReg}:}{Object of class \code{"numeric"}. Row number in \code{ExpressionSet}-object corresponding to the last gene of the region (\code{NULL} in the \code{regional}-analysis). }
    \item{\code{shrinkage}:}{Object of class \code{"numeric"}. Amount of shrinkage applied in the regional analysis (\code{NULL} in the \code{regional}-analysis). }
    \item{\code{p.value}:}{Object of class \code{"numeric"}. P-value for the non-parametric test of the genomic cis-effect on expression levels. }
    \item{\code{adjP.value}:}{Object of class \code{"numeric"}. BH-multiple testing correct p-values. }
    \item{analysisType}{ Indicator whether the test statistic has been "shrunken" within a region. Either \code{"univariate"} (no shrinkage) or \code{"regional"} (shrinkage).}
    \item{testStatistic}{ Test statistic used, either \code{"wcvm"} or \code{"wmw"}, the weighted Cramer-Von Mises and the weighted Mann-Whitney test statistic, respectively. }
    \item{nPerm}{ Number of permutations used for the p-value calculation. }
  }
}
\section{Methods}{
  \describe{
    \item{cisEffectTable}{\code{signature(object = "cisTest")}: Prints the test results. }
	 }
}

\author{Wessel van Wieringen: \email{w.vanwieringen@vumc.nl}}

\seealso{
\code{\link{cisEffectTest}}}
\examples{
showClass("cisTest")
}
\keyword{classes}
