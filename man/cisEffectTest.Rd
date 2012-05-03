\name{cisEffectTest}
\alias{cisEffectTest}
\title{ Nonparametric testing for copy number induced differential gene expression. }

\description{
A nonparametric test for the detection of copy number induced differential gene expression. The test incorporates the uncertainty of the calling of genomic aberrations:
weighted version of well-known test statistics are used. An efficient permutation re-sampling procedure is used for p-value calculation. The test statistics may be "shrunken" to borrow information across neighboring genes that share the same copy number signature. 
}

\usage{
cisEffectTest(CNdata, GEdata, genes2test=NULL, GEchr, analysisType="univariate", testStatistic="wcvm", nPerm = 10000, lowCiThres = 0.1, verbose=TRUE)
}

\arguments{
  \item{CNdata}{ Object of class \code{cghCall}, containing (among others) annotion and call probabilities. Features should be matched with those of the
accompanying \code{ExpressionSet}-object (as may be done using the \code{matchAnn2Ann}-function).  }
  \item{GEdata}{ Object of class \code{ExpressionSet}. Features should be matched with those of the accompanying
\code{cghCall}-object (as may be done using the \code{matchAnn2Ann}-function). }
  \item{genes2test}{ \code{Numeric} indicator vector containing row number of genes for which the DNA copy number cis-effect should be tested. The function \code{\link{cisEffectTune}} yields an optimal selection. }
  \item{GEchr}{ Column in the slot \code{featureData} of the \code{\link{ExpressionSet}}-object \code{GEdata} specifying the chromosome information of the features. }
  \item{analysisType}{ Indicator to determine whether the test statistic should be "shrunken" within a region. Either \code{"univariate"} (no shrinkage) 
or \code{"regional"} (shrinkage).}
  \item{testStatistic}{ Test statistic to be used, either \code{"wcvm"} or \code{"wmw"}, the weighted Cramer-Von Mises and
the weighted Mann-Whitney test statistic, respectively. }
  \item{nPerm}{ Number of permutations used for the p-value calculation. }
  \item{lowCiThres}{ A value between 0 and 1. Determines speed of efficient p-value calculation. 
	Genes with a probability smaller than 0.001 of a p-value smaller than \code{eff.p.val.thres} are discarded from the permutation analysis and their p-value is set equal to 1. 
	Should be chosen in accordance with the FDR-threshold for significance. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}


\value{
Object of class \code{cisTest}. 
}

\note{
This function is a rewritten version of the \code{intCNGEan.test} function of the \code{intCNGEan}-package.
}


\references{ 
Van Wieringen, W.N., Van de Wiel, M.A. (2009), "Non-parametric testing for DNA copy number induced differential mRNA gene expression", \emph{Biometrics}, 65(1), 19-29. 
}

\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{ \code{matchAnn2Ann}, \code{cisEffectTune}, \code{cisEffectTable}, \code{cisEffectPlot} }
\examples{
# load data
data(pollackCN16)
data(pollackGE16)

# test cis-effect of DNA copy number on gene expression levels
cisRes <- cisEffectTest(pollackCN16, pollackGE16, 1:nrow(pollackGE16), 1, nPerm=25)
}
