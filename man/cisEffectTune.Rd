\name{cisEffectTune}
\alias{cisEffectTune}
\title{ Pre-test and tuning. }

\description{
Decides which test to perform: loss vs. no-loss (tumor surpressor) or no-gain vs gain (proto-onco). Followed by a tuning algorithm that enhances the overal power of the FDR procedure by excluding genes
with either unbalanced (many samples having a high call probability of, say, a loss) or imprecise  (many call probabilities close to 0.5) soft calls, which is likely to increase the probability of 
detection for genes with a more favorable call probability distribution. 
}

\usage{
cisEffectTune(CNdata, GEdata, testStatistic, nGenes=250, nPerm=250, minCallProbMass=0.10, verbose=TRUE)
}


\arguments{
  \item{CNdata}{ Object of class \code{cghCall}, containing (among others) annotion and call probabilities. Features should be matched with those of the accompanying \code{ExpressionSet}-object (as may be done using the \code{matchAnn2Ann}-function).  }
  \item{GEdata}{ Object of class \code{ExpressionSet}. Features should be matched with those of the accompanying \code{cghCall}-object (as may be done using the \code{matchAnn2Ann}-function). }
  \item{testStatistic}{ Test statistic to be used, either \code{"wcvm"} or \code{"wmw"}, the weighted Cramer-Von Mises and the weighted Mann-Whitney test statistic, respectively. }
  \item{nGenes}{ Number of genes used for tuning. }
  \item{nPerm}{ Number of permutation used for tuning. }
  \item{minCallProbMass}{ A number inbetween 0 and 1. Genes with a marginal call probabilities in one of the classes smaller than \code{minCallProbMass} are discarded from further analysis. Effectively, this ensures identifiability of copy number effect on expression.}
  \item{verbose}{ Boolean to suppress output, either \code{FALSE} and \code{TRUE}. }
}

\value{
A \code{numeric}-object with the genes selected for testing. Numbering corresponds to genes of the pre-tuned, but matched data set.
}

\note{
This function is a rewritten version of the \code{intCNGEan.tune} function of the \code{intCNGEan}-package.
}

\references{ 
Van Wieringen, W.N., Van de Wiel, M.A. (2009), "Non-parametric testing for DNA copy number induced differential mRNA gene expression", \emph{Biometrics}, 65(1), 19-29. 
}

\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{ \code{matchAnn2Ann}, \code{cisEffectTest} }
\examples{
# load data
data(pollackCN16)
data(pollackGE16)

# select genes that are likely to have a significant genomic cis-effect on expression levels
genes2test <- cisEffectTune(pollackCN16, pollackGE16, "wmw", nGenes=50, nPerm=50)
}
