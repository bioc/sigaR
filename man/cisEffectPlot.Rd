\name{cisEffectPlot}
\alias{cisEffectPlot}

\title{ DNA-mRNA plot }

\description{
A variant on the boxplot, plotting the gene expression against the DNA copy number data. For each individual an open blue circle per call is plotted, all with their centerpoint
at the height of the individual's expression level. The radius of the circles is proportional to the corresponding call
probabilities. Call probabilities equal to zero reduce circles to dots. The red filled circles have a radius proportional to the
estimated expected call probabilities, with their centerpoints at the estimated mean expression for the respective call.
}

\usage{
cisEffectPlot(geneId, CNdata, GEdata, verbose=FALSE)
}

\arguments{
  \item{geneId}{ Indicator of the gene to be plotted. Indicator refers to the row in the \code{\link{ExpressionSet}}-object. }
  \item{CNdata}{ Object of class \code{cghCall}, containing (among others) annotion and call probabilities. Features should be matched with those of the accompanying \code{ExpressionSet}-object (as may be done using the \code{matchAnn2Ann}-function).  }
  \item{GEdata}{ Object of class \code{ExpressionSet}. Features should be matched with those of the accompanying \code{cghCall}-object (as may be done using the \code{matchAnn2Ann}-function). }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}

\note{
This function is a rewritten version of the \code{intCNGEan.plot} function of the \code{intCNGEan}-package.
}

\references{ 
Van Wieringen, W.N., Van de Wiel, M.A. (2009), "Non-parametric testing for DNA copy number induced differential mRNA gene expression", \emph{Biometrics}, 65(1), 19-29. 
}

\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{ \code{boxplot}, \code{cisEffectTune}, \code{cisEffectTest}, \code{matchAnn2Ann}}
\examples{
# load data
data(pollackCN16)
data(pollackGE16)

# plot DNA copy number vs. gene expression.
cisEffectPlot(225, pollackCN16, pollackGE16)
}
