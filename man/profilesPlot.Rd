\name{profilesPlot}
\alias{profilesPlot}
\title{ CN-GE profiles plot }
\description{
Plots a sample's copy number and gene expression data side-by-side.
This visualizes the relation between CN and GE within an individual sample.
}
\usage{
profilesPlot(CNdata, GEdata, sampleNo, chr = 0, verbose=TRUE)
}
\arguments{
  \item{CNdata}{ Object of class \code{\link{cghCall}}, containing (among others) annotion and call probabilities. Features should be matched with those of the
accompanying \code{\link{ExpressionSet}}-object (as may be done using the \code{\link{matchCGHcall2ExpressionSet}}-function).  }
  \item{GEdata}{ Object of class \code{\link{ExpressionSet}}. Features should be matched with those of the accompanying
\code{cghCall}-object (as may be done using the \code{\link{matchCGHcall2ExpressionSet}}-function). }
  \item{sampleNo}{ Sample number of sample to be plotted. Corresponds to the order in which samples appear the \code{CNdata}- and \code{GEdata}-objects. }
  \item{chr}{ Chromosome number for which the profiles are to be plotted. Default \code{chr=0} for whole genome plotting. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\details{
The blue lines in the gene expression profile plot are the median expressions of genes that map to the same copy number segment.
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{ \code{\link{cghCall}}, \code{\link{ExpressionSet}}  }
\examples{
# load data
data(pollackCN16) 
data(pollackGE16) 

# plot CN and GE profiles alongside
profilesPlot(pollackCN16, pollackGE16, 23, 16) 
}
