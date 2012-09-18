\name{matchAnn2Ann}
\alias{matchAnn2Ann}
\title{
Genomic location matching of two sets of features
}
\description{
Genomic location matching of two sets of features
}
\usage{
matchAnn2Ann(chr1, bpstart1, bpend1, chr2, bpstart2, bpend2, method = "distance", maxDist = 10000, minPerc = 0, reference = 1, ncpus = 1, verbose=TRUE)
}
\arguments{
  \item{chr1}{ Object of class \code{numeric} containing chromosome information of features from set 1. }
  \item{bpstart1}{ Object of class \code{numeric} containing start base pair information of features from set 1. Of same length as \code{chr1}. }
  \item{bpend1}{ Object of class \code{numeric} containing end base pair information of features from set 1. Of same length as \code{chr1}. }
  \item{chr2}{ Object of class \code{numeric} containing chromosome information of features from set 2. }
  \item{bpstart2}{ Object of class \code{numeric} containing start base pair information of features from set 2. Of same length as \code{chr2}. }
  \item{bpend2}{ Object of class \code{numeric} containing end base pair information of features from set 2. Of same length as \code{chr2}. }
  \item{method}{ Matching method to be applied, either \code{"distance"} or \code{"overlap"}. See below for details. }
  \item{maxDist}{ Maximum number of bases two features are allowed to be separated for a match. Only used in combination with \code{method="distance"}. }
  \item{minPerc}{ Minimum percentage of overlap between two features required for a match. Only used in combination with \code{method="overlap"}. }
  \item{reference}{ Platform that is taken as a reference in the calculation of the percentage, should equal 1 or two, referring to the platform. }
  \item{ncpus}{ Number of cpus to be used in the computation. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\details{
The features of set 1 (\code{chr1}, \code{bpstart1}, \code{bpend1}) are matched to the features of set 2 (\code{chr2}, \code{bpstart2}, \code{bpend2}).
That is, for every feature in set 2, features in set 1 are sought.

In case \code{method="distance"}, the midpoint of set 1 and set 2 features are calculated and for each feature of set 2
all features of set 1 with midpoints not further than \code{maxDist} are selected. If there are no features in set 1 satisfying this criterion, the feature of set 2 that could not be matched is discarded.

If \code{method="overlap"}, each feature of set 1 is matched to the feature of set 2 on the basis of the percentage of overlap. 
All features of set 1 with a percentage exceeding \code{minPerc} are selected. 
In case no feature in set 1 had any overlap with the features from set 2, the feature of set 2 that could not be matched is discarded.
}
\value{
An object of class \code{list}. Each list item is a three-column matrix with the matched features information. The first column contains feature numbers of set 1 in the order as supplied.
The second column contains feature numbers of set 2 in the order as supplied. Each row thus has two entries. The first entry contains the feature number of set 1 
that has been matched to second entry, representing the feature number of set 2.
The third column contains either the percentage of overlap (\code{method="overlap"}) or the distance between the the midpoints of the two features (\code{method="distance"}).
}
\section{Warning}{
Base pair information of features from both sets should be on the same scale!

Features with incomplete annotation information are removed before matching. For clarity, they are not included in the object with matched features.
}
\references{	
Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2012), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{submitted for publication}. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{\code{matchCGHcall2ExpressionSet}}
\examples{
# load data
data(pollackCN16)
data(pollackGE16)

# extract genomic information from cghCall-object
chr1 <- fData(pollackCN16)[,1]
bpstart1 <- fData(pollackCN16)[,2]
bpend1 <- fData(pollackCN16)[,3]

# extract genomic information from ExpressionSet-object
chr2 <- fData(pollackGE16)[,1]
bpstart2 <- fData(pollackGE16)[,2]
bpend2 <- fData(pollackGE16)[,3]

# match features from both platforms
matchedFeatures <- matchAnn2Ann(chr1, bpstart1, bpend1, chr2, bpstart2, bpend2, method = "distance", maxDist = 10000) 
}
