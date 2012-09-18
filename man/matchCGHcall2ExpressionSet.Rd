\name{matchCGHcall2ExpressionSet}
\alias{matchCGHcall2ExpressionSet}
\title{ Genomic location matching of CN and GE data }
\description{
Integrative CN-GE analysis requires the copy number data of all genes on the expression array to be available.
\code{intCNGEan.match} matches the features of the copy number platform to the genes of the expression array.
This is done using their genomic locations on the basis of either proximity or overlap.
}
\usage{
matchCGHcall2ExpressionSet(CNdata, GEdata, CNchr, CNbpstart, CNbpend, GEchr, GEbpstart, GEbpend, method = "distance", reference=1, ncpus=1, verbose=TRUE)
}
\arguments{
  \item{CNdata}{ Object of class \code{cghCall}, containing (among others) annotion and call probabilities. }
  \item{GEdata}{ Object of class \code{ExpressionSet}. }
  \item{CNchr}{ Column in the slot \code{featureData} of the \code{\link{cghCall}}-object specifying the chromosome information of the features. }
  \item{CNbpstart}{ Column in the slot \code{featureData} of the \code{\link{cghCall}}-object specifying the start basepair information of the features. }
  \item{CNbpend}{ Column in the slot \code{featureData} of the \code{\link{cghCall}}-object specifying the end basepair information of the features. }
  \item{GEchr}{ Column in the slot \code{featureData} of the \code{\link{ExpressionSet}}-object specifying the chromosome information of the features. }
  \item{GEbpstart}{ Column in the slot \code{featureData} of the \code{\link{ExpressionSet}}-object specifying the start basepair information of the features. }
  \item{GEbpend}{ Column in the slot \code{featureData} of the \code{\link{ExpressionSet}}-object specifying the end basepair information of the features. }
  \item{method}{ Matching method to be applied, either \code{"distance"}, \code{"overlap"} or \code{"overlapPlus"}. See below for details. }
  \item{reference}{ Platform that is taken as a reference in the calculation of the percentage, should equal 1 or two, referring to the platform. }
  \item{ncpus}{ Number of cpus to be used in the computation. }
  \item{verbose}{ Logical indicator: should intermediate output be printed on the screen? }
}
\details{
Ideally full annotation information (chromosome number, start base pair, end base pair) for both copy number and gene expression data is available.
In case only start base pair information is available, let \code{CNbpend} and \code{GEbpend} refer to the same columns as \code{CNbpstart} and \code{GEbpstart}.
Base pair information of copy number and expression data should be on the same scale. 

Matching occurs on the basis of genomic locations. In case \code{method="distance"}, the midpoint of CN and GE features are calculated and for each gene on the
expression array the closest feature of the copy number platform is selected. If \code{method="overlap"}, each gene in the \code{ExpressionSet}-object
is matched to the feature from the copy number platform with the maximum percentage of overlap. If the maximum percentage of overlap equals zero, the gene is not included in the 
matched objects. If \code{method="overlapPlus"}, the features are first matched by their percentage of overlap (as with the \code{method="overlap"}-option). For all non-matched GE features
its closest two CN features (one down- and one upstream) are determined. If the copy number signature of these two CN features is identical, intrapolation seems reasonable, and
and the GE feature is matched to the closest of these two CN features. Hence, \code{method="overlapPlus"} makes use of the copy number data, consequently, matching may be different for different data sets. 
}
\value{
  A two-column matrix with the matched features entries.
The first column contains feature numbers of the \code{cghCall}-object.
The second column contains feature numbers of the \code{ExpressionSet}-object.
Each row thus has two entries. The first entry contains the feature number of the \code{cghCall}-object 
that has been matched to second entry, representing the feature number of the \code{ExpressionSet}-object.
}
\references{ 
Van Wieringen, W.N., Belien, J.A.M., Vosse, S.J., Achame, E.M., Ylstra, B. (2006), "ACE-it: a tool for genome-wide integration of gene dosage and RNA expression data", \emph{Bioinformatics}, 22(15), 1919-1920.

Van Wieringen, W.N., Unger, K., Leday, G.G.R., Krijgsman, O., De Menezes, R.X., Ylstra, B., Van de Wiel, M.A. (2012), "Matching of array CGH and gene expression microarray features for the purpose of integrative analysis", \emph{submitted for publication}. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\note{ 
The matching process implemented here is different from the one implemented in the (depreciated) \code{ACEit}-package (Van Wieringen et al., 2006).
}
\section{Warning}{
Features with incomplete annotation information are removed before matching. For clarity, they are not included in the objects with matched features.
}

\seealso{ \code{\link{cghCall}}, \code{\link{ExpressionSet}} }
\examples{
# load data
data(pollackCN16) 
data(pollackGE16) 

# match features from both platforms
featureMatch <- matchCGHcall2ExpressionSet(pollackCN16, pollackGE16, 1, 2, 3, 1, 2, 3) 
}
