\name{CNGEheatmaps}
\alias{CNGEheatmaps}
\title{ Parellel CN and GE heatmap plotting }
\description{
Heatmaps of DNA copy number and gene expression data are plotted together.
}
\usage{
CNGEheatmaps(CNdata, GEdata, location = "mode", colorbreaks = "equiquantiles")
}
\arguments{
  \item{CNdata}{ Object of class \code{cghCall}, containing (among others) annotion and call probabilities. Features should be matched with those of the
accompanying \code{ExpressionSet}-object (as may be done using the \code{matchCGHcall2ExpressionSet}-function).  }
  \item{GEdata}{ Object of class \code{ExpressionSet}. Features should be matched with those of the accompanying
\code{cghCall}-object (as may be done using the \code{matchCGHcall2ExpressionSet}-function). }
  \item{location}{ Parameter (\code{median}, \code{mean}, or \code{mode}) specifying how the center of the gene expression heatmap color-scheme is determined. }
  \item{colorbreaks}{ Parameter specifying how the color distribution of the gene expression heatmap is determined, either \code{equiquantiles} or \code{equidistant}. }
}
\details{
The DNA copy number data heatmap is generated as follows. The DNA copy number data are used to determine the genomic segments 
exhibiting no difference in DNA copy number between the array elements that map to that
segment. This resembles the dimension reduction technique employed in the \code{CGHregions}-package.
Consequently, within a segment the DNA copy number for one sample is constant, but may vary between samples. Note that a 
region may comprise of a whole chromosome, but also of a focal amplication. It is the DNA copy number signature of the segments that is 
depicted in the heatmap of the DNA copy number data. 

For the gene expression heatmap segments as constructed for the array CGH data are adopted. For each segment-sample combination
the expression levels of the genes that map to that segment are averaged. Next, the gene expression data is also collapsed 
to the segment format. It is this collapsed and averaged expression data that is depicted in the corresponding heatmap.
}
\references{ 
Van de Wiel, M.A., Van Wieringen, W.N. (2007), "CGHregions: dimension reduction for array CGH data with minimal 
information loss", \emph{Cancer Informatics}, 2, 55-63. 

Van Wieringen, W.N., Van de Wiel, M.A. (2009), "Non-parametric testing for DNA copy number 
induced differential mRNA gene expression", \emph{Biometrics}, 65(1), 19-29. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{ \code{cghCall}, \code{ExpressionSet}, \code{matchCGHcall2ExpressionSet}, \code{profilesPlot},  }
\examples{
# load data
data(pollackCN16) 
data(pollackGE16) 

# plot heatmaps
CNGEheatmaps(pollackCN16, pollackGE16, location = "mode", colorbreaks = "equiquantiles") 
}
