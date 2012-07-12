\name{hdMI}
\alias{hdMI}
\title{
Mutual information estimation.
}
\description{
The mutual information between two high-dimensional mutivariate random variables is estimated 
from two (high-dimensional matrix) under a normality or k-NN distributional assumption.
}
\usage{
hdMI(Y, X, method = "normal", k = 1, center = TRUE, rescale = TRUE)
}
\arguments{
  \item{Y}{
(High-dimensional) matrix. Rows are assumed to represent the samples, and columns represent the samples' genes or traits.
}
\item{X}{
(High-dimensional) matrix. Rows are assumed to represent the samples, and columns represent the samples' genes or traits. The number of rows of \code{X} must be identical to that of \code{Y}.
}
\item{method}{
Distributional assumption under which mutual information is to be estimated.
}
\item{k}{
k-nearest neighbor parameter.
}
\item{center}{
Logical indicator: should the columns (traits) of \code{Y} and \code{X} be centered at zero? Applied only under the normality assumption.
}
\item{rescale}{
Logical indicator: should \code{Y} and \code{X} be rescaled to have the same scale? Applied only under the k-NN assumption.
}
}
\value{
The mutual information estimate is returned as a \code{numeric}.
}
\references{
Van Wieringen, W.N., Van der Vaart, A.W. (2011), "Statistical analysis of the cancer cell's molecular entropy using high-throughput data", \emph{Bioinformatics}, 27(4), 556-563. 
}
\author{
Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} 
}
\seealso{
\code{\link{mutInfTest}}.
}
\examples{
data(pollackCN16)
data(pollackGE16)
hdMI(t(exprs(pollackGE16)), t(copynumber(pollackCN16)), method="knn")
}
