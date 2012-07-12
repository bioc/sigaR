\name{hdEntropy}
\alias{hdEntropy}
\title{
Entropy estimation.
}
\description{
The (differential) entropy of a high-dimensional multivariate random variable is estimated from a (high-dimensional matrix) under a normality or k-NN distributional assumption.
}
\usage{
hdEntropy(Y, method = "normal", k = 1, center = TRUE, indKnn = TRUE)
}
\arguments{
\item{Y}{
(High-dimensional) matrix. Rows are assumed to represent the samples, and columns represent the samples' genes or traits.
}
\item{method}{
Distributional assumption under which entropy is to be estimated.
}
\item{k}{
k-nearest neighbor parameter.
}
\item{center}{
Logical indicator: should the columns of Y be centered around zero?
}
\item{indKnn}{
Logical indicator: should samples' individual contributions to the k-NN entropy be reported?
}
}
\value{
The entropy estimate is returned as a \code{numeric}.
}
\references{
Van Wieringen, W.N., Van der Vaart, A.W. (2011), "Statistical analysis of the cancer cell's molecular entropy using high-throughput data", \emph{Bioinformatics}, 27(4), 556-563. 
}
\author{
Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} 
}

\seealso{
\code{\link{entropyTest}}.
}
\examples{
data(pollackGE16)
hdEntropy(t(exprs(pollackGE16)), method="knn")
}
