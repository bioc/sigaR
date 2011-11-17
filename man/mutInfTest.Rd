\name{mutInfTest}
\alias{mutInfTest}
\title{
Test for mutual information
}
\description{
A test evaluates the significance of the mutual information between two (high-dimensional) multivariate random variables.
The null distribution is obtained via an efficient permutation resampling algorithm.
}
\usage{
mutInfTest(Y, X, nPerm = 1000, method = "normal", k = 1, center = TRUE, rescale = TRUE, lowCiThres=0.10, ncpus=1, verbose=FALSE)
}
\arguments{
\item{Y}{
(High-dimensional) matrix. Columns are assumed to represent the samples, and rows represent the samples' genes or traits.
}
\item{X}{
(High-dimensional) matrix. Columns are assumed to represent the samples, and rows represent the samples' genes or traits. The number of columns of \code{X} must be identical to that of \code{Y}.
}
\item{nPerm}{
Number of permutations.
}
\item{method}{
Distributional assumption under which mutual information is to be estimated.
}
\item{k}{
k-nearest neighbor parameter.
}
\item{center}{
Logical indicator: should the rows of \code{Y} and \code{X} be centered at zero? Applied only under the normality assumption.
}
\item{rescale}{
Logical indicator: should \code{Y} and \code{X} be rescaled to have the same scale? Applied only under the k-NN assumption.
}
\item{lowCiThres}{
A value between 0 and 1. Determines speed of efficient p-value calculation. 
If the probability of a p-value being below \code{lowCiThres} is smaller than 0.001 (read: the test is unlikely to become significant), the permutation analysis is terminated and a p-value of 1.00 is reported. 
}
\item{ncpus}{Number of cpus used for the permutations.}
\item{verbose}{
Logical indicator: should intermediate output be printed on the screen?
}
}
\value{
Object of \code{miTest}-class.
}
\references{
Van Wieringen, W.N., Van der Vaart, A.W. (2011), "Statistical analysis of the cancer cell's molecular entropy using high-throughput data", \emph{Bioinformatics}, 27(4), 556-563. 

Van Wieringen, W.N., Van de Wiel, M.A., Van der Vaart, A.W. (2008), "A test for partial differential expression", \emph{Journal of the American Statistical Association}, 103(483), 1039-1049.
}
\author{
Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl}
}
\seealso{
\code{\link{hdMI}}
}
\examples{
# load data
data(pollackCN16) 
data(pollackGE16) 
Y <- t(exprs(pollackGE16)) 
X <- t(copynumber(pollackCN16)) 

# perform testing and print test results
testRes <- mutInfTest(Y, X, nPerm = 1000) 
summary(testRes) 
}
