\name{entropyTest}
\alias{entropyTest}
\title{
One-sided two-sample test for entropy comparison
}
\description{
A one-sided two-sample test compares the entropy of a (high-dimensional) multivariate random variable between two groups.
The test is one-sided: one group is a priori suspected to have a larger entropy.
The null distribution is obtained via an efficient permutation resampling algorithm.
}
\usage{
entropyTest(Y, id, nPerm = 1000, method = "normal", k0 = 1, k1 = 1, center = TRUE, lowCiThres=0.10, ncpus=1, verbose=FALSE)
}
\arguments{
\item{Y}{
(High-dimensional) matrix. Columns are assumed to represent the samples, and rows represent the samples' genes or traits.
}
\item{id}{
An indicator variable for the two groups to be compared. The groups should be coded as \code{0} and \code{1}. 
There is an asymmetric interest in the groups: the group indicated by \code{1} is believed to exhibit a larger entropy.
}
\item{nPerm}{
Number of permutations.
}
\item{method}{
Distributional assumption under which entropy is to be estimated.
}
\item{k0}{
k-nearest neighbor parameter for group comprising of samples indicated by a zero in the indicator variable \code{id}.
}
\item{k1}{
k-nearest neighbor parameter for group comprising of samples indicated by a one in the indicator variable \code{id}.
}
\item{center}{
Logical indicator: should the rows of Y be centered at zero?
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
Object of \code{entTest}-class.
}
\references{
Van Wieringen, W.N., Van der Vaart, A.W. (2011), "Statistical analysis of the cancer cell's molecular entropy using high-throughput data", \emph{Bioinformatics}, 27(4), 556-563. 

Van Wieringen, W.N., Van de Wiel, M.A., Van der Vaart, A.W. (2008), "A test for partial differential expression", \emph{Journal of the American Statistical Association}, 103(483), 1039-1049.
}
\author{
Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl}
}
\seealso{
\code{\link{hdEntropy}}
}
\examples{
# load data
data(pollackGE16) 
Y <- exprs(pollackGE16) 

# assign samples to groups
id <- sample(c(0,1), 41, replace=TRUE) 

# perform testing and print test results
testRes <- entropyTest(t(Y), id, nPerm = 5, method="knn") 
summary(testRes) 
}
