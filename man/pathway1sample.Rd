\name{pathway1sample}
\alias{pathway1sample}
\title{
Penalized estimation of a pathyway's regulatory network from DNA copy number and gene expression data (one-sample).
}
\description{
The regulatory relationships between DNA copy number and gene expression within a pathway are modeled by a simulteneous-equations model.
Parameters of this model are fitted by minimizing of a penalized least squares criterion. The employed penalty is that of the lasso, encouraging sparsity.
}
\usage{
pathway1sample(Y, X, lambda1 = 1, constr = TRUE, startCis=numeric(), 
	startTrans=matrix(), verbose = FALSE)
}
\arguments{
   \item{Y}{ \code{matrix}. Rows are assumed to represent the samples, and columns represent the samples' gene expression levels. }
   \item{X}{ \code{matrix}. Rows are assumed to represent the samples, and columns represent the samples' genes or traits. The number of rows and columns of \code{X} must be identical to that of \code{Y}. }
   \item{lambda1}{ \code{numeric} or \code{matrix}. The lasso parameter. In case \code{lambda1} is of class \code{numeric} and its length equals one, the same penalty parameter is applied to all \emph{trans}-effects. 
	In case \code{lambda1} is of class \code{matrix} its column and row dimension equal the number of columns of \code{Y}. A possibly different penalty parameter is applied to each \emph{trans}-effect. 	}
  \item{constr}{ \code{logical}. Should the \emph{cis}-effect (the direct effect of a column of \code{X} on column of \code{Y}) be positive? }
  \item{startCis}{ \code{numeric}. Starting values for the \emph{cis}-effect. }
  \item{startTrans}{ \code{matrix}. Starting values for the \emph{trans}-effect. }
  \item{verbose}{ \code{logical}. Should intermediate output be printed on the screen? }
}
\details{
The model is fitted equation-by-equation. This is warranted by the assumption of independent errors. The expression levels of one gene is regressed on its own DNA copy number data and
the expression levels of all other genes in the pathway.
}
\value{
Object of class \code{pathwayFit}.
}
\references{
Van Wieringen, W.N., Van de Wiel, M.A. (2012), "Modeling the \emph{cis}- and \emph{trans}-effect of DNA copy number aberrations on gene expression levels in a pathway", \emph{submitted for publication}. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
See also \code{pathwayFit} and \code{pathway2sample}.
}
\examples{
# set number of genes (p) and samples (n)
p <- 10
n <- 1000

# sample cis-effects
beta <- abs(rnorm(p))

# sample trans-effects
Theta <- matrix(sample(c(-1,1), p^2, replace=TRUE, prob=c(0.2, 0.8)), ncol=p) * 
		matrix(runif(p^2), ncol=p) / 4
diag(Theta) <- 1

# sample error variances
Sigma <- diag(rchisq(p, df=1)/5 + 0.5)

# sample DNA copy number data
X <- matrix(runif(n*p, min=-2, max=2), ncol=p)

# sample gene expression data
Y <- t(apply(X, 1, function(Y, beta){ Y * beta }, beta=beta)) \%*\% t(solve(Theta)) + 
	rmvnorm(n, sigma=solve(Theta) \%*\% Sigma \%*\% t(solve(Theta)))

# fit model
pFit <- pathway1sample(Y, X, lambda1=1, verbose=TRUE)

# compare fit to "truth" for cis-effects
plot(pFit@Cis ~ beta, pch=20)

# compare fit to "truth" for trans-effects
penFits <- c(pFit@Trans[upper.tri(Theta)], pFit@Trans[lower.tri(Theta)])
truth <- c(Theta[upper.tri(Theta)], Theta[lower.tri(Theta)])
plot(penFits ~ truth, pch=20)
}
