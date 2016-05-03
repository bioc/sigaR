\name{pathwayPlot}
\alias{pathwayPlot}
\title{
Plot of the pathyway topology as reconstructed from from DNA copy number and gene expression data (one-sample only).
}
\description{
Plotting the topology of a pathway's regulatory network as reconstructed from from DNA copy number and gene expression data by the \code{pathway1sample}-function.
}
\usage{
pathwayPlot(pFit, directed = TRUE, tWidth = 1, cWidth = 1, geWidth=10, 
	cnWidth=10, circleDist = 1.5, gNames = NULL, main = "", remove = FALSE)
}

\arguments{
  \item{pFit}{ Object of class \code{pathwayFit} as returned by the function \code{pathway1sample}. }
  \item{directed}{ A \code{logical} indicating whether to plot directed or undirected \emph{trans}-effects. }
  \item{tWidth}{ A \code{numeric} that scales the width of the \emph{trans}-effect edges. }
  \item{cWidth}{ A \code{numeric} that scales the width of the \emph{cis}-effect edges. }
  \item{geWidth}{ A \code{numeric} that scales the width of the gene expression nodes. }
  \item{cnWidth}{ A \code{numeric} that scales the width of the DNA copy number nodes. }
  \item{circleDist}{ A \code{numeric} that scales the inner and outer circle. }
  \item{gNames}{ A \code{character} containing the gene names to written inside the nodes. }
  \item{main}{ The \code{character} to be plotted as plot title. }
  \item{remove}{ A \code{logical} indicating whether to remove genes without \emph{trans}-effects. }
}

\references{
Van Wieringen, W.N., Van de Wiel, M.A. (2012), "Modeling the \emph{cis}- and \emph{trans}-effect of DNA copy number aberrations on gene expression levels in a pathway", \emph{submitted for publication}. 
}

\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }

\seealso{
See also \code{pathway1sample}.
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
pFit <- pathway1sample(Y, X, lambda1=500)

# plot pathway topology
pathwayPlot(pFit, tWidth=5, cWidth=5)
}
