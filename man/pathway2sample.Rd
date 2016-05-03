\name{pathway2sample}
\alias{pathway2sample}
\title{
Penalized estimation of a pathyway's regulatory network from DNA copy number and gene expression data (two-sample).
}
\description{
The regulatory relationships between DNA copy number and gene expression within a pathway are modeled by a simulteneous-equations model.
Parameters of this model are fitted by minimizing of a penalized least squares criterion. The employed penalty is a combination of the lasso and the fused lasso. This combination
encourages within-sample sparsity (lasso), and limits the between-sample differences (fused lasso).
}
\usage{
pathway2sample(Y, X, id, lambda1 = 1, lambdaF = 1, method = "FL", 
	constr = TRUE, startCis = numeric(), startTrans1 =
        matrix(), startTrans2 = matrix(), epsilon = 0, verbose = FALSE)
}
\arguments{
   \item{Y}{ Object of class \code{matrix}. Rows are assumed to represent the samples, and columns represent the samples' gene expression levels. }
   \item{X}{ Object of class \code{matrix}. Rows are assumed to represent the samples, and columns represent the samples' genes or traits. The number of rows and columns of \code{X} must be identical to that of \code{Y}. }
   \item{id}{ An indicator variable of class \code{numeric} for the two groups to be compared. The groups should be coded as \code{0} and \code{1}. }
   \item{lambda1}{ Either a \code{numeric}- or \code{matrix}-object. The lasso parameter. In case \code{lambda1} is of class \code{numeric} its length is one, and the same penalty parameter is applied to all \emph{trans}-effects. 
	In case \code{lambda1} is of class \code{matrix} its column and row dimension equal the number of columns of \code{Y}. A possibly different penalty parameter is applied to each \emph{trans}-effect. 	}
   \item{lambdaF}{ Either a \code{numeric}- or \code{matrix}-object. The fused lasso parameter. In case \code{lambdaF} is of class \code{numeric} and of length one, the same penalty parameter is applied to all differential \emph{trans}-effects. 
	In case \code{lambdaF} is of class \code{matrix} its column and row dimension equal the number of columns of \code{Y}. A possibly different penalty parameter is applied to each differential \emph{trans}-effect. 	}
  \item{method}{ A \code{character}-object. Indicates which penalty to employ (see details). }
  \item{constr}{ \code{logical}. Should the \emph{cis}-effect (the direct effect of a column of \code{X} on column of \code{Y}) be positive? }
  \item{startCis}{ \code{numeric}. Starting values for the \emph{cis}-effect. }
  \item{startTrans1}{ \code{matrix}. Starting values for the \emph{trans}-effect of group 1 (coded as 0). }
  \item{startTrans2}{ \code{matrix}. Starting values for the \emph{trans}-effect of group 2 (coded as 1). }
  \item{epsilon}{ A \code{numeric}. Non-negative positive in the low-dimensional case. \code{epsilon} is to assume a positive value in the high-dimensional case. }
  \item{verbose}{ \code{logical}. Should intermediate output be printed on the screen? }
}
\details{
The model is fitted equation-by-equation. This is warranted by the assumption of independent errors. The expression levels of one gene is regressed on its own DNA copy number data and
the expression levels of all other genes in the pathway.

The \code{method}-option indicates which penalty is combined with the least squares loss function. 
In case \code{methode = FL}, this the fused lasso penalty (as described in Van Wieringen, W.N., Van de Wiel, M.A., 2012):
\deqn{\lambda_1 \| \Theta^{(a)} \|_1 + \lambda_1 \| \Theta^{(b)} \|_1 + \lambda_F \| \Theta^{(a)} - \Theta^{(b)} \|_1.  }
When \code{methode = FLs}, this penalty is simplified to:
\deqn{\lambda_1 \| \Theta^{(a)} + \Theta^{(b)} \|_1 + \lambda_F \| \Theta^{(a)} - \Theta^{(b)} \|_1.  }
The use of this penalty may be motivated as follows. 
The two samples used to share a common network architecture. 
One expects only a relatively limited number of edges to have changed.
Hence, the majority of edges will have the same sign, resulting in equality of the two penalties. 
An other motivation for this second penalty arises from the the observation that it is computationally faster. And, as
\deqn{
\lambda_1 \| \Theta^{(a)} \|_1 + \lambda_1 \| \Theta^{(b)} \|_1 \geq \lambda_1 \| \Theta^{(a)} + \Theta^{(b)} \|_1,  
}
it penalizes less. As such, the resulting \code{FLs} penalized estimates may be used as starting values for fitting the model with the \code{FL} penalty.
}
\value{
Object of class \code{pathwayFit}.
}
\references{
Van Wieringen, W.N., Van de Wiel, M.A. (2012), "Modeling the \emph{cis}- and \emph{trans}-effect of DNA copy number aberrations on gene expression levels in a pathway", \emph{submitted for publication}. 
}
\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{
See also \code{pathwayFit} and \code{pathway1sample}.
}
\examples{
# set number of genes (p) and samples (n)
p <- 10
n <- 1000

# sample cis-effects
beta <- abs(rnorm(p))

# sample trans-effects for first sample
Theta1 <- matrix(sample(c(-1,1), p^2, replace=TRUE, prob=c(0.2, 0.8)), ncol=p) * 
	matrix(runif(p^2), ncol=p) / 4
diag(Theta1) <- 1

# sample trans-effects for second sample
idDiff <- sample(which(Theta1 != 1), 10)
Theta2 <- Theta1
Theta2[idDiff] <- -Theta1[idDiff]

# sample error variances
Sigma <- diag(rchisq(p, df=1)/5 + 0.5)

# sample DNA copy number data of sample 1
X1 <- matrix(runif(n*p, min=-2, max=2), ncol=p)

# sample gene expression data
Y1 <- t(apply(X1, 1, function(Y, beta){ Y * beta }, beta=beta)) \%*\% t(solve(Theta1)) + 
	rmvnorm(n, sigma=solve(Theta1) \%*\% Sigma \%*\% t(solve(Theta1)))

# sample DNA copy number data of sample 1
X2 <- matrix(runif(n*p, min=-2, max=2), ncol=p)

# sample gene expression data
Y2 <- t(apply(X2, 1, function(Y, beta){ Y * beta }, beta=beta)) \%*\% t(solve(Theta2)) + 
	rmvnorm(n, sigma=solve(Theta2) \%*\% Sigma \%*\% t(solve(Theta2)))

# construct id-vector
id <- c(rep(0, n), rep(1, n))

# fit model
pFit <- pathway2sample(Y=rbind(Y1, Y2), X=rbind(X1, X2), id=id, lambda1=0, lambdaF=0.01)

# compare fit to "truth" for cis-effects
plot(pFit@Cis ~ beta, pch=20)

# compare fit to "truth" for differential trans-effects
penFits1 <- c(pFit@Trans1[upper.tri(Theta1)], pFit@Trans1[lower.tri(Theta1)])
penFits2 <- c(pFit@Trans2[upper.tri(Theta2)], pFit@Trans2[lower.tri(Theta2)])
truth1 <- c(Theta1[upper.tri(Theta1)], Theta1[lower.tri(Theta1)])
truth2 <- c(Theta2[upper.tri(Theta2)], Theta2[lower.tri(Theta2)])
plot(penFits1 - penFits2, truth1 - truth2, pch=20)
cor(penFits1 - penFits2, truth1 - truth2, m="s")
}
