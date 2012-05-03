\name{cisEffectTable}
\alias{cisEffectTable}

\title{ Table of cis-effect test results }

\description{
Function to display the results of \code{\link{cisEffectTest}}-function in a table-format. 
Table may be restricted to a specified number of genes and sorted by relevant statistics.
}

\usage{
cisEffectTable(testRes, number=10, sort.by=NULL)
}

\arguments{
  \item{testRes}{ Object of class \code{cisTest} as produced by the \code{\link{cisEffectTest}}-function. }
  \item{number}{ Number of genes whose results are to be included in the table. }
  \item{sort.by}{ \code{character} indicating how the table is to sorted: \code{NULL} no sorting (genomic order), \code{p.value}, \code{R2} or \code{effect} sort the table by the corresponding statistic.  } 
}

\references{ 
Van Wieringen, W.N., Van de Wiel, M.A. (2009), "Non-parametric testing for DNA copy number induced differential mRNA gene expression", \emph{Biometrics}, 65(1), 19-29. 
}

\author{ Wessel N. van Wieringen: \email{w.vanwieringen@vumc.nl} }
\seealso{ \code{cisEffectTest} }
\examples{
# load data
data(pollackCN16)
data(pollackGE16)

# test cis-effect of DNA copy number on gene expression levels
cisRes <- cisEffectTest(pollackCN16, pollackGE16, 1:nrow(pollackGE16), 1, nPerm=25)

# display top results
cisEffectTable(cisRes, number=10, sort.by="R2")
}
