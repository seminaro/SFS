\name{SFS-package}
\alias{SFS-package}
\alias{SFS}
\docType{package}
\title{
\packageTitle{SFS}
}
\description{
\packageDescription{SFS}
}
\details{

\packageDESCRIPTION{SFS}
\packageIndices{SFS}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{SFS}

Maintainer: \packageMaintainer{SFS}
}

\references{
  
  M. Seminaroti. (2016): Combinatorial algorithms for the seriation
  problem. \emph{PhD thesis}, 1--209.
  
  M. Laurent, M. Seminaroti. (2015): A Lex-BFS-based recognition
  algorithm for Robinsonian matrices. \emph{Algorithms and Complexity:
  9th International Conference, CIAC 2015, Paris, France, May 20-22,
  2015}. 325--338.

  M. Laurent, M. Seminaroti. (2014) The quadratic assignment problem is
  easy for Robinsonian matrices with Toeplitz
  structure. \emph{Operations Research Letters}.

}

\keyword{ package }
\keyword{ cluster }
\seealso{
  \code{\link[seriation]{seriate}}
%  \code{\link{seriate.matrix}}
%  \code{\link{seriate.dist}}
  \code{\link{dist}}
}
\examples{
 ## we assume you already have done the equivalent of
 # install.packages("SFS_0.1.tar.gz")
 ## or
 # install.packages("SFS")
 ## load it
 library(SFS)
 ## create a 5x5 matrix with successive entries
 s <- matrix(1:25,5)
 ## and a random one
 r <- matrix(runif(100,1,1000),10)
 ## symmetrize both
 s[lower.tri(s)] = t(s)[lower.tri(s)]
 r[lower.tri(r)] = t(r)[lower.tri(r)]
 ## s is actually a Robinson matrix by construction
 ## run SFS, returns the permutation
 ## (expected here: identity or reverse identity)
 SFS::sfs(s)
 ## and r might give some more interesting result
 SFS::sfs(r)

}
