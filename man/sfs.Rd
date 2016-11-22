\name{Rcpp-SFS-Functions}
\alias{rcp_sfs}
\title{Functions in RcppSFS package}
\description{
  The Similarity-First Search package exports one function.
}
\usage{
rcpp_sfs(m)
}
\arguments{
  \item{m}{a matrix}
}
\value{
  \code{rcp_sfs()} returns a numeric vector specifying the 
  permutation of the input matrix \code{m} computed.
}
\details{
  The 
  These are example functions which should be largely
  self-explanatory. Their main benefit is to demonstrate how to write a
  function using the Armadillo C++ classes, and to have to such a
  function accessible from R.
}
\references{
  See the paper of Laurent and Seminaroti.
}
\examples{
  # install.packages("SFS_0.1.tar.gz")
  library(SFS)
  mtxt <- c("11 2 9 0 5 0 5 5 2 0 5 0 5 6 0 0 2 0 5",
         "2 11 2 0 9 0 8 5 10 0 5 0 5 2 0 0 10 0 8",
         "9 2 11 0 5 0 5 5 2 0 5 0 5 10 0 0 2 0 5",
         "0 0 0 11 0 3 0 0 0 3 0 3 0 0 10 3 0 9 0",
         "5 9 5 0 11 0 8 7 9 0 7 0 7 5 0 0 9 0 10",
         "0 0 0 3 0 11 0 0 0 10 0 6 0 0 5 8 0 5 0",
         "5 8 5 0 8 0 11 7 8 0 7 0 7 5 0 0 8 0 9",
         "5 5 5 0 7 0 7 11 6 0 10 0 8 7 0 0 6 0 7",
         "2 10 2 0 9 0 8 6 11 0 6 0 5 2 0 0 10 0 8",
         "0 0 0 3 0 10 0 0 0 11 0 6 0 0 4 9 0 5 0",
         "5 5 5 0 7 0 7 10 6 0 11 0 9 7 0 0 6 0 7",
         "0 0 0 3 0 6 0 0 0 6 0 11 0 0 9 6 0 10 0",
         "5 5 5 0 7 0 7 8 5 0 9 0 11 7 0 0 5 0 7",
         "6 2 10 0 5 0 5 7 2 0 7 0 7 11 0 0 2 0 5",
         "0 0 0 10 0 5 0 0 0 4 0 9 0 0 11 4 0 10 0",
         "0 0 0 3 0 8 0 0 0 9 0 6 0 0 4 11 0 4 0",
         "2 10 2 0 9 0 8 6 10 0 6 0 5 2 0 0 11 0 8",
         "0 0 0 9 0 5 0 0 0 5 0 10 0 0 10 4 0 11 0",
         "5 8 5 0 10 0 9 7 8 0 7 0 7 5 0 0 8 0 11")
  M <- as.matrix(read.table(textConnection(mtxt)))
  SFS::sfs(M)
}
\author{Matteo Seminaroti (sfs) and Utz-Uwe Haus (R wrapping)}

