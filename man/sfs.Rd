\name{SFS-Functions}
\alias{sfs}
\title{Functions in SFS package}
\description{
  The Similarity-First Search package exports one function.
}
\usage{
sfs(dissimilarity)
}
\arguments{
  \item{dissimilarity}{a matrix, a dist object, a sparse Matrix, or a
    3-column data frame}
}
\value{
  \code{sfs()} returns a (row) vector of 0-based integers specifying the 
  permutation computed to bring the dissimilarity matrix
  \code{dissimilarity} to Robinson form, or close to it if it is not
  Robinsonian.
}
\details{
  \code{sfs()} is the only function currently exported to R.
}
\references{
  M. Seminaroti. (2016): Combinatorial algorithms for the seriation
  problem. \emph{PhD thesis}, 1--209.
  
  M. Laurent, M. Seminaroti. (2015): A Lex-BFS-based recognition
  algorithm for Robinsonian matrices. \emph{Algorithms and Complexity:
  9th International Conference, CIAC 2015, Paris, France, May 20-22,
  2015}. 325--338.

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
  m <- as.matrix(read.table(textConnection(mtxt)))
  SFS::sfs(m)

  # if you have the seriation package at hand:
  library(seriation)
  data("iris");
  x <- as.matrix(iris[-5]);
  x <- x[sample(1:nrow(x)),];
  d <- dist(x)
  SFS::sfs(d)

  # you can also use 3-column data-frames, with (row, col, value) triples 
  data <- c("1 1 1.0",
            "2 1 1.5",
            "3 1 1.9",
            "1 2 2.0",
            "2 2 2.5",
            "3 2 2.9",
            "1 3 3.0",
            "2 3 3.5",
            "3 3 3.9")
  m <- read.table(textConnection(data))
  SFS::sfs(m)
}
\author{Matteo Seminaroti (sfs) and Utz-Uwe Haus (R wrapping)}

