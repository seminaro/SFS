\name{Read}
\alias{read}
\title{Read the input data}
\description{
  Read the similarity (or dissimilarity) information between the objects that want to be ordered and build a 3-columns \code{data frame} with no symmetric entries, where each row \eqn{(i, j, A_{ij})} represents the similarity \eqn{A_{ij}} between objects \eqn{i} and \eqn{j}.
}
\usage{
read(data, zero_epsilon = 1e-200, symmetric = TRUE, identical_val = FALSE)
}

\arguments{
  \item{data}{a representation of the (similarity or dissimilarity) between pairs of objects.}
  \item{zero_epsilon}{a numeric value that determines when values in \code{data} are to be considered \code{0}.}
  \item{symmetric}{a boolean value equal to \code{TRUE} if the input data is symmetric.}
  \item{identical_val}{a boolean value equal to \code{TRUE} if the data is a 3-columns \code{data frame} and symmetric edges are included.}
}
\value{
  Returns a 3-columns \code{data frame} representation of the original data listing all the (not repeated) pairwise (dis)similarities \eqn{(i, j, A_{ij})} between objects.
}
\details{
  The input data can be a weighted adjacency matrix (represented by objects \code{matrix} a \code{dist} and \code{data frame}) or a list of all the edges of the directed weighted graph (represented by a 3-col \code{data frame}) where each row \eqn{(i, j, A_{ij})} represents the similarity \eqn{A_{ij}} between objects \eqn{i} and \eqn{j}. If not specified, the data is assumed to be symmetric. Since by defualt the data is assumed to be symmetric, if it is represented by a 3-columns \code{data frame}, then it is assumed that symmetric edges are not listed, and thus by default \code{repeated_val = FALSE}. The reason for this choice is that for large symmetric data, it is more efficient to list the entries once. 
 However, note that if \code{symmetric = FALSE} then \code{identical_val = TRUE} automatically.
 }.  
 
\author{Matteo Seminaroti (sfs) and Utz-Uwe Haus (R wrapping)}

