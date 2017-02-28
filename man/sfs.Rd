\name{SFS}
\alias{sfs}
\title{Similarity-First Search multisweep algorithm}
\description{
  Return a ranking of the objects such that similar objects are ordered close to each other. If the matrix is \emph{Robinsonian}, then the ranking returned is a \emph{Robinson ordering}. 
}
\usage{
sfs(matrix, sfs_epsilon = 0, dissimilarity = FALSE, Robinsonian = FALSE, num_sweeps = 4)
}

\arguments{
  \item{matrix}{a 3-columns \code{data frame} with no symmetric entries, representing the list of all similarities (or dissimilarities) \eqn{(i, j, A_{ij})} between the pairs of objects to reorder.}
  \item{sfs_epsilon}{a numeric value that determines the threshold for entries of being considered the same.}
  \item{dissimilarity}{a boolean value equal to \code{TRUE} is the input data is a dissimilarity.}
  \item{Robinsonian}{a boolean value equal to \code{TRUE} is one wants to recognize a Robinsonian matrix}
  \item{num_sweeps}{an integer value that determines how many iterations of SFS shall be repeated.}
}
\value{
  Return a (row) vector of integers representing the rank of the objects, which is  0-based or 1-based accordingly with the input \code{matrix}.
}
\details{
	Given a a 3-columns \code{data frame} \eqn{(i, j, A_{ij})} listing all the similarities (or dissimilarities) among the objects, this function builds a \code{spMat} object in \emph{Armadillo} and computes a finite number of repeated SFS iterations (each called \emph{sweep}). The user may decide the threshold for which two entries are considered equivalent, meaning that if \eqn{|A_{ij} - A_{ik}| \leq} \code{sfs_epsilon}, then objects \eqn{j} and \eqn{k} have the same similarity (or dissimilarity) with respect to obejct \eqn{i}. By default, this threshold is set to \eqn{0}. 
	If not specified, the \code{matrix} represents the similarity information between objects. If \code{dissimilarity = TRUE}, then the SFS algorithm is modified by sorting the neighborhood of a visited vertex for decreasing values (instead than for increasing values).	
	The parameter \eqn{k=}\code{num_sweeps} sets the number of sweeps performed by \code{SFS()}, which directly affects the complexity of the function, as \code{SFS()} runs in \eqn{(k(n+m\log n))} time. By default, \code{num_sweeps=4}, as experiments show that for most of data fours sweeps are enough to find a good ranking.
If \code{Robinsonian = TRUE}, then \code{sfs()} also checks if the returend permutation is a Robinson ordering.
In this case, the number of sweeps is automatically set to the number of objects to rank \eqn{- 1}. In any case, efficient measures are implemented in order to avoid time consuming loops between consecutive SFS iterations. Note that at the moment the Robinsonian recognition is implemented only for similarity matrices.
Finally, the object returned by \code{SFS()} is a vector of integers, where the entry at position \eqn{i} represents the \eqn{i}-th object in the rank.
If the \code{matrix} is 0-based, then the returned vector is 0-based. If \code{matrix} is 1-based, then the returned vector is 1-based.
}
\examples{
	##install package in R
 	install.packages("SFS_0.1.tar.gz")
 	
 	##load package in R
 	library(SFS)
 
 	##invoke SFS on a R Matrix
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
  	A <- SFS::read(M)
  	SFS::sfs(A, Robinsonian = TRUE)
  	
  	##invoke SFS on a data-frame with 3-column data-frames with 0-based vertices, with (row, col, value) triples containing symmetric values
  	data <- c("0 0 1.0",
            "1 0 1.5",
            "2 0 1.9",
            "0 1 2.0",
            "1 1 2.5",
            "2 1 2.9",
            "0 2 3.0",
            "1 2 3.5",
            "2 2 3.9")
  	M <- read.table(textConnection(data))
  	A <- SFS::read(M, identical_val = TRUE)
  	SFS::sfs(A)

  	## invoke SFS on a \code{dist} from seriation package:
    library(seriation)
    data("iris");
    x <- as.matrix(iris[-5]);
    x <- x[sample(1:nrow(x)),];
    M <- dist(x)
    D <- SFS::read(M)
    SFS::sfs(D, sfs_epsilon = 0.01, dissimilarity = TRUE)

	## invoke sfs reading from file a Robinsonian matrix
	M <- read.table("examples/tmp_R_19_sim.txt")
	A <- SFS::read(M)
	SFS::sfs(A, Robinsonian = TRUE)
	
	## invoke sfs reading from file containing 3-columns (row, col, value) entries of an asymmetric non-Robinsonian similarity matrix with 1-based vertices
	M <- read.table("examples/list_130.txt")
	A <- SFS::read(M, symmetric = FALSE)
	SFS::sfs(A, num_sweeps = 7)
  
}
\author{Matteo Seminaroti (sfs) and Utz-Uwe Haus (R wrapping)}

