## SFS.R -- R-level functionality for SFS
## Currently only used to let roxygen (and devtools::check) generate the
## NAMESPACE file
#' @useDynLib SFS
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL
