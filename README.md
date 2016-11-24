# SFS
Similarity-First Search


-------------------------------------------------------------------------------

# Developer notes #

If you change any interface wrapping C++ code for R, regenerate
the wrapper code by 

	library(Rcpp)
	compileAttributes() 


and if you change the set of exported functions (or import from other packages)

	library(devtools)
	devtools::document("../SFS")

To prepare a CRAN upload, do

	library(devtools)
	devtools::check("../SFS")
