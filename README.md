# SFS
Similarity-First Search

Documentation is available for the R functions `SFS::read` and `SFS::sfs`.

-------------------------------------------------------------------------------

# Developer notes #

If you change any interface wrapping C++ code for R, regenerate
the wrapper code by 

```R
	R> library(Rcpp)
	R> compileAttributes() 
```

and if you change the set of exported functions (or import from other packages)

```R
	R> library(devtools)
	R> devtools::document("../SFS")
```

To prepare a CRAN upload, do

```R
	R> library(devtools)
	R> devtools::check("../SFS")
```

and also consider uploading to
[https://win-builder.r-project.org/upload.aspx](Win-Builder) to see if the code
works on Windows.
