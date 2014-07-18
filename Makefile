# build package documentation
doc:
	R -e 'library(roxygen2);roxygenize()'

compile:
	R -e 'library(Rcpp);compileAttributes(".")'
