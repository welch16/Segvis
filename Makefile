# build package documentation
doc:
	R -e 'library(roxygen2);roxygenize()'

# compile cpp code with Rcpp
#compile:
#	R -e 'library(Rcpp);compileAttributes(".")'
