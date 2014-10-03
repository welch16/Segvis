# build package documentation
doc:
	R -e 'roxygen2::roxygenize()'

clean:
	rm -f *~
	rm -f */*~
	rm -f .*~

# knit the vignettes
vignettes/%.pdf:vignettes/%.Rnw
	cd vignettes;R CMD Sweave --engine=knitr::knitr --pdf vignette.Rnw;cd ..

