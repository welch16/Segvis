# build package documentation
doc:
	R -e 'roxygen2::roxygenize()'

clean:
	rm -f *~
	rm -f */*~
	rm -f .*~
	rm -f inst/rscripts/*~

# knit the vignettes
vignettes/%.pdf:vignettes/%.Rnw
	cd vignettes;R CMD Sweave --engine=knitr::knitr --pdf $(<F);cd ..

# create the poster
poster/%.pdf:poster/%.tex poster/figs/fig0.pdf poster/figs/fig3.pdf poster/figs/fig2.pdf
	cd poster; pdflatex $(<F);pdflatex $(<F);pdflatex $(<F);cd ..

# figs poster
poster/figs/%.pdf:poster/figs/codes/%.R
	cd poster/figs/codes; R CMD BATCH $(<F);cd ../../..
